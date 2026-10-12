# barcode-forge — plan

## Context
You need a small Python CLI that designs **N sets of inline/cell barcodes** for sequencing. All barcodes must be unique across sets, separated by a minimum **edit (Levenshtein) distance**. That distance is either set by the user or, by default, made as large as possible for the requested size. Barcodes must also follow wet-lab design rules: PCR-friendly GC and Tm, no long homopolymers, no hairpins or self-complementarity, and enough base diversity for the first sequencing cycles. The project is a new greenfield project at `/Users/hdash/code/barcode-forge`. It uses **uv** for env and packaging, **pytest** with test-driven development, and **git**.

**Prerequisite:** `uv` is not installed on this machine. Step 0 installs it with `brew install uv`, or with `curl -LsSf https://astral.sh/uv/install.sh | sh` if Homebrew is unavailable. Python 3.12.9 is already present.

## Framework choice
| Concern | Choice | Why |
|---|---|---|
| CLI | **Typer** | Options are defined from type hints, so it needs little boilerplate. It is built on Click, so it is mature and `typer.testing.CliRunner` makes CLI tests easy. It generates `--help` automatically. (argparse needs too much boilerplate; plain Click is fine but more verbose.) |
| Config | `dataclass` + optional `--config file.toml` (stdlib `tomllib`) | CLI flags override the TOML values, which override the preset. No pydantic dependency. |
| Edit distance | **rapidfuzz** (`rapidfuzz.distance.Levenshtein`, `Hamming`) | It is C-fast and supports `score_cutoff` for early exit. A plain Python DP reference implementation sits alongside it so tests can check the two agree. |
| Thermodynamics | **primer3-py** (`calc_tm`, `calc_hairpin`, `calc_homodimer`) | This is the standard tool for Tm and secondary-structure ΔG. Wheels exist for macOS. |
| Tests | pytest, pytest-cov, **hypothesis** | Property tests fit distance metrics well (symmetry, triangle inequality, d=0 ⇔ equal) and fit the guarantee checks on generator output. |
| Lint | ruff | Fast and needs one tool only. |
| Build | `uv init --package` (uv_build, `src/` layout, `[project.scripts] barcode-forge = "barcode_forge.cli:app"`) | |

## Design
**Algorithm: a global pool, then partitioning.** The generator builds *one* pool of `n_sets × set_size` barcodes in which **every pair** (within a set and across sets) has distance ≥ d. It then splits the pool into sets. This makes sets collision-free and also safe to pool together in one run. Steps:
1. Seeded random candidate generation → **biological filters** (cheap filters first, primer3 checks last).
2. Greedy acceptance: a candidate is accepted if `Levenshtein(c, b) ≥ d` for every accepted b (using `score_cutoff` to exit early).
3. **"Maximise" mode (default when `--min-distance` is not given):** search d from high to low (binary search over `[1, length]`), keeping the largest d at which the greedy pass reaches the target count within `--max-attempts`. The stats file reports the d it reached.
4. Partitioning: deal barcodes round-robin, shuffled by seed. Optionally, a balancing pass swaps barcodes so each set's per-position base composition (needed for low-diversity inline reads) is as close to 25% as possible.
5. If the target cannot be reached, the tool exits with a clear error that suggests a longer barcode or smaller d.

**Default biological rules for inline barcodes** (each one configurable or toggled off):
- GC content 40–60%
- max homopolymer run ≤ 2
- no dinucleotide repeats ≥ 3 (e.g. `ATATAT`)
- not its own reverse complement, and no palindromic stem ≥ 4 nt
- hairpin and homodimer ΔG > −3 kcal/mol (primer3)
- Tm within a window (default: computed, then within ±5 °C of the set median)
- exclude user-supplied motifs: restriction sites (`--avoid GAATTC,...`) and adapter substrings
- optional `--flank5/--flank3`: run the hairpin and homopolymer checks on the barcode *in context* (adapter + barcode + linker)
- distance metric is `levenshtein` by default; `hamming` and `seqlev` (Sequence-Levenshtein, Buschmann & Bystrykh 2013, which handles indels that shift into the downstream sequence) can be chosen with `--metric`. Seqlev is the most correct metric for inline barcodes; it comes in a later milestone.

**CLI (draft)**
```
barcode-forge generate -n 4 --set-size 96 --length 10 [--min-distance 3 | (default: maximise)]
    [--gc 0.4-0.6] [--max-homopolymer 2] [--metric levenshtein] [--avoid GAATTC]
    [--flank5 ...] [--flank3 ...] [--seed 42] [--config cfg.toml] -o outdir/
```
**Output:** `outdir/set_01.tsv … set_NN.tsv` are single-column files of barcodes with no header, one file per set. `outdir/stats.txt` records the version, full parameters, seed, and the reached min distance. It also gives min and mean pairwise distance within and across sets, GC and Tm ranges, the maximum homopolymer, worst hairpin ΔG, a per-position base composition table per set, and filter rejection counts.

## Layout: one concept per file, simple names, tests mirror the source
```
barcode-forge/
  pyproject.toml  .python-version  .gitignore  README.md
  src/barcode_forge/
    __init__.py
    cli.py                 # Typer app only: parse args, call core code, no logic
    config.py              # Config dataclass, defaults, TOML loading + merging
    sequence.py            # DNA alphabet, reverse_complement, gc_fraction, validation
    distance/
      __init__.py          # get_metric(name) -> callable
      hamming.py
      levenshtein.py
      seqlev.py
    rules/                 # one biological rule per file, each a pure `check(seq, cfg) -> bool`
      __init__.py          # Rule registry + RuleSet (runs rules cheapest first, counts rejections)
      gc.py
      homopolymer.py
      repeats.py           # dinucleotide/simple repeats
      palindrome.py        # self-reverse-complement and palindromic stems
      thermo.py            # primer3 Tm, hairpin, homodimer
      motifs.py            # avoided motifs (restriction sites, adapters)
      context.py           # rules evaluated with flank5/flank3
    design/
      __init__.py
      candidates.py        # seeded random candidate stream
      pool.py              # greedy distance-constrained selection
      search.py            # maximise min distance (binary search over d)
      partition.py         # split pool into sets
      balance.py           # per-position base-composition balancing
    report/
      __init__.py
      stats.py             # compute stats (pure data)
      writer.py            # write set_XX.tsv and stats.txt
  tests/                   # mirrors src: tests/distance/test_levenshtein.py, tests/rules/test_gc.py, ...
    conftest.py
```
Conventions: small pure functions, type hints everywhere, no I/O outside `report/writer.py` and `cli.py`, and the random seed is passed in explicitly (`random.Random(seed)`) instead of using global state.

## TDD milestones (red → green → refactor, one git commit per milestone)
0. Install uv; run `uv init --package barcode-forge`; `uv add typer rapidfuzz primer3-py`; `uv add --dev pytest pytest-cov hypothesis ruff`; `git init`; add `.gitignore`; make the initial commit.
1. **sequence.py + config.py:** test the DNA helpers and config defaults and merging first, then implement them.
2. **distance/:** start with tests for Hamming and Levenshtein against a reference DP and hypothesis properties, then implement. Later add seqlev.
3. **rules/:** write table-driven tests for each rule file (known good/bad sequences), then implement each rule. After that, test the RuleSet's ordering and rejection counting.
4. **design/candidates + pool + search:** tests cover the invariants for any seed: all unique, all pairwise ≥ d, all pass rules, determinism under one seed, a clean failure when infeasible, and search returning a d ≥ the d of any fixed-d run that succeeds.
5. **design/partition + balance:** tests check set sizes, no overlap between sets, that the union equals the pool, and that balancing improves per-position entropy.
6. **report/:** tests check file format (no header, one column), that the stats content parses, and that values match the recomputed values.
7. **cli.py:** CliRunner end-to-end tests cover flags, TOML config precedence, errors, and `--help`.
8. README with usage and the reasoning behind each design rule; `ruff check`; coverage ≥ 90%.

## Verification
- `uv run pytest --cov=barcode_forge` passes, with coverage ≥ 90%.
- `uv run barcode-forge generate -n 4 --set-size 96 --length 10 --seed 1 -o /tmp/out`: check that `set_01..04.tsv` each hold 96 lines with no header, and that `stats.txt` shows reached d ≥ 3 and cross-set min distance = within-set min distance.
- An independent check script recomputes all pairwise distances across every output file.
- Running the same seed twice gives byte-identical output.

## Out of scope for now (possible later)
- BK-tree or indexing for pools of 10k+ barcodes (O(N²) checks are fine up to a few thousand).
- Illumina i7/i5 presets and Nanopore presets.
- Pushing to a GitHub remote (this happens only if you ask).

## Changes made during implementation (after code and biology review)
- Default metric is now **seqlev** rather than levenshtein: the Levenshtein correction guarantee does not hold for inline reads. `design/index.py` makes it C-fast with prefix scans.
- The **Tm rule is opt-in**, because an embedded barcode never anneals alone. **GC bounds round to whole bases.**
- Added **sequencing chemistry**:
  - 2-colour by default, with a dark-start rule, poly-G as a reserved phantom barcode, and per-set channel-balance validation plus repair from spare candidates.
  - Optional reverse-complement separation (`--revcomp`).
- The context hairpin baseline uses the least-structured homopolymer spacer.
- Build backend is `uv_build`, not hatchling.
- Deferred: primer heterodimer checks; codon/frame rules for translated regions.
