# barcode-forge

Design **N sets of inline / cell DNA barcodes** for short-read sequencing.

- All barcodes are unique across **all** sets, and every pair (within and across sets) is at least `d` edits apart, so sets can be pooled in one run.
- `d` is either given with `-d` or, by default, **maximised** for the requested number and length.
- Every barcode passes PCR and sequencing design rules (below), and each set is channel-balanced for the instrument chemistry.

## Install

```bash
uv sync
uv run barcode-forge --help
```

## Usage

```bash
# 4 sets of 96 ten-mers, maximal Sequence-Levenshtein distance, reproducible
uv run barcode-forge generate -n 4 -s 96 -l 10 --seed 1 -o out/

# fixed distance, avoid EcoRI/BamHI, check junctions with the adapter context
uv run barcode-forge generate -n 2 -s 48 -l 12 -d 4 \
    --avoid GAATTC --avoid GGATCC \
    --flank5 ACACTCTTTCCCTACACGACGCTCTTCCGATCT --flank3 ATCACG -o out/

# options from a TOML file; command-line flags override it
uv run barcode-forge generate --config design.toml -o out/
```

`design.toml` keys are the `Config` field names (see `src/barcode_forge/config.py`), plus `disable = [...]`:

```toml
n_sets = 4
set_size = 96
length = 10
avoid_motifs = ["GAATTC"]
chemistry = "2colour"     # "4colour" or "none"
disable = ["homodimer"]
```

**Output** (`-o out/`):
- `set_01.tsv … set_NN.tsv`: one barcode per line, single column, no header.
- `stats.txt`: every parameter (including the resolved seed, so any run can be reproduced), achieved distances within and across sets, GC/Tm/ΔG ranges, rule rejection counts, per-position base composition and channel balance per set.

## Design rules and why

| Rule | Default | Rationale |
|---|---|---|
| Distance metric | `seqlev` | Sequence-Levenshtein (Buschmann & Bystrykh 2013). An indel in an inline barcode shifts downstream bases into the read window; plain Levenshtein's correction guarantee does not hold there. `levenshtein` / `hamming` (fixed-position indexes only) are available. |
| GC content | 40–60%, rounded to whole bases | Even PCR amplification; at 8 nt the strict window would admit only 4/8 GC. |
| Homopolymer | ≤ 2 | Polymerase slippage and sequencing indel errors. |
| Dinucleotide repeat | ≤ 2 copies | Low-complexity slippage and mispriming. |
| Palindrome | ≤ 4 nt | Self-annealing; restriction-site-like. |
| Hairpin / homodimer ΔG | > −3 / > −6 kcal/mol | primer3 at 50 mM Na⁺, 1.5 mM Mg²⁺, 0.2 mM dNTP, 250 nM oligo. |
| Flank context | on if `--flank5/3` given | Junction homopolymers, motifs formed across junctions, barcode-induced hairpins with the adapter. |
| Avoided motifs | none | `--avoid`, checked on both strands. |
| Dark start | ≤ 1 leading G | 2-colour SBS: G is dark; GG at read start hurts cluster registration. |
| Poly-G distance | on (2-colour) | Failed/empty clusters read as poly-G, so every barcode is ≥ d from it. |
| Channel balance | per set, per position | 2-colour (NextSeq/NovaSeq): red (A,C) and green (A,T) ≥ 25%, dark (G) ≤ 50%. 4-colour: red (A,C), green (G,T) ≥ 25%. |
| Reverse complements | off (`--revcomp`) | For non-directional libraries, where a barcode may be read on either strand. |
| Tm window | off (`--tm-tolerance`) | An embedded barcode never anneals alone; enable when the barcode itself primes. |

Disable a rule with `--disable <rule>` (repeatable): `homopolymer`, `repeats`, `palindrome`, `hairpin`, `homodimer`, `tm`, `dark_start`.

**Error correction:** with minimum distance `d` (decoded with the same metric), up to `⌊(d−1)/2⌋` errors can be corrected and `d−1` detected. For barcodes that follow variable-length sequence (e.g. later split-pool rounds), anchor on the linker before decoding.

## How it works

1. Generate unique random candidates (or enumerate all of them for length ≤ 10). Screen them with the rules, cheapest first.
2. Greedily build **one global pool** where every pair is ≥ `d` apart. In maximise mode, binary-search the largest feasible `d`. Neighbour checks run in C (rapidfuzz), including a prefix-scan formulation of Sequence-Levenshtein.
3. Split the pool into sets, swap barcodes between sets to balance base composition, and repair any channel shortfall from spare candidates.
4. Fail loudly, with a suggestion, if the request is infeasible.

## Limitations

- Greedy construction gives a good, but not provably maximal, `d`.
- Not checked: heterodimers with your primers (use `--avoid` and flanks), and start/stop codons for barcodes in translated regions.
- primer3 rules require barcodes ≤ 60 nt; disable `hairpin`, `homodimer` and `tm` for longer barcodes.

## Development

See [AGENTS.md](AGENTS.md) for conventions and [docs/PLAN.md](docs/PLAN.md) for the original design.

```bash
uv run pytest -q --cov
uv run ruff check . && uv run ruff format --check .
```
