# AGENTS.md

barcode-forge is a CLI that designs N sets of **inline/cell DNA barcodes** for sequencing. All barcodes are unique across sets, every pair is ≥ d apart by edit distance (d is given or maximised), and each barcode passes PCR and sequencing design rules. Full design: `docs/PLAN.md`.

## Commands
```bash
uv sync                          # install
uv run pytest -q                 # tests (add --cov for coverage)
uv run ruff check . && uv run ruff format --check .
uv run barcode-forge --help
```

## Workflow
- **TDD.** Write a failing test, make it pass, then refactor. Never weaken a test to make it pass.
- **Atomic commits:** one concept per commit (e.g. one rule module + its tests). Use imperative subjects.
- Add dependencies only with `uv add` (or `uv add --dev`), and only when needed.
- Run the full test suite and ruff before every commit.

## Conventions
- **One concept per file**, simple names. `tests/` mirrors `src/barcode_forge/`.
- Use pure functions and type hints everywhere. I/O happens only in `cli.py` and `report/writer.py`.
- Randomness always goes through an explicit `random.Random` built from the seed. Never use global `random`. The same seed must give byte-identical output.
- `Config` (`config.py`) is the single source of parameters. It is frozen and validated, and setting a rule parameter to `None` disables that rule.
- Sequences are uppercase `ACGT` strings. Validate at the boundaries with `sequence.validate`.
- Keep comments short and explain the *biological why*, with a citation where a threshold comes from literature or vendor guidance.

## Module contracts
| Module | Contract |
|---|---|
| `sequence.py` | `DNA`, `validate`, `reverse_complement`, `gc_fraction` |
| `composition.py` | `positional_composition(seqs)`, `imbalance(seqs)` |
| `distance/<metric>.py` | `distance(a, b, *, score_cutoff: int \| None = None) -> int`. If the true distance > cutoff, it returns `cutoff + 1` (rapidfuzz semantics). |
| `distance/__init__.py` | `get_metric(name) -> Metric` for `hamming`, `levenshtein`, `seqlev` |
| `rules/<rule>.py` | Each module exports `RULE = Rule(name, cost, check)` where `check(seq, cfg) -> bool` and True means the sequence passes. Each rule also exports its own measurement function (e.g. `longest_homopolymer`) for reuse by stats. |
| `rules/ruleset.py` | `Rule` dataclass; `RuleSet(rules, cfg)` runs rules cheapest-first, `.passes(seq)`, `.rejections: Counter[str]` |
| `config.py` | Frozen `Config`; `from_layers(*dicts, disable=...)` merges defaults < TOML < CLI in one validation; `without(*rules)` |
| `design/` | `candidates` → `pool` (greedy, all pairs ≥ d, `reserved`, `revcomp`) → `search` (maximise d) → `partition` → `balance` → `repair` (channel balance from spares); `pipeline.run(cfg) -> Design` orchestrates them |
| `design/index.py` | `CloseIndex(metric, d)`: `add`, `has_close`. Picks a C path per metric (rapidfuzz, Hamming with `pad=False`, seqlev prefix scans). Add new metrics here. |
| `design/channels.py` | `channel_signal`, `channel_violations` per set for `2colour` / `4colour` |
| `report/` | `stats.compute(...)` is pure; `writer` writes `set_NN.tsv` (one column, no header) and `stats.txt` |

## Biological invariants (tests must enforce)
- Every barcode in the output, across **all** sets, is unique and pairwise ≥ the reported min distance under the chosen metric.
- Every barcode passes every enabled rule, including rules evaluated in flank context.
- Defaults are GC 40–60%, homopolymer ≤ 2, hairpin ΔG > −3 and homodimer ΔG > −6 kcal/mol (primer3; 50 mM Na⁺, 1.5 mM Mg²⁺, 0.2 mM dNTP), Tm window opt-in (`tm_tolerance`).
- The default metric is `seqlev`, because indels shift inline barcodes. Levenshtein is the alternative; Hamming is for fixed-position indexes only.
- On 2-colour chemistry (the default): no leading GG, every barcode ≥ d from poly-G, and every set channel-balanced at every position. Otherwise the design raises.
- Seqlev is not a metric (the triangle inequality fails), so never rely on metric-space indexing (e.g. BK-trees).
- If a design is infeasible, fail loudly with an actionable message. Never silently return fewer barcodes.
