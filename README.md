# barcode-forge

Design N sets of inline/cell DNA barcodes for sequencing. Every barcode is unique across all sets and separated by a minimum edit distance (or the maximum achievable). Barcodes pass PCR and sequencing design rules: GC content, homopolymers, repeats, palindromes, hairpin/dimer ΔG, Tm and avoided motifs.

## Development

```bash
uv sync
uv run pytest --cov
uv run ruff check .
uv run barcode-forge --help
```

Status: scaffold only. See [docs/PLAN.md](docs/PLAN.md) for the design and milestones.
