# Reproducible MMseqs2 result summaries

This small, dependency-free utility turns MMseqs2 12-column tabular output into per-query top-hit annotations and a count table. It records identity, alignment length, bit score, E-value, and query coverage when the contig header encodes `length_N`.

It does **not** decide whether a sequence is prey, identify diet, resolve taxonomy, or assess contamination. Target-header organism names are unverified similarity annotations and must be interpreted with alignment coverage, database composition, controls, and independent evidence.

## Requirements

- Python 3.9+ (standard library only)
- MMseqs2 tabular output with columns: `query target fident alnlen mismatch gapopen qstart qend tstart tend evalue bits`

## Run

```bash
python analysis/summarize_mmseqs.py mmseqs_hits.tsv \
  --sample F734 \
  --outdir results/mmseqs
```

Outputs are `<sample>.per_query.tsv`, `<sample>.top_hit_annotations.tsv`, and `<sample>.run_summary.json`. The command keeps the strongest hit per query by bit score; ties are resolved deterministically by E-value, alignment length, and target ID. MMseqs2 fractional identity is converted to percent; percentage-form identity is also accepted.

For the uploaded F734 example, the original large inputs are not stored in GitHub. Keep assemblies, raw reads, full search tables, and reference databases in an appropriate data store; commit only small, reviewed summaries and code.

## Tests

```bash
python -m unittest discover -s analysis/tests -v
```

GitHub Actions runs these tests when analysis utility files change.
