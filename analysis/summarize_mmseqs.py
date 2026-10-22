#!/usr/bin/env python3
"""Summarize MMseqs2 12-column tabular similarity results without taxonomic overclaiming.

Expected input columns (no header required): query, target, fident, alnlen,
mismatch, gapopen, qstart, qend, tstart, tend, evalue, bits.
The `fident` column follows MMseqs2 easy-search convention (fraction 0..1),
but percentage values (0..100) are also accepted.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import statistics
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Optional, TextIO, Tuple

FIELDS = [
    "query", "target", "fident", "alnlen", "mismatch", "gapopen",
    "qstart", "qend", "tstart", "tend", "evalue", "bits",
]


def normalize_identity(value: float) -> float:
    """Return percent identity from MMseqs fraction or percentage formats."""
    if value < 0 or value > 100:
        raise ValueError(f"identity must be between 0 and 100, got {value}")
    return value * 100.0 if value <= 1.0 else value


def taxon_label(target: str) -> str:
    """Extract the free-text OS annotation; never treat it as a verified taxon."""
    match = re.search(r"(?:^|_)OS=(.*?)(?:_OX=|_GN=|_PE=|_SV=|$)", target)
    if not match:
        return "(no OS annotation)"
    return match.group(1).replace("_", " ").strip() or "(empty OS annotation)"


def parse_rows(handle: TextIO, delimiter: str = "\t") -> Tuple[Dict[str, dict], int, int]:
    """Read a 12-column MMseqs table and retain each query's strongest hit.

    Ties are resolved deterministically: higher bitscore, lower E-value, longer
    alignment, then lexical target ID. The top hit is a similarity annotation,
    not a dietary or species-level identification.
    """
    best: Dict[str, dict] = {}
    valid_rows = 0
    skipped_rows = 0
    for line_number, raw in enumerate(handle, 1):
        if not raw.strip() or raw.startswith("#"):
            continue
        parts = raw.rstrip("\r\n").split(delimiter)
        if len(parts) < 12:
            skipped_rows += 1
            continue
        try:
            fident = normalize_identity(float(parts[2]))
            alnlen = int(parts[3])
            evalue = float(parts[10])
            bits = float(parts[11])
            qstart, qend = int(parts[6]), int(parts[7])
        except (ValueError, IndexError):
            skipped_rows += 1
            continue
        if alnlen < 0 or not parts[0] or not parts[1]:
            skipped_rows += 1
            continue
        query = parts[0].lstrip(">")
        row = {
            "query": query,
            "target": parts[1],
            "taxon_label": taxon_label(parts[1]),
            "identity_percent": fident,
            "alignment_length": alnlen,
            "query_start": qstart,
            "query_end": qend,
            "evalue": evalue,
            "bitscore": bits,
            "query_coverage_percent": None,
        }
        # The query length is sometimes encoded in assembler headers as length_N.
        length_match = re.search(r"(?:^|[_:])length[_:](\d+)", query)
        if length_match:
            qlen = int(length_match.group(1))
            if qlen > 0:
                row["query_coverage_percent"] = round(100.0 * abs(qend - qstart + 1) / qlen, 3)
        rank = (bits, -evalue, alnlen, parts[1])
        previous = best.get(query)
        if previous is None or rank > previous["_rank"]:
            row["_rank"] = rank
            best[query] = row
        valid_rows += 1
    for row in best.values():
        row.pop("_rank", None)
    return best, valid_rows, skipped_rows


def summarize(input_path: Path, output_dir: Path, sample: str) -> dict:
    with input_path.open("r", encoding="utf-8", errors="replace") as handle:
        best, valid_rows, skipped_rows = parse_rows(handle)
    if not best:
        raise ValueError("No valid MMseqs rows found; expected at least 12 tab-delimited columns.")
    output_dir.mkdir(parents=True, exist_ok=True)
    ordered = sorted(best.values(), key=lambda row: row["query"])
    per_query_path = output_dir / f"{sample}.per_query.tsv"
    columns = ["query", "target", "taxon_label", "identity_percent", "alignment_length",
               "query_start", "query_end", "query_coverage_percent", "evalue", "bitscore"]
    with per_query_path.open("w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(ordered)

    taxon_counts = Counter(row["taxon_label"] for row in ordered)
    taxon_path = output_dir / f"{sample}.top_hit_annotations.tsv"
    with taxon_path.open("w", newline="", encoding="utf-8") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["similarity_annotation_from_target_header", "queries_with_top_hit"])
        writer.writerows(taxon_counts.most_common())

    summary = {
        "sample_label": sample,
        "input_file": input_path.name,
        "valid_alignment_rows": valid_rows,
        "skipped_rows": skipped_rows,
        "queries_with_hits": len(best),
        "distinct_top_hit_annotations": len(taxon_counts),
        "median_top_hit_identity_percent": round(statistics.median(r["identity_percent"] for r in ordered), 3),
        "median_top_hit_alignment_length_bp": round(statistics.median(r["alignment_length"] for r in ordered), 3),
        "interpretation": "Top-hit labels are similarity annotations only; they do not establish prey identity or consumption.",
        "outputs": [per_query_path.name, taxon_path.name],
    }
    with (output_dir / f"{sample}.run_summary.json").open("w", encoding="utf-8") as out:
        json.dump(summary, out, indent=2)
        out.write("\n")
    return summary


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="MMseqs2 12-column tabular output")
    parser.add_argument("--outdir", type=Path, default=Path("results/mmseqs"))
    parser.add_argument("--sample", default="sample", help="Short label used in generated output filenames")
    args = parser.parse_args(argv)
    try:
        summary = summarize(args.input, args.outdir, args.sample)
    except (OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
