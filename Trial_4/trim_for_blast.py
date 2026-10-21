#!/usr/bin/env python3
"""Trim each prey contig to just its COX1 region so NCBI BLAST can handle it.
Usage: python3 trim_for_blast.py results/prey_unique_for_BLAST.fasta MIDORI2_COI.fasta results/prey_COX1_trimmed_for_BLAST.fasta
"""
import sys, subprocess

inp, ref, out = sys.argv[1:4]
PAD = 50  # extra bases kept on each side of the COX1 match

seqs, name = {}, None
for line in open(inp):
    line = line.rstrip()
    if line.startswith(">"):
        name = line[1:]; seqs[name] = []
    elif name:
        seqs[name].append(line)
seqs = {k: "".join(v) for k, v in seqs.items()}
ids = {k.split()[0]: k for k in seqs}

print("Mapping prey contigs to MIDORI2 to find the COX1 region (a few minutes)...")
paf = subprocess.run(["minimap2", "-c", "-t", "8", "-k", "13", "-w", "5", ref, inp],
                     capture_output=True, text=True, check=True).stdout

best = {}
for l in paf.splitlines():
    f = l.split("\t")
    q, qs, qe, matches = f[0], int(f[2]), int(f[3]), int(f[9])
    if q not in best or matches > best[q][2]:
        best[q] = (qs, qe, matches)

n = 0
with open(out, "w") as o:
    for q, (qs, qe, _) in best.items():
        full = ids[q]; s = seqs[full]
        a, b = max(0, qs - PAD), min(len(s), qe + PAD)
        o.write(f">{q.split('|')[0]}|{q.split('|')[1] if '|' in q else q}_{a+1}-{b}\n{s[a:b]}\n")
        n += 1
print(f"Wrote {n} of {len(seqs)} sequences to {out} (each now ~100-700 bp instead of whole contigs).")
missing = len(seqs) - n
if missing:
    print(f"{missing} sequences had no COX1 match and were left out.")
