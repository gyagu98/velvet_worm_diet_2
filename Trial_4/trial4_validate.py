#!/usr/bin/env python3
"""
Trial 4 - validation of prey COX1 hits (velvet worm diet, Peripatus Cuenca)

Takes the *_diet_clean.reads.tsv output of prey_cox1_screen_v2.py for every worm and:
  1. classifies each hit: prey / lab contaminant / implausible (marine) / non-animal / unresolved
  2. pulls the prey contig sequences from each sample's mmseqs2 contig FASTA
  3. de-duplicates overlapping contigs (same transcript rebuilt by TransLiG / Trans-ABySS / rnaSPAdes)
     with minimap2 all-vs-all (>=97% identity over >=50% of the shorter contig)
  4. assigns a confidence tier and the taxonomic rank that tier is allowed to support
       strong   : alignment >=150 bp and identity >=95%  -> reported to FAMILY
       standard : alignment >= 80 bp and identity >=90%  -> reported to ORDER only
  5. writes one FASTA of unique prey sequences for an independent BLAST check (NCBI nt)
  6. optional: --blast-remote runs that BLAST for you (needs blastn; slow) and flags
     anything whose best hit is a velvet worm (host) or disagrees with the MIDORI2 order

Needs: python3, minimap2 (already installed).

Example (run from the folder holding everything):
  python3 trial4_validate.py \
    --sample F734 F734_diet_clean.reads.tsv "F734_Peripatus_1_midsection_[whole_library_all]_mmseqs2_results.fasta" \
    --sample F739 F739_diet_clean.reads.tsv F739_*mmseqs2_results.fasta \
    --sample F744 F744_diet_clean.reads.tsv F744_*mmseqs2_results.fasta \
    -o Trial_4 -t 8
"""
import argparse, os, subprocess, sys, shutil
from collections import defaultdict, Counter

# ---------------------------------------------------------------- filtering rules
LAB_CONTAM_CLASSES = {"Amphibia", "Mammalia", "Reptilia", "Aves"}           # frog, human, mouse, snake...
LAB_CONTAM_ORDERS = {"Squamata", "Primates", "Rodentia", "Anura", "Caudata"}
IMPLAUSIBLE_CLASSES = {"Actinopteri", "Ostracoda", "Cephalopoda", "Bivalvia", "Chondrichthyes"}
IMPLAUSIBLE_ORDERS = {"Decapoda", "Halocyprida", "Perciformes", "Pleurobranchida",
                      "Centropomidae", "Calanoida", "Euphausiacea"}
ANIMAL_PHYLA = {"Arthropoda", "Mollusca", "Annelida", "Nematoda", "Chordata",
                "Onychophora", "Platyhelminthes", "Tardigrada", "Nemertea"}
HOST_WORDS = ("peripat", "onychophora", "epiperipatus", "oroperipatus", "euperipatoides")


def log(*a):
    print("[trial4]", *a, file=sys.stderr, flush=True)


def norm_id(x):
    return x.split()[0].lstrip(">").split("__")[0]


def parse_lineage(s):
    p = [x.strip() for x in s.split(";")] + [""] * 8
    g = lambda i: p[i] if p[i] else None
    return dict(phylum=g(1), cls=g(2), order=g(3), family=g(4), genus=g(5), species=g(6))


def classify(t):
    if not t["phylum"]:
        return "unresolved (no phylum)"
    if t["phylum"] not in ANIMAL_PHYLA:
        return "non-animal (fungus etc.)"
    if t["cls"] in LAB_CONTAM_CLASSES or t["order"] in LAB_CONTAM_ORDERS:
        return "lab contaminant (vertebrate)"
    if t["cls"] in IMPLAUSIBLE_CLASSES or t["order"] in IMPLAUSIBLE_ORDERS:
        return "implausible (marine)"
    if t["phylum"] == "Onychophora":
        return "host (velvet worm)"
    if not t["order"]:
        return "unresolved (no order)"
    return "prey"


def read_hits(path):
    rows = []
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        for line in f:
            v = dict(zip(hdr, line.rstrip("\n").split("\t")))
            if not v.get("read"):
                continue
            v["id"] = norm_id(v["read"])
            v["identity"] = float(v["identity"]); v["aln_len"] = int(float(v["aln_len"]))
            v.update(parse_lineage(v.get("lineage", "")))
            v["status"] = classify(v)
            rows.append(v)
    return rows


def fetch_seqs(fasta, wanted):
    out, name, buf = {}, None, []
    with open(fasta) as f:
        for line in f:
            if line.startswith(">"):
                if name in wanted:
                    out[name] = "".join(buf)
                name, buf = norm_id(line), []
            elif name in wanted:
                buf.append(line.strip())
        if name in wanted:
            out[name] = "".join(buf)
    return out


def tier(h):
    if h["aln_len"] >= 150 and h["identity"] >= 0.95:
        return "strong"
    if h["aln_len"] >= 80 and h["identity"] >= 0.90:
        return "standard"
    return "weak"


def cluster(seqs, threads, tmp):
    """union-find on minimap2 all-vs-all; returns {id: cluster_rep}"""
    ids = list(seqs)
    parent = {i: i for i in ids}
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    if len(ids) > 1 and shutil.which("minimap2"):
        fa = os.path.join(tmp, "prey_contigs.fa")
        with open(fa, "w") as o:
            for i in ids:
                o.write(f">{i}\n{seqs[i]}\n")
        res = subprocess.run(["minimap2", "-X", "-c", "-x", "asm20", "-t", str(threads), fa, fa],
                             capture_output=True, text=True)
        for line in res.stdout.splitlines():
            c = line.split("\t")
            q, ql, t, tl, match, blk = c[0], int(c[1]), c[5], int(c[6]), int(c[9]), int(c[10])
            if q == t or blk == 0:
                continue
            if match / blk >= 0.97 and blk / min(ql, tl) >= 0.5:
                parent[find(q)] = find(t)
    elif len(ids) > 1:
        log("WARNING: minimap2 not found - no de-duplication done")
    groups = defaultdict(list)
    for i in ids:
        groups[find(i)].append(i)
    rep = {}
    for members in groups.values():
        r = max(members, key=lambda m: len(seqs[m]))
        for m in members:
            rep[m] = r
    return rep


def remote_blast(fa, out):
    if not shutil.which("blastn"):
        log("blastn not installed - skipping remote BLAST (upload the FASTA to NCBI web BLAST instead)")
        return {}
    log("running remote BLAST against NCBI nt (can take 10-60 min)...")
    subprocess.run(["blastn", "-remote", "-db", "nt", "-query", fa, "-max_target_seqs", "5",
                    "-outfmt", "6 qseqid sseqid pident length evalue bitscore sscinames",
                    "-out", out], check=False)
    best = {}
    if os.path.exists(out):
        for line in open(out):
            c = line.rstrip("\n").split("\t")
            if c[0] not in best:
                best[c[0]] = c
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample", nargs=3, action="append", required=True,
                    metavar=("NAME", "READS_TSV", "CONTIG_FASTA"))
    ap.add_argument("-o", "--outdir", default="Trial_4")
    ap.add_argument("-t", "--threads", type=int, default=4)
    ap.add_argument("--blast-remote", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    all_rows, log_counts = [], defaultdict(Counter)
    fasta_out = open(os.path.join(a.outdir, "prey_unique_for_BLAST.fasta"), "w")

    for name, tsv, fasta in a.sample:
        log(f"--- {name}")
        hits = read_hits(tsv)
        for h in hits:
            log_counts[name][h["status"]] += 1
        prey = [h for h in hits if h["status"] == "prey"]
        log(f"hits: {len(hits)} | prey after filtering: {len(prey)}")
        seqs = fetch_seqs(fasta, {h["id"] for h in prey})
        missing = [h["id"] for h in prey if h["id"] not in seqs]
        if missing:
            log(f"WARNING: {len(missing)} prey contigs not found in FASTA (kept without de-dup)")
        for m in missing:
            seqs[m] = "N" * 100
        rep = cluster(seqs, a.threads, a.outdir)
        by_id = {h["id"]: h for h in prey}
        clusters = defaultdict(list)
        for h in prey:
            clusters[rep[h["id"]]].append(h)
        log_counts[name]["prey: duplicate contigs merged"] = len(prey) - len(clusters)
        for r, members in clusters.items():
            best = max(members, key=lambda m: (m["aln_len"], m["identity"]))
            tr = tier(best)
            orders = {m["order"] for m in members}
            fams = {m["family"] for m in members if m["family"]}
            row = dict(sample=name, seq_id=f"{name}|{r}", n_contigs_merged=len(members),
                       contig_len=len(seqs[r]), identity=best["identity"], aln_len=best["aln_len"],
                       tier=tr, phylum=best["phylum"], cls=best["cls"],
                       order=best["order"] if len(orders) == 1 else "conflict:" + "/".join(sorted(o or "?" for o in orders)),
                       family=(best["family"] if tr == "strong" and len(fams) == 1 else ""),
                       midori_lineage=best.get("lineage", ""), blast_top_hit="", blast_pident="",
                       blast_flag="",
                       contig_flag=("long contig (>5 kb): check it is not a host nuclear gene" if len(seqs[r]) > 5000 else ""))
            all_rows.append(row)
            if not seqs[r].startswith("NNNN"):
                fasta_out.write(f">{row['seq_id']} {best.get('lineage','')}\n{seqs[r]}\n")
        log_counts[name]["prey: unique sequences"] = len(clusters)
    fasta_out.close()

    if a.blast_remote:
        best = remote_blast(os.path.join(a.outdir, "prey_unique_for_BLAST.fasta"),
                            os.path.join(a.outdir, "blast_nt_results.tsv"))
        for row in all_rows:
            b = best.get(row["seq_id"])
            if not b:
                continue
            row["blast_top_hit"], row["blast_pident"] = b[6], b[2]
            if any(w in b[6].lower() for w in HOST_WORDS):
                row["blast_flag"] = "HOST - remove"

    cols = ["sample", "seq_id", "n_contigs_merged", "contig_len", "identity", "aln_len", "tier",
            "phylum", "cls", "order", "family", "midori_lineage", "blast_top_hit", "blast_pident", "blast_flag", "contig_flag"]
    with open(os.path.join(a.outdir, "all_samples.validated_prey.tsv"), "w") as o:
        o.write("\t".join(c if c != "cls" else "class" for c in cols) + "\n")
        for r in all_rows:
            o.write("\t".join(str(r[c]) for c in cols) + "\n")

    steps = sorted({k for c in log_counts.values() for k in c})
    names = [s[0] for s in a.sample]
    with open(os.path.join(a.outdir, "table_filtering_steps.tsv"), "w") as o:
        o.write("step\t" + "\t".join(names) + "\n")
        for s in steps:
            o.write(s + "\t" + "\t".join(str(log_counts[n][s]) for n in names) + "\n")

    print("\n=== Trial 4 summary (unique prey sequences after de-duplication) ===")
    for n in names:
        rs = [r for r in all_rows if r["sample"] == n]
        t = Counter(r["tier"] for r in rs)
        orders = Counter(r["order"] for r in rs)
        print(f"{n}: {len(rs)} unique prey sequences (strong {t['strong']}, standard {t['standard']}) | "
              + ", ".join(f"{k} {v}" for k, v in orders.most_common()))
    print(f"\nOutputs in {a.outdir}/ :")
    print("  all_samples.validated_prey.tsv   main table (one row per unique prey sequence)")
    print("  table_filtering_steps.tsv        how many hits were removed at each step (Methods/Supp table)")
    print("  prey_unique_for_BLAST.fasta      upload to https://blast.ncbi.nlm.nih.gov (blastn, nt)")


if __name__ == "__main__":
    main()
