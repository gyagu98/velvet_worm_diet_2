#!/usr/bin/env python3
"""
prey_cox1_screen.py — fast diet barcoding from velvet worm (Onychophora) transcriptomes via prey COX1.

Strategy (quickest first, stops early when signal is clear):
  1. Pool candidate reads: mitochondrial reads + reads annotated as COX1/COI (+ optionally raw reads).
  2. Remove host COX1 reads by k-mer screen against your onychophoran COX1 (pure Python, fast).
  3. Map remaining reads to a COI reference DB (MIDORI2 / BOLD export) with minimap2
     (falls back to BLASTN if minimap2 is absent).
  4. Filter hits (identity, aligned length, MAPQ-agnostic best-hit with LCA over ties),
     drop Onychophora, summarise by taxon, and call PREY SIGNAL yes/no.

Requirements: Python >=3.8, minimap2 (recommended) or blastn on PATH. No Python packages needed.

Reference DB: FASTA whose headers contain taxonomy, e.g. MIDORI2:
  >MW123456.1.<1-1536>root_1;Eukaryota_2759;Arthropoda_6656;Insecta_50557;Diptera_7147;Muscidae_7367;Musca_7369;Musca_domestica_7370
  or BOLD-style:  >ID|Arthropoda,Insecta,Diptera,Muscidae,Musca,Musca domestica
Any header with ';' or ',' or '|' separated ranks works; last field = species.

Examples:
  # fastest: mito + annotated reads only
  python prey_cox1_screen.py --mito mito.fq.gz --annot-reads cox1_annotated.fa \
      --host-cox1 host_COX1.fa --ref MIDORI2_COI.fasta -t 16 -o diet

  # annotation is a table (read_id <tab> ... description) + raw reads to pull from
  python prey_cox1_screen.py --mito mito.fq.gz --annot-table reads_annot.tsv --raw R1.fq.gz R2.fq.gz \
      --host-cox1 host_COX1.fa --ref MIDORI2_COI.fasta -o diet

  # deep mode: also screen ALL raw reads (minimap2 is fast enough; ~minutes per 10M reads)
  python prey_cox1_screen.py --mito mito.fq.gz --raw R1.fq.gz R2.fq.gz --screen-raw ...
"""
import argparse, collections, gzip, os, re, shutil, subprocess, sys, tempfile

COX1_RE = re.compile(r"(?<![a-z0-9])(cox1|coxi|co1|coi|mt-co1|cytochrome[ _]c[ _]oxidase[ _]subunit[ _](1|i))(?![a-z0-9])", re.I)
HOST_TAXA = {"onychophora", "peripatidae", "peripatopsidae"}
# Co-extraction contamination: snakes were extracted alongside the worms.
CONTAM_TAXA = {"serpentes", "squamata", "reptilia", "ophidia"}


def log(*a):
    print("[prey-cox1]", *a, file=sys.stderr, flush=True)


def opn(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def read_seqs(path):
    """Yield (id, seq) from FASTA or FASTQ (gz ok)."""
    with opn(path) as f:
        first = f.read(1)
        if not first:
            return
        if first == ">":
            name, buf = None, []
            line = first + f.readline()
            while line:
                if line.startswith(">"):
                    if name is not None:
                        yield name, "".join(buf)
                    name, buf = line[1:].split()[0], []
                else:
                    buf.append(line.strip())
                line = f.readline()
            if name is not None:
                yield name, "".join(buf)
        else:
            hdr = first + f.readline()
            while hdr:
                s = f.readline().strip(); f.readline(); f.readline()
                yield hdr[1:].split()[0], s
                hdr = f.readline()


COMP = str.maketrans("ACGTN", "TGCAN")


def kmers(seq, k):
    seq = seq.upper()
    return {seq[i:i + k] for i in range(len(seq) - k + 1)}


def build_host_index(path, k):
    idx = set()
    for _, s in read_seqs(path):
        idx |= kmers(s, k)
        idx |= kmers(s.upper().translate(COMP)[::-1], k)
    log(f"host COX1 k-mers: {len(idx):,}")
    return idx


def is_host(seq, idx, k, frac):
    km = kmers(seq, k)
    if not km:
        return False
    return sum(1 for x in km if x in idx) / len(km) >= frac


def norm_id(x):
    # contig headers in the mmseqs2 fasta carry the annotation after "__"
    return x.split("__")[0]


def annotated_ids(table):
    ids = set()
    with opn(table) as f:
        for line in f:
            if COX1_RE.search(line):
                ids.add(norm_id(line.split("\t")[0].split()[0].lstrip(">@")))
    log(f"COX1-annotated read IDs from table: {len(ids):,}")
    return ids


def pool_reads(a, host_idx, out_fa):
    seen, n_in, n_host, n_short = set(), 0, 0, 0
    want_ids = annotated_ids(a.annot_table) if a.annot_table else set()

    def take(rid, s, out):
        nonlocal n_in, n_host, n_short
        if rid in seen:
            return
        seen.add(rid); n_in += 1
        if len(s) < a.min_len:
            n_short += 1; return
        if host_idx and is_host(s, host_idx, a.k, a.host_frac):
            n_host += 1; return
        out.write(f">{rid}\n{s}\n")

    with open(out_fa, "w") as out:
        for p in (a.mito or []) + (a.annot_reads or []):
            for rid, s in read_seqs(p):
                take(rid, s, out)
        if a.raw and (want_ids or a.screen_raw):
            for p in a.raw:
                tag = os.path.basename(p).split(".")[0]
                for rid, s in read_seqs(p):
                    rid = norm_id(rid)
                    if a.screen_raw or rid in want_ids:
                        take(f"{rid}/{tag}" if a.screen_raw else rid, s, out)
    kept = n_in - n_host - n_short
    if want_ids and n_in == 0:
        log("WARNING: no contig names matched the table; check headers")
    log(f"pooled reads: {n_in:,} | host COX1 removed: {n_host:,} | too short: {n_short:,} | kept: {kept:,}")
    return {"pooled": n_in, "host": n_host, "short": n_short, "kept": kept}


def parse_tax(header):
    h = header.split(None, 1)[-1] if " " in header and ";" not in header.split()[0] else header
    for sep in (";", ",", "|"):
        if h.count(sep) >= 2:
            parts = [re.sub(r"_\d+$", "", x).replace("_", " ").strip() for x in h.split(sep)]
            parts = [x for x in parts if x and not re.match(r"^(root|[A-Z]{1,3}\d+[\.\d<>\-]*.*)$", x)]
            return parts
    return [h]


def load_taxonomy(ref):
    tax = {}
    with opn(ref) as f:
        for line in f:
            if line.startswith(">"):
                h = line[1:].strip()
                tax[h.split()[0]] = parse_tax(h)
    log(f"reference sequences: {len(tax):,}")
    return tax


def run_minimap2(ref, reads, threads, out):
    cmd = ["minimap2", "-x", "sr", "-c", "--secondary=yes", "-N", "20", "-t", str(threads), ref, reads]
    log("running:", " ".join(cmd))
    with open(out, "w") as fh:
        subprocess.run(cmd, stdout=fh, stderr=subprocess.DEVNULL, check=True)
    hits = collections.defaultdict(list)
    for line in open(out):
        c = line.split("\t")
        alen, matches = int(c[10]), int(c[9])
        if alen:
            hits[c[0]].append((c[5], matches / alen, alen))
    return hits


def run_blast(ref, reads, threads, tmp):
    db = os.path.join(tmp, "refdb")
    subprocess.run(["makeblastdb", "-in", ref, "-dbtype", "nucl", "-out", db], check=True,
                   stdout=subprocess.DEVNULL)
    out = os.path.join(tmp, "hits.tsv")
    subprocess.run(["blastn", "-task", "megablast", "-query", reads, "-db", db, "-num_threads", str(threads),
                    "-max_target_seqs", "20", "-outfmt", "6 qseqid sseqid pident length", "-out", out], check=True)
    hits = collections.defaultdict(list)
    for line in open(out):
        q, s, pid, ln = line.rstrip().split("\t")
        hits[q].append((s, float(pid) / 100, int(ln)))
    return hits


def lca(lineages):
    out = []
    for ranks in zip(*lineages):
        if len(set(ranks)) == 1:
            out.append(ranks[0])
        else:
            break
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mito", nargs="+", help="mitochondrial reads (fa/fq, gz ok)")
    ap.add_argument("--annot-reads", nargs="+", help="reads already annotated as COX1 (fa/fq)")
    ap.add_argument("--annot-table", help="annotation table; rows mentioning COX1/COI select read IDs from --raw")
    ap.add_argument("--raw", nargs="+", help="raw reads (fa/fq, gz ok)")
    ap.add_argument("--screen-raw", action="store_true", help="map ALL raw reads, not only annotated IDs")
    ap.add_argument("--host-cox1", help="host (velvet worm) COX1 fasta for subtraction (strongly recommended)")
    ap.add_argument("--ref", required=True, help="COI reference fasta with taxonomy in headers (e.g. MIDORI2)")
    ap.add_argument("-o", "--out", default="prey_cox1", help="output prefix")
    ap.add_argument("-t", "--threads", type=int, default=os.cpu_count() or 4)
    ap.add_argument("--k", type=int, default=25, help="k-mer size for host subtraction")
    ap.add_argument("--host-frac", type=float, default=0.5, help="fraction of host k-mers to call a read host")
    ap.add_argument("--min-len", type=int, default=60)
    ap.add_argument("--min-id", type=float, default=0.90, help="min identity for a prey hit (0.97≈species, 0.90≈genus/family)")
    ap.add_argument("--min-aln", type=int, default=80)
    ap.add_argument("--min-reads", type=int, default=3, help="reads needed to report a taxon")
    ap.add_argument("--drop-taxa", nargs="*", default=set(), type=lambda s: s.lower(),
                    help="extra contaminant taxa to drop (e.g. --drop-taxa Homo Mus Drosophila)")
    a = ap.parse_args()

    if not (a.mito or a.annot_reads or a.raw):
        ap.error("give at least one of --mito / --annot-reads / --raw")
    if not a.host_cox1:
        log("WARNING: no --host-cox1; host reads will only be removed by taxonomy (Onychophora) after mapping")

    tmp = tempfile.mkdtemp(prefix="preycox1_")
    host_idx = build_host_index(a.host_cox1, a.k) if a.host_cox1 else None
    reads_fa = os.path.join(tmp, "candidates.fa")
    stats = pool_reads(a, host_idx, reads_fa)
    if stats["kept"] == 0:
        print("PREY SIGNAL: NO (no non-host candidate reads)"); return

    tax = load_taxonomy(a.ref)
    if shutil.which("minimap2"):
        hits = run_minimap2(a.ref, reads_fa, a.threads, os.path.join(tmp, "hits.paf"))
    elif shutil.which("blastn"):
        hits = run_blast(a.ref, reads_fa, a.threads, tmp)
    else:
        sys.exit("need minimap2 (recommended: conda install -c bioconda minimap2) or blastn on PATH")

    per_taxon = collections.Counter(); best_id = collections.defaultdict(list)
    n_host_tax = n_contam = n_lowq = 0
    with open(f"{a.out}.reads.tsv", "w") as rf:
        rf.write("read\tidentity\taln_len\tn_ties\tlineage\n")
        for q, hs in hits.items():
            hs = [h for h in hs if h[1] >= a.min_id and h[2] >= a.min_aln]
            if not hs:
                n_lowq += 1; continue
            top = max(h[1] for h in hs)
            ties = [h for h in hs if h[1] >= top - 0.005]
            lin = lca([tax.get(h[0], [h[0]]) for h in ties]) or ["unresolved"]
            if any(x.lower() in HOST_TAXA for x in lin):
                n_host_tax += 1; continue
            if any(x.lower() in CONTAM_TAXA or x.lower() in a.drop_taxa for x in lin):
                n_contam += 1; continue
            # cap resolution by identity
            depth = len(lin) if top >= 0.97 else max(1, len(lin) - 1) if top >= 0.93 else max(1, len(lin) - 2)
            name = ";".join(lin[:depth])
            per_taxon[name] += 1; best_id[name].append(top)
            rf.write(f"{q}\t{top:.3f}\t{max(h[2] for h in ties)}\t{len(ties)}\t{name}\n")

    reported = [(t, n) for t, n in per_taxon.most_common() if n >= a.min_reads]
    with open(f"{a.out}.taxa.tsv", "w") as tf:
        tf.write("lineage\treads\tmean_identity\tmax_identity\n")
        for t, n in reported:
            ids = best_id[t]
            tf.write(f"{t}\t{n}\t{sum(ids)/len(ids):.3f}\t{max(ids):.3f}\n")

    prey_reads = sum(n for _, n in reported)
    print("\n=== Prey COX1 screen ===")
    print(f"candidate reads pooled     : {stats['pooled']:,}")
    print(f"host removed (k-mer)       : {stats['host']:,}")
    print(f"host removed (taxonomy)    : {n_host_tax:,}")
    print(f"contaminant reads removed  : {n_contam:,} (snakes + --drop-taxa)")
    print(f"mapped, below thresholds   : {n_lowq:,}")
    print(f"prey-assigned reads        : {prey_reads:,} across {len(reported)} taxa (>= {a.min_reads} reads)")
    for t, n in reported[:15]:
        print(f"  {n:>7,}  {max(best_id[t]):.3f}  {t}")
    print(f"\nPREY SIGNAL: {'YES' if reported else 'NO'}")
    print(f"tables: {a.out}.taxa.tsv, {a.out}.reads.tsv")
    print("Note: check for lab/reagent contaminants (e.g. Homo, Drosophila, Mus) and other samples run on the same lane.")
    shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    main()
