#!/usr/bin/env python3
"""Trial 5: decide which COX1 sequences are real prey using the NCBI BLAST hit table.

Every sequence is judged by its best BLAST hit (highest bit score) in core_nt:
  * no hit                          -> dropped (unverifiable)
  * best hit is a velvet worm       -> HOST (Onychophora), dropped
  * best hit is human/mouse/frog/snake -> lab contamination, dropped
  * best hit is a marine-only group -> flagged implausible, dropped
  * anything else                   -> PREY, reported at order level
                                       (family too if aln >= 150 bp and id >= 95%)
Sequences hitting both a velvet worm and a prey taxon equally well count as host.

Usage:
  python3 trial5_blast_verify.py HIT_TABLE.txt [more hit tables...] -o trial5_results
Needs internet (looks up species names at NCBI). No extra packages.
"""
import sys, os, json, time, argparse, urllib.request, urllib.parse
from collections import defaultdict, Counter

CONTAM = {"Homo", "Mus", "Rattus", "Amphibia", "Serpentes", "Squamata", "Bos", "Sus", "Canis", "Felis", "Gallus"}
MARINE = {"Cnidaria", "Echinodermata", "Porifera", "Ctenophora", "Bryozoa", "Brachiopoda",
          "Cephalopoda", "Bivalvia", "Actinopteri", "Chondrichthyes", "Maxillopoda", "Thecostraca",
          "Copepoda", "Hexanauplia", "Polychaeta", "Ascidiacea", "Thaliacea"}
RANKS = ["phylum", "class", "order", "family", "genus", "species"]
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"

def get(url, tries=4):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return r.read().decode()
        except Exception:
            time.sleep(2 + i * 3)
    raise SystemExit("Could not reach NCBI - check your internet connection and rerun.")

def read_tables(paths):
    best, allhits, queries = {}, defaultdict(list), []
    for p in paths:
        for line in open(p):
            if line.startswith("# Query:"):
                q = line.split(":", 1)[1].strip()
                if q not in queries:
                    queries.append(q)
                continue
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            q, acc, pid, aln, bits = f[0], f[1], float(f[2]), int(f[3]), float(f[11])
            allhits[q].append((bits, acc, pid, aln))
            if q not in queries:
                queries.append(q)
    for q, hs in allhits.items():
        best[q] = max(hs)
    return queries, best, allhits

def lineages(accs):
    accs = sorted(set(accs)); acc2tax, out = {}, {}
    for i in range(0, len(accs), 150):
        chunk = accs[i:i + 150]
        d = json.loads(get(EUTILS + "esummary.fcgi?db=nuccore&retmode=json&id=" + ",".join(chunk)))["result"]
        for uid in d.get("uids", []):
            r = d[uid]; acc2tax[r.get("accessionversion", "")] = (str(r.get("taxid", "")), r.get("title", ""))
        time.sleep(0.4)
    taxids = sorted({t for t, _ in acc2tax.values() if t})
    tax = {}
    import xml.etree.ElementTree as ET
    for i in range(0, len(taxids), 150):
        x = ET.fromstring(get(EUTILS + "efetch.fcgi?db=taxonomy&id=" + ",".join(taxids[i:i + 150])))
        for t in x.findall("Taxon"):
            lin = {"species": t.findtext("ScientificName")}
            names = [t.findtext("ScientificName")]
            for a in t.findall("LineageEx/Taxon"):
                rk = a.findtext("Rank"); nm = a.findtext("ScientificName"); names.append(nm)
                if rk in RANKS:
                    lin[rk] = nm
            lin["_all"] = set(names) | {a.findtext("ScientificName") for a in t.findall("LineageEx/Taxon")}
            tax[t.findtext("TaxId")] = lin
        time.sleep(0.4)
    for acc, (tid, title) in acc2tax.items():
        out[acc] = dict(tax.get(tid, {"_all": set()}), title=title)
    return out

def verdict(lin):
    a = lin.get("_all", set())
    if "Onychophora" in a: return "host (velvet worm)"
    if a & CONTAM: return "lab contamination"
    if a & MARINE: return "implausible (marine)"
    if "Metazoa" not in a: return "non-animal"
    return "PREY"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tables", nargs="+"); ap.add_argument("-o", "--out", default="trial5_results")
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    queries, best, allhits = read_tables(a.tables)
    print(f"[trial5] {len(queries)} sequences in BLAST table(s), {len(best)} with hits. Looking up names at NCBI...")
    accs = {h[1] for hs in allhits.values() for h in sorted(hs, reverse=True)[:10]}
    lin = lineages(accs)
    rows = []
    for q in queries:
        sample = q.split("|")[0]
        if q not in best:
            rows.append(dict(sample=sample, seq=q, verdict="no BLAST hit")); continue
        top = sorted(allhits[q], reverse=True)
        b = top[0]; L = lin.get(b[1], {"_all": set()}); v = verdict(L)
        # host wins any tie with a velvet worm hit of equal bit score
        if v == "PREY" and any(h[0] >= b[0] and "Onychophora" in lin.get(h[1], {}).get("_all", set()) for h in top):
            v = "host (velvet worm)"
        level = "family" if (b[3] >= 150 and b[2] >= 95) else "order"
        rows.append(dict(sample=sample, seq=q, verdict=v, top_hit=b[1], top_hit_name=L.get("species", ""),
                         pident=b[2], aln_len=b[3], bitscore=b[0], phylum=L.get("phylum", ""),
                         **{"class": L.get("class", "")}, order=L.get("order", ""), family=L.get("family", ""),
                         report_level=level if v == "PREY" else ""))
    cols = ["sample", "seq", "verdict", "top_hit", "top_hit_name", "pident", "aln_len", "bitscore",
            "phylum", "class", "order", "family", "report_level"]
    with open(f"{a.out}/trial5_all_sequences.tsv", "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in rows: o.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    prey = [r for r in rows if r["verdict"] == "PREY"]
    with open(f"{a.out}/trial5_verified_prey.tsv", "w") as o:
        o.write("\t".join(cols) + "\n")
        for r in prey: o.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    with open(f"{a.out}/trial5_filtering_table.tsv", "w") as o:
        o.write("sample\tverdict\tn_sequences\n")
        for (s, v), n in sorted(Counter((r["sample"], r["verdict"]) for r in rows).items()):
            o.write(f"{s}\t{v}\t{n}\n")
    print("\nTrial 5 summary")
    for s in sorted({r["sample"] for r in rows}):
        c = Counter(r["verdict"] for r in rows if r["sample"] == s)
        print(f"  {s}: " + ", ".join(f"{k} {v}" for k, v in c.most_common()))
    print("\nVerified prey:")
    for r in prey:
        print(f"  {r['sample']}: {r['class']} > {r['order']} > {r['family']} | best hit {r['top_hit_name']} "
              f"({r['pident']:.1f}% over {r['aln_len']} bp) -> report at {r['report_level']} level")
    if not prey: print("  none")
    print(f"\nFiles written to {a.out}/")

if __name__ == "__main__":
    main()
