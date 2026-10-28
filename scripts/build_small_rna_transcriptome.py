#!/usr/bin/env python3
# Small-RNA transcriptome for size-selected libraries, from the filtered GTF; from the repo root: build_small_rna_transcriptome.py <organism> <release> [species FASTA]
import gzip, io, os, re, sys, urllib.request
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from filter_rnacentral_gff3 import _WORST, _best_rank, _length_groups, _overlap_clusters, _unique_count

# classes kept from the filtered GTF; rRNA records only when the organism lists them
BIOTYPES = {"tRNA", "snRNA", "scaRNA", "snoRNA", "pre_miRNA", "Y_RNA", "vault_RNA", "SRP_RNA",
            "RNase_P_RNA", "RNase_MRP_RNA", "telomerase_RNA"}
MAX_LEN = 500
MAX_MISMATCH = 2
FASTA_URL = "https://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/{}.0/sequences/rnacentral_species_specific_ids.fasta.gz"
ORGANISMS = {
    "homo_sapiens": {
        "gtf": "homo_sapiens.GRCh38.filtered.gtf.gz",
        # the two chr21 rDNA units rrna_common sits on, the 5S array, mt 16S and mt 12S (biotype ncRNA in the GTF)
        "rrna": ["URS0000ABD879_9606", "URS0000ABD87F_9606", "URS00000F9D45_9606", "URS000080E357_9606", "URS000044DFF6_9606"],
    },
}


def gtf_records(path, rrna):
    recs = []
    for line in gzip.open(path, "rt"):
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "transcript":
            continue
        a = dict(re.findall(r'(\S+) "([^"]*)"', f[8]))
        urs = a["transcript_id"].rsplit(".", 1)[0]
        if a.get("transcript_biotype") not in BIOTYPES and urs not in rrna:
            continue
        prov = _best_rank(a.get("source_db", ""))
        recs.append({"id": a["transcript_id"], "urs": urs, "biotype": a["transcript_biotype"], "chrom": f[0],
                     "strand": f[6], "start": int(f[3]), "end": int(f[4]), "databases": a.get("databases", ""),
                     "source_db": a.get("source_db", ""), "db_count": _unique_count(a.get("databases", "")),
                     "rank": prov if prov < _WORST else _best_rank(a.get("databases", ""))})
    return recs


def best_per_locus(recs):
    # same rule as the GTF filter's overlap thinning, applied across the kept classes
    loci = defaultdict(list)
    for r in recs:
        loci[(r["chrom"], r["strand"])].append(r)
    kept = []
    for txs in loci.values():
        for cluster in _overlap_clusters(txs):
            for group in _length_groups(cluster):
                kept.append(min(group, key=lambda x: (x["rank"], -x["db_count"], -(x["end"] - x["start"]), x["id"])))
    return kept


def read_sequences(ids, source, cache):
    seqs = {}
    if os.path.exists(cache):
        seqs = {k: v for k, v in read_fasta(cache).items() if k in ids}
    if ids - seqs.keys():
        stream = gzip.open(source, "rt") if os.path.exists(source) else io.TextIOWrapper(gzip.GzipFile(fileobj=urllib.request.urlopen(source)))
        name = None
        for line in stream:
            if line.startswith(">"):
                name = line[1:].split()[0]
                name = name if name in ids else None
                if name:
                    seqs[name] = ""
            elif name:
                seqs[name] += line.strip()
        with open(cache, "w") as fh:
            fh.writelines(f">{k}\n{v}\n" for k, v in sorted(seqs.items()))
    return seqs


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]; seqs[name] = []
        elif name:
            seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def merge_near_identical(order, seqs):
    # each record joins the first better-ranked representative within MAX_MISMATCH substitutions, else becomes one
    reps, members, index = [], {}, defaultdict(list)
    for urs in order:
        s = seqs[urs]
        cuts = [0, len(s) // 3, 2 * len(s) // 3, len(s)]
        segments = [(len(s), i, s[cuts[i]:cuts[i + 1]]) for i in range(3)]
        # with at most 2 substitutions in 3 segments, one segment matches exactly
        candidates = sorted({r for seg in segments for r in index[seg]}, key=reps.index)
        rep = next((r for r in candidates if sum(a != b for a, b in zip(s, seqs[r])) <= MAX_MISMATCH), None)
        if rep:
            members[rep].append(urs)
            continue
        reps.append(urs); members[urs] = []
        for seg in segments:
            index[seg].append(urs)
    return reps, members


def main():
    name, release = sys.argv[1], sys.argv[2]
    org, orgdir = ORGANISMS[name], os.path.join(f"release{release}", name)
    source = sys.argv[3] if len(sys.argv) > 3 else FASTA_URL.format(release)
    os.makedirs(os.path.join(orgdir, "work"), exist_ok=True)
    rrna = set(org["rrna"])

    recs = gtf_records(os.path.join(orgdir, "reference", org["gtf"]), rrna)
    print(f"{len(recs)} GTF records in the kept classes")
    kept = best_per_locus(recs)
    print(f"{len(kept)} after keeping the best record per locus")

    copies = defaultdict(list)
    for r in kept:
        copies[r["urs"]].append(r)
    best = {u: min(rs, key=lambda x: (x["rank"], -x["db_count"], x["id"])) for u, rs in copies.items()}
    for u in rrna - best.keys():
        best[u] = {"urs": u, "biotype": "rRNA", "rank": 0, "db_count": 0, "databases": "", "source_db": "", "chrom": None}
        copies[u] = []
    print(f"{len(best)} distinct sequences")

    seqs = read_sequences(set(best), source, os.path.join(orgdir, "work", "small_rna.seqs.fa"))
    missing = set(best) - seqs.keys()
    if missing:
        sys.exit(f"no sequence in the release FASTA for {len(missing)} records, e.g. {sorted(missing)[:3]}")
    seqs = {u: re.sub("[^ACGTN]", "N", s.upper().replace("U", "T")) for u, s in seqs.items()}
    too_long = {u for u in best if len(seqs[u]) > MAX_LEN and u not in rrna}
    print(f"{len(too_long)} dropped as longer than {MAX_LEN} nt")

    order = sorted((u for u in best if u not in too_long and u not in rrna), key=lambda u: (best[u]["rank"], -best[u]["db_count"], u))
    reps, members = merge_near_identical(order, seqs)
    print(f"{len(order) - len(reps)} merged into a representative within {MAX_MISMATCH} substitutions")
    records = sorted(reps + sorted(rrna))

    fa = os.path.join(orgdir, "reference", f"{name}.small_rna.transcriptome.fa")
    with open(fa, "w") as fh:
        for u in records:
            fh.write(f">{u}\n" + "\n".join(seqs[u][i:i + 60] for i in range(0, len(seqs[u]), 60)) + "\n")
    with open(os.path.join(orgdir, "MANIFEST-small_rna.tsv"), "w") as fh:
        fh.write("record\tbiotype\tlength\tsource_db\tdatabases\tgenomic_copies\tfirst_locus\tmerged\n")
        for u in records:
            b, locs = best[u], sorted(copies[u], key=lambda x: (x["chrom"], x["start"]))
            locus = f"{locs[0]['chrom']}:{locs[0]['start']}-{locs[0]['end']}:{locs[0]['strand']}" if locs else "not in GTF"
            fh.write(f"{u}\t{b['biotype']}\t{len(seqs[u])}\t{b['source_db']}\t{b['databases']}\t{len(locs)}\t{locus}\t{','.join(members.get(u, []))}\n")
    counts = Counter(best[u]["biotype"] for u in records)
    print(f"{fa}: {len(records)} records (" + ", ".join(f"{k} {v}" for k, v in sorted(counts.items())) + ")")


if __name__ == "__main__":
    main()
