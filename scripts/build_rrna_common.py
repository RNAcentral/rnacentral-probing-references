#!/usr/bin/env python3
"""Build a common rRNA rf-eval reference, same shape as human/rrna_common.

The cytoplasmic rRNAs (18S/5.8S/25S-28S/5S, or the E. coli rrn operon's
16S/23S/5S), found by cmsearch --cut_ga against the
transcripts a run is keyed on, one locus per molecule, SS_cons projected.

    python3 build_rrna_common.py yeast|mouse|ecoli
"""
import json, os, subprocess, sys, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
from project_ss_cons import pair_columns, project, read_pfam_sto

IMG = "rnacentral/r2dt:latest"
CMS = os.path.join(HERE, "rrna_cm")
GTFS = os.path.join(os.path.dirname(HERE), "rnac-release27")
EUK = {"RF01960": "18S_rRNA", "RF00002": "5.8S_rRNA", "RF02543": "28S_rRNA", "RF00001": "5S_rRNA"}
ORGANISMS = {
    "yeast": {"species": "saccharomyces_cerevisiae", "gtf": "saccharomyces_cerevisiae.R64-1-1.release27.filtered_v5.gtf",
              "families": {**EUK, "RF02543": "25S_rRNA"}},
    # mouse 5S genes are a separate cluster, not in the rDNA unit
    "mouse": {"species": "mus_musculus", "gtf": "mus_musculus.GRCm39.release27.filtered_v5.gtf", "families": EUK,
              "unlinked": {"RF00001"}},
    # windows are moved onto the RNAcentral records the run's GTF places on the genome
    "ecoli": {"ncbi": "NC_000913.3", "remap": "ECOLI_MG1655.gtf", "families": {"RF00177": "16S_rRNA", "RF02541": "23S_rRNA", "RF00001": "5S_rRNA"}},
}


def post_json(url, body):
    req = urllib.request.Request(url, json.dumps(body).encode(), {"Content-Type": "application/json", "Accept": "application/json"})
    return json.load(urllib.request.urlopen(req))


def gtf_rrna_exons(path):
    """rRNA-biotype transcripts >= 100 nt: {transcript_id: [(chrom, start, end, strand), ...]} in transcript order."""
    exons = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon" or 'transcript_biotype "rRNA"' not in f[8]:
            continue
        tx = f[8].split('transcript_id "', 1)[1].split('"', 1)[0]
        exons.setdefault(tx, []).append((f[0], int(f[3]), int(f[4]), f[6]))
    for tx, ex in exons.items():
        ex.sort(key=lambda e: e[1], reverse=ex[0][3] == "-")
    return {tx: ex for tx, ex in exons.items() if sum(e[2] - e[1] + 1 for e in ex) >= 100}


def fetch_candidates(org, fa):
    """Write the candidate sequences; Ensembl REST spells the same assembly the pipeline downloads."""
    if "ncbi" in org:
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={org['ncbi']}&rettype=fasta&retmode=text"
        seq = "".join(urllib.request.urlopen(url).read().decode().splitlines()[1:])
        open(fa, "w").write(f">{org['ncbi']}\n{seq}\n")
        return
    exons = gtf_rrna_exons(os.path.join(GTFS, org["gtf"]))
    regions = sorted({f"{c}:{s}..{e}:{1 if st == '+' else -1}" for ex in exons.values() for c, s, e, st in ex})
    seqs = {}
    for i in range(0, len(regions), 50):
        for r in post_json(f"https://rest.ensembl.org/sequence/region/{org['species']}", {"regions": regions[i:i + 50]}):
            seqs[r["query"]] = r["seq"].upper()
    with open(fa, "w") as fh:
        for tx, ex in sorted(exons.items()):
            fh.write(f">{tx}\n" + "".join(seqs[f"{c}:{s}..{e}:{1 if st == '+' else -1}"] for c, s, e, st in ex) + "\n")


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]; seqs[name] = []
        elif name:
            seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def main():
    name = sys.argv[1]
    org, out = ORGANISMS[name], os.path.join(HERE, name, "rrna_common")
    work = os.path.join(out, "work")
    os.makedirs(work, exist_ok=True)
    docker = lambda *cmd: subprocess.run(["docker", "run", "--rm", "-u", f"{os.getuid()}:{os.getgid()}", "-v", f"{work}:/w",
                                          "-v", f"{CMS}:/cm:ro", "-w", "/w", IMG, *cmd], check=True, capture_output=True)
    fa = os.path.join(work, "candidates.fa")
    if not os.path.exists(fa):
        fetch_candidates(org, fa)
    seqs = read_fasta(fa)
    print(f"{name}: {len(seqs)} candidate sequences")
    exons = gtf_rrna_exons(os.path.join(GTFS, org["gtf"])) if "gtf" in org else {}
    locus = lambda tx, s, e: (exons[tx][0][0], min(x[1] for x in exons[tx]), max(x[2] for x in exons[tx])) if exons else (tx, s, e)

    remap = {}
    if "remap" in org:
        ids = {l.split('transcript_id "')[1].split('"')[0] for l in open(os.path.join(HERE, name, org["remap"]))}
        remap = {k: v for k, v in read_fasta(os.path.join(HERE, name, "records.fa")).items() if k in ids}
    db, windows, anchor = [], [], None
    for acc, molecule in org["families"].items():
        # --toponly: only sense-strand hits are measured on a transcript
        tbl = os.path.join(work, f"{acc}.tbl")
        if not os.path.exists(tbl):
            docker("cmsearch", "--cpu", "8", "--cut_ga", "--nohmmonly", "--toponly", "--tblout", f"{acc}.tbl", f"/cm/{acc}.cm", "candidates.fa")
        hits = [l.split() for l in open(tbl) if not l.startswith("#")]
        hits = [(float(h[14]), len(seqs[h[0]]), h[0], int(h[7]), int(h[8])) for h in hits]
        for h in sorted(hits, key=lambda h: (-h[0], h[1], h[2])):
            print(f"  {molecule:10s} {h[2]:32s} {h[3]}-{h[4]} {h[0]:.1f} bits (transcript {h[1]} nt)")
        if anchor and acc not in org.get("unlinked", ()):
            # the rest must come from the SSU's own rDNA unit/operon, not a pseudogene or the mitochondrion
            hits = [h for h in hits if (l := locus(*h[2:])) and l[0] == anchor[0] and max(l[1], anchor[1]) - min(l[2], anchor[2]) < 10000]
        if not hits:
            print(f"  {molecule}: no GA hit in the rDNA unit"); continue
        # identical copies tie on score; the shortest transcript is the mature molecule itself
        score, _, tx, start, end = min(hits, key=lambda h: (-h[0], h[1], h[2]))
        anchor = anchor or locus(tx, start, end)
        open(os.path.join(work, f"{acc}.fa"), "w").write(f">{tx}\n{seqs[tx][start - 1:end]}\n")
        docker("cmalign", "--cpu", "8", "--outformat", "Pfam", "-o", f"{acc}.sto", f"/cm/{acc}.cm", f"{acc}.fa")
        aligned, ss_cons = read_pfam_sto(os.path.join(work, f"{acc}.sto"))
        seq, dot, lo, hi = project(aligned[tx], pair_columns(ss_cons))
        seq, dot = seq[lo:hi + 1], dot[lo:hi + 1]
        # a pair whose partner fell outside the trimmed span is no longer a pair
        stack, fixed = [], list(dot)
        for i, ch in enumerate(dot):
            if ch == "(":
                stack.append(i)
            elif ch == ")":
                if stack: stack.pop()
                else: fixed[i] = "."
        for i in stack:
            fixed[i] = "."
        dot, sid = "".join(fixed), f"{tx}__{molecule}"
        print(f"  -> {sid} {start + lo}-{start + hi}: {dot.count('(')} bp over {len(seq)} nt ({dot.count('(') / len(seq):.2f} bp/nt)")
        db.append(f">{sid}\n{seq}\n{dot}\n")
        ref, ref_start = tx, start + lo
        if remap:
            dna = seq.upper().replace("U", "T")
            (ref, ref_start), = [(k, v.upper().find(dna) + 1) for k, v in remap.items() if dna in v.upper()]
        windows.append(f"{ref}\t{ref_start}\t{ref_start + len(seq) - 1}\t{sid}\n")

    open(os.path.join(out, "rrna_common.db"), "w").writelines(db)
    open(os.path.join(out, "rrna_common.windows.tsv"), "w").writelines(["# ref_seq_id\tstart\tend\tstructure_id\n"] + windows)


if __name__ == "__main__":
    main()
