#!/usr/bin/env python3
# rf-eval set for genome-route organisms: the rrna_common structures plus every GA hit of each small-RNA family; from the repo root: rfeval_genome_families.py <organism> <release>
import gzip, os, re, subprocess, sys, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from project_ss_cons import pair_columns, project, read_pfam_sto
from rfeval_rrna_common import CMS, IMG, MIN_BP_PER_NT, unpair_orphans

ORGANISMS = {
    "saccharomyces_cerevisiae": {
        "genome": "https://ftp.ensembl.org/pub/current_fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz",
        "gtf": "saccharomyces_cerevisiae.R64-1-1.filtered.gtf.gz",
        "families": {"RF00005": "tRNA", "RF00003": "U1_snRNA", "RF00004": "U2_snRNA", "RF00015": "U4_snRNA", "RF00020": "U5_snRNA",
                     "RF00026": "U6_snRNA", "RF00009": "RNase_P_RNA", "RF00030": "RNase_MRP_RNA"},
    },
}


def read_genome(path):
    seqs, name, buf = {}, None, []
    for line in gzip.open(path, "rt"):
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip().upper())
    seqs[name] = "".join(buf)
    return seqs


def transcripts(gtf, genome):
    exons, strand = {}, {}
    for line in gzip.open(gtf, "rt"):
        f = line.rstrip("\n").split("\t")
        if len(f) > 8 and f[2] == "exon":
            tid = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
            exons.setdefault(tid, []).append((f[0], int(f[3]), int(f[4])))
            strand[tid] = f[6]
    comp = str.maketrans("ACGTN", "TGCAN")
    out = {}
    for tid, ex in exons.items():
        seq = "".join(genome[c][s - 1:e] for c, s, e in sorted(ex, key=lambda x: x[1]))
        out[tid] = seq.translate(comp)[::-1] if strand[tid] == "-" else seq
    return out


def main():
    name, release = sys.argv[1], sys.argv[2]
    org, orgdir = ORGANISMS[name], os.path.join(f"release{release}", name)
    work, out = os.path.abspath(os.path.join(orgdir, "work", "rfeval_families")), os.path.join(orgdir, "rfeval")
    os.makedirs(work, exist_ok=True)
    gpath = os.path.join(work, "genome.fa.gz")
    if not os.path.exists(gpath):
        urllib.request.urlretrieve(org["genome"], gpath)
    seqs = transcripts(os.path.join(orgdir, "reference", org["gtf"]), read_genome(gpath))
    # every copy of a URS gets the same multimapped reads on the genome route, so score its lowest-numbered copy only
    reps = {}
    for tid in sorted(seqs, key=lambda t: (t.rsplit(".", 1)[0], int(t.rsplit(".", 1)[1]))):
        reps.setdefault(tid.rsplit(".", 1)[0], tid)
    with open(os.path.join(work, "records.fa"), "w") as fh:
        fh.writelines(f">{t}\n{seqs[t]}\n" for t in reps.values())
    docker = lambda *cmd: subprocess.run(["docker", "run", "--rm", "-u", f"{os.getuid()}:{os.getgid()}", "-v", f"{work}:/w",
                                          "-v", f"{CMS}:/cm:ro", "-w", "/w", IMG, *cmd], check=True, capture_output=True)

    db = [f">{b}" for b in open(os.path.join(out, "rrna_common.db")).read().split(">")[1:]]
    windows = [l for l in open(os.path.join(out, "rrna_common.windows.tsv")) if not l.startswith("#")]
    manifest = []
    for acc, molecule in org["families"].items():
        tbl = os.path.join(work, f"{acc}.tbl")
        if not os.path.exists(tbl):
            docker("cmsearch", "--cpu", "8", "--cut_ga", "--nohmmonly", "--toponly", "--tblout", f"{acc}.tbl", f"/cm/{acc}.cm", "records.fa")
        hits = [l.split() for l in open(tbl) if not l.startswith("#")]
        best = {}
        for h in sorted(((float(h[14]), h[0], int(h[7]), int(h[8])) for h in hits), key=lambda h: (-h[0], h[1])):
            best.setdefault(h[1], h)
        if not best:
            print(f"  {molecule:14s} no GA hit"); continue
        with open(os.path.join(work, f"{acc}.fa"), "w") as fh:
            fh.writelines(f">{rec}\n{seqs[rec][start - 1:end]}\n" for _, rec, start, end in best.values())
        docker("cmalign", "--cpu", "8", "--outformat", "Pfam", "-o", f"{acc}.sto", f"/cm/{acc}.cm", f"{acc}.fa")
        aligned, ss_cons = read_pfam_sto(os.path.join(work, f"{acc}.sto"))
        kept = 0
        for score, rec, start, end in best.values():
            seq, dot, lo, hi = project(aligned[rec], pair_columns(ss_cons))
            if lo is None:
                continue
            seq, dot = seq[lo:hi + 1], unpair_orphans(dot[lo:hi + 1])
            sid, bp = f"{rec}__{molecule}", dot.count("(")
            if bp / len(seq) < MIN_BP_PER_NT:
                continue
            kept += 1
            db.append(f">{sid}\n{seq}\n{dot}\n")
            windows.append(f"{rec}\t{start + lo}\t{start + hi}\t{sid}\n")
            manifest.append(f"{sid}\t{acc}\t{rec}\t{start + lo}\t{start + hi}\t{len(seq)}\t{bp}\t{score:.1f}\n")
        print(f"  {molecule:14s} {len(best)} records hit, {kept} structures kept")

    open(os.path.join(out, "ncrna_families.db"), "w").writelines(db)
    open(os.path.join(out, "ncrna_families.windows.tsv"), "w").writelines(["# ref_seq_id\tstart\tend\tstructure_id\n"] + windows)
    open(os.path.join(orgdir, "MANIFEST-rfeval-ncrna_families.tsv"), "w").writelines(
        ["structure\trfam\tref_seq_id\tstart\tend\tlength\tbase_pairs\tbits\n"] + manifest)
    print(f"{len(db)} structures in all ({len(db) - len(manifest)} from rrna_common)")


if __name__ == "__main__":
    main()
