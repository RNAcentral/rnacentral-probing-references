#!/usr/bin/env python3
# rf-eval set for the small-RNA transcriptome, one representative structure per family; from the repo root: rfeval_small_rna.py <organism> <release>
import os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from project_ss_cons import pair_columns, project, read_pfam_sto
from rfeval_rrna_common import CMS, IMG, MIN_BP_PER_NT, read_fasta, unpair_orphans

ORGANISMS = {
    "homo_sapiens": {"RF00005": "tRNA", "RF00001": "5S_rRNA", "RF01960": "18S_rRNA",
                     "RF00003": "U1_snRNA", "RF00017": "7SL_RNA"},
}


def main():
    name, release = sys.argv[1], sys.argv[2]
    orgdir = os.path.join(f"release{release}", name)
    work, out = os.path.abspath(os.path.join(orgdir, "work", "rfeval_small_rna")), os.path.join(orgdir, "rfeval")
    os.makedirs(work, exist_ok=True)
    fa = os.path.join(orgdir, "reference", f"{name}.small_rna.transcriptome.fa")
    seqs = read_fasta(fa)
    with open(os.path.join(work, "records.fa"), "w") as fh:
        fh.writelines(f">{k}\n{v}\n" for k, v in seqs.items())
    docker = lambda *cmd: subprocess.run(["docker", "run", "--rm", "-u", f"{os.getuid()}:{os.getgid()}", "-v", f"{work}:/w",
                                          "-v", f"{CMS}:/cm:ro", "-w", "/w", IMG, *cmd], check=True, capture_output=True)

    db, windows, manifest = [], [], []
    for acc, molecule in ORGANISMS[name].items():
        tbl = os.path.join(work, f"{acc}.tbl")
        if not os.path.exists(tbl):
            # --toponly: only sense-strand hits are measured on a transcript
            docker("cmsearch", "--cpu", "8", "--cut_ga", "--nohmmonly", "--toponly", "--tblout", f"{acc}.tbl", f"/cm/{acc}.cm", "records.fa")
        hits = [l.split() for l in open(tbl) if not l.startswith("#")]
        if not hits:
            print(f"  {molecule}: no GA hit"); continue
        # the representative is the best-scoring hit; ties go to the shorter record, then the ID
        score, _, rec, start, end = min(((float(h[14]), len(seqs[h[0]]), h[0], int(h[7]), int(h[8])) for h in hits),
                                        key=lambda h: (-h[0], h[1], h[2]))
        print(f"  {molecule:9s} {len(hits)} GA hits, best {rec} {start}-{end} {score:.1f} bits")
        with open(os.path.join(work, f"{acc}.fa"), "w") as fh:
            fh.write(f">{rec}\n{seqs[rec][start - 1:end]}\n")
        docker("cmalign", "--cpu", "8", "--outformat", "Pfam", "-o", f"{acc}.sto", f"/cm/{acc}.cm", f"{acc}.fa")
        aligned, ss_cons = read_pfam_sto(os.path.join(work, f"{acc}.sto"))
        seq, dot, lo, hi = project(aligned[rec], pair_columns(ss_cons))
        seq, dot = seq[lo:hi + 1], unpair_orphans(dot[lo:hi + 1])
        sid, bp = f"{rec}__{molecule}", dot.count("(")
        print(f"  -> {sid} {start + lo}-{start + hi}: {bp} bp over {len(seq)} nt ({bp / len(seq):.2f} bp/nt)")
        if bp / len(seq) < MIN_BP_PER_NT:
            print(f"  {molecule}: below {MIN_BP_PER_NT} bp/nt, dropped"); continue
        db.append(f">{sid}\n{seq}\n{dot}\n")
        windows.append(f"{rec}\t{start + lo}\t{start + hi}\t{sid}\n")
        manifest.append(f"{sid}\t{acc}\t{rec}\t{start + lo}\t{start + hi}\t{len(seq)}\t{bp}\t{score:.1f}\t{len(hits)}\n")

    open(os.path.join(out, "small_rna.db"), "w").writelines(db)
    open(os.path.join(out, "small_rna.windows.tsv"), "w").writelines(["# ref_seq_id\tstart\tend\tstructure_id\n"] + windows)
    open(os.path.join(orgdir, "MANIFEST-rfeval-small_rna.tsv"), "w").writelines(
        ["structure\trfam\tref_seq_id\tstart\tend\tlength\tbase_pairs\tbits\tga_hits\n"] + manifest)


if __name__ == "__main__":
    main()
