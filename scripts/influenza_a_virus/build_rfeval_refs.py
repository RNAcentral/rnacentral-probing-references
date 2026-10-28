#!/usr/bin/env python3
"""Build the rf-eval reference for A/WSN/1933 from Rfam.

Rfam has a single influenza A family, RF01099 PK-IAV (a 48-nt pseudoknot). It is
searched against the eight segments, GA-included hits are cut out, cmaligned back
to the model with the seed structure mapped in (--mapali/--mapstr, so the pseudoknot
survives) and the SS_cons projected onto the hit: nested pairs as (), pseudoknot
pairs as []. Unlike the dengue build the pseudoknot is kept, because for this
family it is the structure, and rf-eval runs with -kp.
"""
import os, subprocess
from project_ss_cons import pair_columns, project, read_pfam_sto

HERE = os.path.dirname(os.path.abspath(__file__))
R27 = os.path.dirname(HERE)
IMG = "rnacentral/r2dt:latest"
FAM = "RF01099"
STRAIN = "IAV_WSN1933"


def docker(*cmd):
    subprocess.run(["docker", "run", "--rm", "-u", f"{os.getuid()}:{os.getgid()}",
                    "-v", f"{R27}:/w", "-w", "/w/influenza", IMG, *cmd], check=True)


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]; seqs[name] = []
        else:
            seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def read_tbl(path):
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.split()
        yield dict(ref=f[0], start=int(f[7]), end=int(f[8]), score=float(f[14]),
                   evalue=float(f[15]), inc=f[16] == "!")


def pk_columns(ss_cons):
    """Pseudoknot pairs from WUSS letters: an uppercase run opens, its lowercase run closes."""
    opens, pairs = {}, {}
    for col, ch in enumerate(ss_cons):
        if ch.isupper():
            opens.setdefault(ch, []).append(col)
        elif ch.islower():
            i = opens[ch.upper()].pop()
            pairs[i], pairs[col] = col, i
    return pairs


def main():
    for d in ("search", "in", "sto"):
        os.makedirs(os.path.join(HERE, d), exist_ok=True)
    fasta = read_fasta(os.path.join(HERE, f"{STRAIN}.fa"))
    segment = {row.split()[3]: row.split()[2] for row in list(open(os.path.join(HERE, "MANIFEST.tsv")))[1:]}
    docker("cmsearch", "--cpu", "4", "--max", "--toponly", "--cut_ga", "--tblout", f"search/{FAM}.tbl",
           "-o", f"search/{FAM}.log", f"cm/{FAM}.cm", f"{STRAIN}.fa")

    db, windows, manifest = [], [], []
    for h in sorted((h for h in read_tbl(os.path.join(HERE, "search", f"{FAM}.tbl")) if h["inc"]),
                    key=lambda h: (h["ref"], h["start"])):
        sid = f"{STRAIN}_{segment[h['ref']]}_PK"
        with open(os.path.join(HERE, "in", f"{sid}.fa"), "w") as fh:
            fh.write(f">{sid}\n{fasta[h['ref']][h['start'] - 1:h['end']]}\n")
        docker("cmalign", "--outformat", "Pfam", "--mapali", f"cm/{FAM}.seed.sto", "--mapstr",
               "-o", f"sto/{sid}.sto", f"cm/{FAM}.cm", f"in/{sid}.fa")
        seqs, ss_cons = read_pfam_sto(os.path.join(HERE, "sto", f"{sid}.sto"))
        nested = "".join("." if c.isalpha() else c for c in ss_cons)
        seq, dot, lo, hi = project(seqs[sid], pair_columns(nested))
        _, pk, _, _ = project(seqs[sid], pk_columns(ss_cons))
        dot = "".join(p.translate(str.maketrans("()", "[]")) if p != "." else d for d, p in zip(dot, pk))
        seq, dot = seq[lo:hi + 1], dot[lo:hi + 1]
        start, end = h["start"] + lo, h["start"] + hi
        db.append(f">{sid}\n{seq}\n{dot}\n")
        windows.append((h["ref"], start, end, sid))
        manifest.append((sid, FAM, h["ref"], start, end, len(seq),
                         dot.count("("), dot.count("["), f"{h['evalue']:.2g}"))
        print(f"{sid} {h['ref']}:{start}-{end} {len(seq)} nt {dot.count('(')} bp + {dot.count('[')} pk bp E={h['evalue']:.2g}")
        print(f"  {seq}\n  {dot}")

    with open(os.path.join(HERE, "iav_elements.db"), "w") as fh:
        fh.writelines(db)
    with open(os.path.join(HERE, "iav_elements.windows.tsv"), "w") as fh:
        fh.write("# ref_seq_id\tstart\tend\tstructure_id\n")
        for w in windows:
            fh.write("\t".join(map(str, w)) + "\n")
    with open(os.path.join(HERE, "MANIFEST-rfeval.tsv"), "w") as fh:
        fh.write("structure\trfam\tref_seq_id\tstart\tend\tlength\tbase_pairs\tpk_pairs\tEvalue\n")
        for m in manifest:
            fh.write("\t".join(map(str, m)) + "\n")


if __name__ == "__main__":
    main()
