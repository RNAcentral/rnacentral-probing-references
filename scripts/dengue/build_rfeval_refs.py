#!/usr/bin/env python3
"""Build rf-eval reference structures for the release27 dengue serotypes from Rfam.

Per serotype: cmsearch the flavivirus CMs against the EDEN genome, keep one
GA-included hit per element, cut each hit span out, cmalign it back to its model
and project the family SS_cons onto it (project_ss_cons.py). Every element is
windowed on the genome record (_ncbi); 3' UTR elements additionally on the
RNAcentral sfRNA transcript (_rnac).
The xrRNA families (RF01415/RF03547) never reach GA on dengue and are left out.
"""
import os, subprocess, sys
from project_ss_cons import pair_columns, project, read_pfam_sto

HERE = os.path.dirname(os.path.abspath(__file__))
R27 = os.path.dirname(HERE)
IMG = "rnacentral/r2dt:latest"

# family -> element labels in genome order
ELEMENTS = {
    "RF02340": ["SLA"], "RF03546": ["5UTR"], "RF00617": ["cHP"],
    "RF00525": ["DB1", "DB2"], "RF00185": ["3SL"],
}


def docker(*cmd):
    subprocess.run(["docker", "run", "--rm", "-u", f"{os.getuid()}:{os.getgid()}",
                    "-v", f"{R27}:/w", "-w", "/w/dengue", IMG, *cmd], check=True)


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
        yield dict(fam=f[3], start=int(f[7]), end=int(f[8]), score=float(f[14]),
                   evalue=float(f[15]), inc=f[16] == "!")


def min_loop(dot, n=3):
    """Unpair the innermost pairs of any hairpin whose loop is shorter than n nt."""
    out, stack = list(dot), []
    for i, c in enumerate(dot):
        if c == "(":
            stack.append(i)
        elif c == ")":
            j = stack.pop()
            if all(x == "." for x in out[j + 1:i]) and i - j - 1 < n:
                out[j] = out[i] = "."
    return "".join(out)


def pick_hits(hits):
    """One GA-included hit per element; families with two copies take the two best by E."""
    chosen = {}
    for fam, labels in ELEMENTS.items():
        fam_hits = sorted((h for h in hits if h["fam"] == fam and h["inc"]), key=lambda h: h["evalue"])
        for label, h in zip(labels, sorted(fam_hits[:len(labels)], key=lambda h: h["start"])):
            chosen[label] = h
    return chosen


def main():
    for d in ("search", "in", "sto"):
        os.makedirs(os.path.join(HERE, d), exist_ok=True)
    cms = sorted(f"{fam}.cm" for fam in ELEMENTS)
    with open(os.path.join(HERE, "flavi.cm"), "w") as out:
        for f in cms:
            out.write(open(os.path.join(HERE, "cm", f)).read())
    docker("cmpress", "-F", "flavi.cm")

    db, windows, manifest = [], [], []
    for row in list(open(os.path.join(HERE, "MANIFEST-genome.tsv")))[1:]:
        dataset, sero, strain, acc, glen, urs, ulen, ident, coords = row.split()
        utr_start = int(coords.split("-")[0])
        fasta = read_fasta(os.path.join(HERE, f"{sero}_{strain}.fa"))
        genome = fasta[acc]
        with open(os.path.join(HERE, "in", f"{sero}_genome.fa"), "w") as fh:
            fh.write(f">{acc}\n{genome}\n")
        docker("cmsearch", "--cpu", "4", "--max", "--toponly", "--cut_ga", "--tblout", f"search/{sero}.tbl",
               "-o", f"search/{sero}.log", "flavi.cm", f"in/{sero}_genome.fa")
        chosen = pick_hits(list(read_tbl(os.path.join(HERE, f"search/{sero}.tbl"))))

        for label, h in sorted(chosen.items(), key=lambda kv: kv[1]["start"]):
            sid = f"{sero}_{label}"
            with open(os.path.join(HERE, "in", f"{sid}.fa"), "w") as fh:
                fh.write(f">{sid}\n{genome[h['start'] - 1:h['end']]}\n")
            docker("cmalign", "--outformat", "Pfam", "-o", f"sto/{sid}.sto",
                   f"cm/{h['fam']}.cm", f"in/{sid}.fa")
            seqs, ss_cons = read_pfam_sto(os.path.join(HERE, "sto", f"{sid}.sto"))
            seq, dot, lo, hi = project(seqs[sid], pair_columns(ss_cons))
            seq, dot = seq[lo:hi + 1], min_loop(dot[lo:hi + 1])
            gstart, gend = h["start"] + lo, h["start"] + hi
            # Every element is scored on the genome record; 3' UTR elements also on the
            # sfRNA transcript, so the two records can be compared on the same structure
            targets = [("ncbi", acc, gstart, gend)]
            if gstart >= utr_start:
                targets.append(("rnac", urs, gstart - utr_start + 1, gend - utr_start + 1))
            for source, ref, start, end in targets:
                db.append(f">{sid}_{source}\n{seq}\n{dot}\n")
                windows.append((ref, start, end, f"{sid}_{source}"))
                manifest.append((f"{sid}_{source}", sero, h["fam"], acc, gstart, gend, ref, start, end,
                                 len(seq), dot.count("("), f"{h['evalue']:.2g}"))
            print(f"[{sero}] {label:7s} {h['fam']} {acc}:{gstart}-{gend} "
                  f"{len(seq)} nt {dot.count('(')} bp E={h['evalue']:.2g}")

    with open(os.path.join(HERE, "dengue_elements.db"), "w") as fh:
        fh.writelines(db)
    with open(os.path.join(HERE, "dengue_elements.windows.tsv"), "w") as fh:
        fh.write("# ref_seq_id\tstart\tend\tstructure_id\n")
        for w in windows:
            fh.write("\t".join(map(str, w)) + "\n")
    with open(os.path.join(HERE, "MANIFEST.tsv"), "w") as fh:
        fh.write("structure\tserotype\trfam\tgenome\tgenome_start\tgenome_end\tref_seq_id\t"
                 "start\tend\tlength\tbase_pairs\tEvalue\n")
        for m in manifest:
            fh.write("\t".join(map(str, m)) + "\n")


if __name__ == "__main__":
    main()
