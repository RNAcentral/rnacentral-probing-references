#!/usr/bin/env python3
"""Build rf-eval reference structures for the three Zika references from Rfam.

Per strain: cmsearch the flavivirus CMs against the genome at GA (--cut_ga), cut every
included hit out, cmalign it back to its model with the seed structure mapped in
(--mapali/--mapstr) and project the SS_cons: nested pairs as (), pseudoknot pairs as [].
Every element is windowed on the genome (_ncbi); an element inside one of the RNAcentral
records appended to the FASTA (none qualified for Zika) would be windowed on it too.
"""
import os, subprocess
from project_ss_cons import pair_columns, project, read_pfam_sto

HERE = os.path.dirname(os.path.abspath(__file__))
R27 = os.path.dirname(HERE)
IMG = "rnacentral/r2dt:latest"
FAMILIES = {"RF02340": "SLA", "RF03546": "5UTR", "RF00617": "cHP", "RF03547": "xrRNA",
            "RF01415": "SLIV", "RF00525": "DB", "RF00185": "3SL"}
RECORDS = {}


def docker(*cmd):
    subprocess.run(["docker", "run", "--rm", "-u", f"{os.getuid()}:{os.getgid()}",
                    "-v", f"{R27}:/w", "-w", "/w/zika", IMG, *cmd], check=True)


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
        yield dict(fam=f[3], start=int(f[7]), end=int(f[8]), evalue=float(f[15]), inc=f[16] == "!")


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


def clip(dot):
    """Unpair any bracket whose partner was cut off with the flank."""
    out, stacks = list(dot), {"(": [], "[": []}
    close = {")": "(", "]": "["}
    for i, c in enumerate(dot):
        if c in stacks:
            stacks[c].append(i)
        elif c in close:
            if stacks[close[c]]:
                stacks[close[c]].pop()
            else:
                out[i] = "."
    for st in stacks.values():
        for i in st:
            out[i] = "."
    return "".join(out)


def main():
    for d in ("search", "in", "sto"):
        os.makedirs(os.path.join(HERE, d), exist_ok=True)
    with open(os.path.join(HERE, "flavi.cm"), "w") as out:
        for fam in FAMILIES:
            out.write(open(os.path.join(HERE, "cm", f"{fam}.cm")).read())
    docker("cmpress", "-F", "flavi.cm")

    # strain -> (accession, [(record, start, end)]) from the reference build's manifest
    strains = {}
    for row in list(open(os.path.join(HERE, "MANIFEST-genome.tsv")))[1:]:
        _, label, acc, rec, _, _, start, end, _ = row.rstrip("\n").split("\t")
        strains.setdefault(label, (acc, []))
        if rec != acc:
            strains[label][1].append((rec, int(start), int(end)))

    db, windows, manifest = [], [], []
    for label, (acc, records) in strains.items():
        fasta = read_fasta(os.path.join(HERE, f"ZIKV_{label}.fa"))
        genome = fasta[acc]
        with open(os.path.join(HERE, "in", f"{label}_genome.fa"), "w") as fh:
            fh.write(f">{acc}\n{genome}\n")
        docker("cmsearch", "--cpu", "4", "--max", "--toponly", "--cut_ga", "--tblout", f"search/{label}.tbl",
               "-o", f"search/{label}.log", "flavi.cm", f"in/{label}_genome.fa")
        hits = sorted((h for h in read_tbl(os.path.join(HERE, "search", f"{label}.tbl")) if h["inc"]),
                      key=lambda h: h["start"])
        # a family hitting more than once is numbered in genome order
        counts = {}
        for h in hits:
            counts[h["fam"]] = counts.get(h["fam"], 0) + 1
        seen = {}
        for h in hits:
            fam = h["fam"]
            seen[fam] = seen.get(fam, 0) + 1
            element = FAMILIES[fam] + (str(seen[fam]) if counts[fam] > 1 else "")
            sid = f"{label}_{element}"
            with open(os.path.join(HERE, "in", f"{sid}.fa"), "w") as fh:
                fh.write(f">{sid}\n{genome[h['start'] - 1:h['end']]}\n")
            docker("cmalign", "--outformat", "Pfam", "--mapali", f"cm/{fam}.seed.sto", "--mapstr",
                   "-o", f"sto/{sid}.sto", f"cm/{fam}.cm", f"in/{sid}.fa")
            seqs, ss_cons = read_pfam_sto(os.path.join(HERE, "sto", f"{sid}.sto"))
            nested = "".join("." if c.isalpha() else c for c in ss_cons)
            seq, dot, lo, hi = project(seqs[sid], pair_columns(nested))
            _, pk, _, _ = project(seqs[sid], pk_columns(ss_cons))
            dot = "".join(p.translate(str.maketrans("()", "[]")) if p != "." else d for d, p in zip(dot, pk))
            seq, dot = seq[lo:hi + 1], dot[lo:hi + 1]
            gstart, gend = h["start"] + lo, h["start"] + hi
            # (id, seq, dot, genome span, [(source, record, start, end)])
            outputs = [(sid, seq, dot, gstart, gend, [("ncbi", acc, gstart, gend)])]
            for rec, rstart, rend in records:
                if rstart <= gstart and gend <= rend:
                    outputs[0][5].append(("rnac", rec, gstart - rstart + 1, gend - rstart + 1))
                elif gstart < rstart and rend < gend:
                    lo2, hi2 = rstart - gstart, rend - gstart
                    outputs.append((f"{sid}_{RECORDS[rec]}", seq[lo2:hi2 + 1], clip(dot[lo2:hi2 + 1]), rstart, rend,
                                    [("ncbi", acc, rstart, rend), ("rnac", rec, 1, rend - rstart + 1)]))
            for oid, oseq, odot, ostart, oend, targets in outputs:
                for source, ref, start, end in targets:
                    db.append(f">{oid}_{source}\n{oseq}\n{odot}\n")
                    windows.append((ref, start, end, f"{oid}_{source}"))
                    manifest.append((f"{oid}_{source}", label, fam, acc, ostart, oend, ref, start, end,
                                     len(oseq), odot.count("("), odot.count("["), f"{h['evalue']:.2g}"))
                print(f"[{label}] {oid[len(label) + 1:]:9s} {fam} {acc}:{ostart}-{oend} {len(oseq)} nt "
                      f"{odot.count('(')} bp + {odot.count('[')} pk E={h['evalue']:.2g} on {len(targets)} record(s)")

    with open(os.path.join(HERE, "zikv_elements.db"), "w") as fh:
        fh.writelines(db)
    with open(os.path.join(HERE, "zikv_elements.windows.tsv"), "w") as fh:
        fh.write("# ref_seq_id\tstart\tend\tstructure_id\n")
        for w in windows:
            fh.write("\t".join(map(str, w)) + "\n")
    with open(os.path.join(HERE, "MANIFEST.tsv"), "w") as fh:
        fh.write("structure\tstrain\trfam\tgenome\tgenome_start\tgenome_end\tref_seq_id\t"
                 "start\tend\tlength\tbase_pairs\tpk_pairs\tEvalue\n")
        for m in manifest:
            fh.write("\t".join(map(str, m)) + "\n")


if __name__ == "__main__":
    main()
