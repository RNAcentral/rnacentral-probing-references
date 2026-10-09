#!/usr/bin/env python3
"""Write reference/<stem>.transcriptome.fa from each work/<stem>.gtf + .fa under the given organism dirs.

Each GTF transcript is cut from the genome FASTA. Contigs not spanned whole by a
transcript are dropped (bacterial chromosome); contigs a transcript spans whole are
kept once under the contig ID (viral genomes, so rf-eval windows still match).
"""
import re
import sys
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq


def transcripts(gtf):
    exons = defaultdict(list)
    for line in open(gtf):
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        tid = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
        exons[tid].append((f[0], int(f[3]), int(f[4]), f[6]))
    return exons


def build(gtf):
    stem = gtf.with_suffix("")
    fa = next((p for p in (stem.with_suffix(".fa"), stem.with_suffix(".fasta")) if p.exists()), None)
    if fa is None:
        print(f"SKIP {gtf}: no matching FASTA")
        return
    genome = {r.id: str(r.seq) for r in SeqIO.parse(fa, "fasta")}
    out, whole = [], set()
    for tid, ex in transcripts(gtf).items():
        contig, strand = ex[0][0], ex[0][3]
        ex.sort(key=lambda e: e[1])
        seq = "".join(genome[contig][s - 1:e] for _, s, e, _ in ex)
        if strand == "-":
            seq = str(Seq(seq).reverse_complement())
        if seq == genome[contig]:
            whole.add(contig)
        else:
            out.append((tid, seq))
    records = [(c, genome[c]) for c in genome if c in whole] + out
    dest = gtf.parent.parent / "reference" / f"{stem.name}.transcriptome.fa"
    dest.parent.mkdir(exist_ok=True)
    with open(dest, "w") as fh:
        for rid, seq in records:
            fh.write(f">{rid}\n" + "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\n")
    dropped = [c for c in genome if c not in whole]
    print(f"{dest}: {len(whole)} whole-contig + {len(out)} records; dropped contigs: {dropped}")


for d in sys.argv[1:]:
    for gtf in sorted(Path(d, "work").glob("*.gtf")):
        build(gtf)
