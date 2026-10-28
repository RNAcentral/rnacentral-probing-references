"""Build the A/WSN/1933(H1N1) reference for rnastruct00024: the eight NCBI segment records.

Same accessions as conf/viral_genomes.config, pinned locally so the run does not depend
on an NCBI fetch. RNAcentral holds no full-length influenza A transcript (only PDB
promoter fragments), so unlike dengue there is no second record to add.
"""
import os, sys, urllib.request

STRAIN = "RVA_RF"
SEGMENTS = [
    (1, "VP1", "KF729687.1"), (2, "VP2", "KF729688.1"), (3, "VP3", "KF729689.1"),
    (4, "VP4", "KF729690.1"), (5, "NSP1", "KF729691.1"), (6, "VP6", "KF729692.1"),
    (7, "VP7", "KF729693.1"), (8, "NSP2", "KF729694.1"), (9, "NSP3", "KF729695.1"),
    (10, "NSP4", "KF729696.1"), (11, "NSP5", "KF729697.1"),
]


def fetch(accs):
    """One efetch for all accessions; per-accession calls trip NCBI's rate limit."""
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={','.join(accs)}&rettype=fasta&retmode=text"
    with urllib.request.urlopen(url) as r:
        text = r.read().decode()
    out = {}
    for block in text.strip().split(">")[1:]:
        lines = block.splitlines()
        out[lines[0].split()[0]] = (lines[0], "".join(lines[1:]).upper())
    return out


fasta = fetch([acc for _, _, acc in SEGMENTS])
records = [(seg, name, acc, *fasta[acc]) for seg, name, acc in SEGMENTS]
os.makedirs("work", exist_ok=True)
with open(f"work/{STRAIN}.fa", "w") as fa, open(f"work/{STRAIN}.gtf", "w") as gtf, open("MANIFEST.tsv", "w") as man:
    man.write("dataset\tsegment\tprotein\tncbi_accession\tlength\n")
    for seg, name, acc, header, seq in records:
        fa.write(f">{header}\n" + "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\n")
        attrs = f'gene_id "{acc}"; transcript_id "{acc}";'
        for feat in ("transcript", "exon"):
            gtf.write(f"{acc}\tncbi\t{feat}\t1\t{len(seq)}\t.\t+\t.\t{attrs}\n")
        man.write(f"rnastruct00024\t{seg}\t{name}\t{acc}\t{len(seq)}\n")
        print(f"seg {seg} {name:4s} {acc} {len(seq)} nt", file=sys.stderr)
