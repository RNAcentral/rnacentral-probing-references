"""Build the A/WSN/1933(H1N1) reference for rnastruct00013: the eight NCBI segment records.

Same accessions as conf/viral_genomes.config, pinned locally so the run does not depend
on an NCBI fetch. RNAcentral holds no full-length influenza A transcript (only PDB
promoter fragments), so unlike dengue there is no second record to add.
"""
import os, sys, urllib.request

STRAIN = "IAV_WSN1933"
SEGMENTS = [
    (1, "PB2", "CY034139.1"), (2, "PB1", "CY034138.1"), (3, "PA", "CY034137.1"),
    (4, "HA", "CY034132.1"), (5, "NP", "CY034135.1"), (6, "NA", "CY034134.1"),
    (7, "M", "CY034133.1"), (8, "NS", "CY034136.1"),
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
        man.write(f"rnastruct00013\t{seg}\t{name}\t{acc}\t{len(seq)}\n")
        print(f"seg {seg} {name:4s} {acc} {len(seq)} nt", file=sys.stderr)
