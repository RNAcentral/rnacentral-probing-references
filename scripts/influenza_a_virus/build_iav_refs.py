"""Build the influenza A references: the eight NCBI segment records of each strain.

A/WSN/1933(H1N1) serves rnastruct00013; A/Puerto Rico/8/1934(H1N1) serves rnastruct00014, 00134
and 00135. Same accessions as conf/viral_genomes.config, pinned locally so the run does not depend
on an NCBI fetch. RNAcentral holds no full-length influenza A transcript (only PDB
promoter fragments), so unlike dengue there is no second record to add.
"""
import os, sys, time, urllib.request

STRAINS = [
    # stem, datasets, [(segment, protein, accession)]
    ("IAV_WSN1933", "rnastruct00013", [
        (1, "PB2", "CY034139.1"), (2, "PB1", "CY034138.1"), (3, "PA", "CY034137.1"),
        (4, "HA", "CY034132.1"), (5, "NP", "CY034135.1"), (6, "NA", "CY034134.1"),
        (7, "M", "CY034133.1"), (8, "NS", "CY034136.1")]),
    ("IAV_PR8", "rnastruct00014,rnastruct00134,rnastruct00135", [
        (1, "PB2", "NC_002023.1"), (2, "PB1", "NC_002021.1"), (3, "PA", "NC_002022.1"),
        (4, "HA", "NC_002017.1"), (5, "NP", "NC_002019.1"), (6, "NA", "NC_002018.1"),
        (7, "M", "NC_002016.1"), (8, "NS", "NC_002020.1")]),
]


def get(url, tries=6):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code == 429 or e.code >= 500:
                time.sleep(3 * (i + 1)); continue
            raise
        except (urllib.error.URLError, TimeoutError):
            time.sleep(3 * (i + 1))
    raise RuntimeError(f"gave up on {url}")


def fetch(accs):
    """One efetch per strain; per-accession calls trip NCBI's rate limit."""
    text = get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={','.join(accs)}&rettype=fasta&retmode=text").decode()
    out = {}
    for block in text.strip().split(">")[1:]:
        lines = block.splitlines()
        out[lines[0].split()[0]] = (lines[0], "".join(lines[1:]).upper())
    return out


os.makedirs("work", exist_ok=True)
with open("MANIFEST.tsv", "w") as man:
    man.write("dataset\tsegment\tprotein\tncbi_accession\tlength\n")
    for stem, datasets, segments in STRAINS:
        fasta = fetch([acc for _, _, acc in segments])
        with open(f"work/{stem}.fa", "w") as fa, open(f"work/{stem}.gtf", "w") as gtf:
            for seg, name, acc in segments:
                header, seq = fasta[acc]
                fa.write(f">{header}\n" + "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\n")
                attrs = f'gene_id "{acc}"; transcript_id "{acc}";'
                for feat in ("transcript", "exon"):
                    gtf.write(f"{acc}\tncbi\t{feat}\t1\t{len(seq)}\t.\t+\t.\t{attrs}\n")
                man.write(f"{datasets}\t{seg}\t{name}\t{acc}\t{len(seq)}\n")
                print(f"{stem} seg {seg} {name:4s} {acc} {len(seq)} nt", file=sys.stderr)
        time.sleep(1)
