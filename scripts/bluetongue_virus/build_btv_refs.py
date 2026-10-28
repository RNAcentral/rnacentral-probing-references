"""Build the bluetongue virus serotype 1 (South Africa) reference for rnastruct00077: the ten NCBI segment records.

The ten dsRNA segments were probed as T7 transcripts, so each is a whole contig and transcript
(no UTR split), as for influenza and rotavirus. FJ969724.1 is the NS1 gene (segment 5) although
GenBank labels it "segment 6"; FJ969723.1 is VP5 (segment 6). RNAcentral holds no full-length BTV transcript.
"""
import os, sys, time, urllib.request

STRAIN = "BTV1"
SEGMENTS = [
    (1, "VP1", "FJ969719.1"), (2, "VP2", "FJ969720.1"), (3, "VP3", "FJ969721.1"),
    (4, "VP4", "FJ969722.1"), (5, "NS1", "FJ969724.1"), (6, "VP5", "FJ969723.1"),
    (7, "VP7", "FJ969725.1"), (8, "NS2", "FJ969726.1"), (9, "VP6", "FJ969727.1"),
    (10, "NS3", "FJ969728.1"),
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
    """One efetch for all accessions; per-accession calls trip NCBI's rate limit."""
    text = get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={','.join(accs)}&rettype=fasta&retmode=text").decode()
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
        man.write(f"rnastruct00077\t{seg}\t{name}\t{acc}\t{len(seq)}\n")
        print(f"seg {seg} {name:4s} {acc} {len(seq)} nt", file=sys.stderr)
