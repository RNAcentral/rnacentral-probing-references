"""Build per-serotype dengue references: NCBI EDEN genome + best-matching RNAcentral sfRNA.

For each serotype, every full-length sfRNA1 (or sfRNA2) entry in RNAcentral for that taxid is
locally aligned to the EDEN genome's 3' end and the closest is placed on the genome at its
aligned coordinates. The FASTA is the genome only; the GTF has a NCBI_GTF-style transcript/exon
pair for the genome, the sfRNA and the 5' and 3' UTRs either side of the GenBank CDS
(<accession>_5UTR/_3UTR: a trial of folding UTRs as their own records before RNAcentral has them).
"""
import json, os, re, sys, time, urllib.request, urllib.parse
from Bio import Align

SEROTYPES = [
    # dataset, serotype, strain, NCBI accession, RNAcentral taxid
    ("rnastruct00046", "DENV1", "EDEN2402", "EU081230.1", 11053),
    ("rnastruct00047", "DENV2", "EDEN3295", "EU081177.1", 11060),
    ("rnastruct00048", "DENV3", "EDEN863",  "EU081190.1", 11069),
    ("rnastruct00049", "DENV4", "EDEN2270", "GQ398256.1", 11070),
]

aligner = Align.PairwiseAligner(mode="local", match_score=1, mismatch_score=-1,
                                open_gap_score=-2, extend_gap_score=-0.5)

def get(url, tries=6):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code == 429:
                time.sleep(3 * (i + 1))
                continue
            raise
    raise RuntimeError(f"gave up on {url}")

def fetch_genome(acc):
    raw = get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={acc}&rettype=fasta&retmode=text").decode()
    lines = raw.strip().splitlines()
    return lines[0], "".join(lines[1:]).upper()

def sfrna_candidates(taxid):
    # DENV-4 has no sfRNA1 entries in RNAcentral; its full-length sfRNA is annotated as sfRNA2
    for name in ("sfRNA1", "sfRNA2"):
        q = urllib.parse.quote(f"TAXONOMY:{taxid} AND NOT partial AND {name}")
        ids, start = [], 0
        while True:
            d = json.loads(get(f"https://www.ebi.ac.uk/ebisearch/ws/rest/rnacentral?query={q}&fields=id&size=100&start={start}&format=json"))
            ids += [e["id"] for e in d["entries"]]
            start += 100
            if start >= d["hitCount"]:
                break
        if ids:
            return ids
    raise RuntimeError(f"no full-length sfRNA for taxid {taxid}")

def fetch_cds(acc):
    gb = get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={acc}&rettype=gb&retmode=text").decode()
    start, end = re.search(r"^ +CDS +(\d+)\.\.(\d+)$", gb, re.M).groups()
    return int(start), int(end)

def best_match(genome, taxid):
    tail = genome[-1500:]
    best = None
    for urs_tax in sfrna_candidates(taxid):
        time.sleep(0.7)
        urs = urs_tax.split("_")[0]
        seq = json.loads(get(f"https://rnacentral.org/api/v1/rna/{urs}?format=json"))["sequence"].upper().replace("U", "T")
        a = aligner.align(tail, seq)[0]
        ident = sum(1 for x, y in zip(*a) if x == y) / len(seq)
        if best is None or ident > best[0]:
            best = (ident, urs_tax, seq)
    return best

def write_fasta(path, records):
    with open(path, "w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")

def write_gtf(path, records):
    with open(path, "w") as fh:
        for rid, seqname, start, end in records:
            attrs = f'gene_id "{rid}"; transcript_id "{rid}";'
            for feat in ("transcript", "exon"):
                fh.write(f"{seqname}\tncbi\t{feat}\t{start}\t{end}\t.\t+\t.\t{attrs}\n")

os.makedirs("work", exist_ok=True)
summary = ["dataset\tserotype\tstrain\tncbi_accession\tgenome_len\trnacentral_id\tsfrna_len\tidentity\tgenome_coords"]
for dataset, serotype, strain, acc, taxid in SEROTYPES:
    header, genome = fetch_genome(acc)
    ident, urs_tax, seq = best_match(genome, taxid)
    a = aligner.align(genome, seq)[0]
    coords = f"{a.aligned[0][0][0] + 1}-{a.aligned[0][-1][1]}"
    sfrna_place = (acc, a.aligned[0][0][0] + 1, a.aligned[0][-1][1])
    cds_start, cds_end = fetch_cds(acc)
    utrs = [(f"{acc}_5UTR", 1, cds_start - 1), (f"{acc}_3UTR", cds_end + 1, len(genome))]
    write_fasta(f"work/{serotype}_{strain}.fa", [(header[1:], genome)])
    write_gtf(f"work/{serotype}_{strain}.gtf", [(acc, acc, 1, len(genome)), (urs_tax, *sfrna_place)]
              + [(rid, acc, start, end) for rid, start, end in utrs])
    summary.append(f"{dataset}\t{serotype}\t{strain}\t{acc}\t{len(genome)}\t{urs_tax}\t{len(seq)}\t{ident:.4f}\t{coords}")
    print(summary[-1], file=sys.stderr)

open("MANIFEST-genome.tsv", "w").write("\n".join(summary) + "\n")
