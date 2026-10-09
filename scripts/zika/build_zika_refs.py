"""Build the three Zika references: NCBI genome + every RNAcentral record that sits in it.

MR766 (AY632535.2) serves rnastruct00009/00051, BeH815744 (KU365780.1) rnastruct00052,
H/PF/2013 (KJ776791.2) rnastruct00053. RNAcentral's Zika records are partial sfRNAs (ENA)
and PDB xrRNA constructs; those of >= 60 nt that align to the genome at >= 99 % identity
become transcripts of their own (dengue's sfRNA rule) so an element can be scored on both;
overlapping records keep only the best by identity then length. The FASTA is the genome only;
the GTF has a NCBI_GTF-style transcript/exon pair for the genome, each record at its aligned
span, and the 5' and 3' UTRs either side of the GenBank CDS (<accession>_5UTR/_3UTR).
"""
import json, os, re, sys, time, urllib.request
from Bio import Align

STRAINS = [
    # datasets, label, NCBI accession
    ("rnastruct00009,rnastruct00051", "MR766", "AY632535.2"),
    ("rnastruct00052", "BeH815744", "KU365780.1"),
    ("rnastruct00053", "HPF2013", "KJ776791.2"),
]
TAXID, MIN_LEN, MIN_IDENT = 64320, 60, 0.99

aligner = Align.PairwiseAligner(mode="local", match_score=1, mismatch_score=-1,
                                open_gap_score=-2, extend_gap_score=-0.5)


def get(url, tries=6):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code == 429:
                time.sleep(3 * (i + 1)); continue
            raise
    raise RuntimeError(f"gave up on {url}")


def fetch_genomes(accs):
    text = get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={','.join(accs)}&rettype=fasta&retmode=text").decode()
    out = {}
    for block in text.strip().split(">")[1:]:
        lines = block.splitlines()
        out[lines[0].split()[0]] = (lines[0], "".join(lines[1:]).upper())
    return out


def fetch_cds(acc):
    gb = get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={acc}&rettype=gb&retmode=text").decode()
    start, end = re.search(r"^ +CDS +(\d+)\.\.(\d+)$", gb, re.M).groups()
    return int(start), int(end)


def rnacentral_records():
    d = json.loads(get(f"https://www.ebi.ac.uk/ebisearch/ws/rest/rnacentral?query=TAXONOMY:{TAXID}&fields=id,length,description&size=100&format=json"))
    for e in d["entries"]:
        if int(e["fields"]["length"][0]) < MIN_LEN:
            continue
        time.sleep(0.5)
        seq = json.loads(get(f"https://rnacentral.org/api/v1/rna/{e['id'].split('_')[0]}?format=json"))["sequence"]
        yield e["id"], seq.upper().replace("U", "T"), e["fields"]["description"][0]


def write_fasta(path, records):
    with open(path, "w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n" + "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\n")


def write_gtf(path, records):
    with open(path, "w") as fh:
        for rid, seqname, start, end in records:
            attrs = f'gene_id "{rid}"; transcript_id "{rid}";'
            for feat in ("transcript", "exon"):
                fh.write(f"{seqname}\tncbi\t{feat}\t{start}\t{end}\t.\t+\t.\t{attrs}\n")


genomes = fetch_genomes([acc for _, _, acc in STRAINS])
candidates = list(rnacentral_records())
os.makedirs("work", exist_ok=True)
rows = ["datasets\tstrain\taccession\trecord\trecord_len\tidentity\tgenome_start\tgenome_end\tdescription"]
for datasets, label, acc in STRAINS:
    header, genome = genomes[acc]
    gtf = [(acc, acc, 1, len(genome))]
    rows.append(f"{datasets}\t{label}\t{acc}\t{acc}\t{len(genome)}\t1.0000\t1\t{len(genome)}\tgenome")
    hits = []
    for urs, seq, desc in candidates:
        a = aligner.align(genome, seq)[0]
        ident = sum(x == y for x, y in zip(*a)) / len(seq)
        if ident >= MIN_IDENT:
            hits.append((ident, len(seq), urs, seq, desc, a.aligned[0][0][0] + 1, a.aligned[0][-1][1]))
    kept = []
    for ident, n, urs, seq, desc, start, end in sorted(hits, key=lambda h: (-h[0], -h[1])):
        if any(start <= k[6] and k[5] <= end for k in kept):
            continue
        kept.append((ident, n, urs, seq, desc, start, end))
    for ident, n, urs, seq, desc, start, end in sorted(kept, key=lambda k: k[5]):
        gtf.append((urs, acc, start, end))
        rows.append(f"{datasets}\t{label}\t{acc}\t{urs}\t{n}\t{ident:.4f}\t{start}\t{end}\t{desc}")
    cds_start, cds_end = fetch_cds(acc)
    gtf += [(f"{acc}_5UTR", acc, 1, cds_start - 1), (f"{acc}_3UTR", acc, cds_end + 1, len(genome))]
    write_fasta(f"work/ZIKV_{label}.fa", [(header, genome)])
    write_gtf(f"work/ZIKV_{label}.gtf", gtf)
    print("\n".join(rows[-(len(gtf) - 2):]), file=sys.stderr)
open("MANIFEST-genome.tsv", "w").write("\n".join(rows) + "\n")
