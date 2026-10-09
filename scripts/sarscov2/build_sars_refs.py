"""Build the SARS-CoV-2 references: genome + every RNAcentral record that sits in it.

Leiden-0002 (MT510999.1) serves rnastruct00025/00028/00030, USA-WA1/2020 (MT576563.1)
rnastruct00031. GSE279203's WT (rnastruct00072) is Wuhan-Hu-1 base for base; its Alpha,
Beta, Delta and Omicron XBB (rnastruct00101-00104) are GISAID-only, so their genomes are read
from the authors' ShapeMapper .map files on GEO. RNAcentral's SARS-CoV-2 records are PDB constructs and ENA miRNAs;
those of >= 60 nt that align to the genome at >= 99 % identity become transcripts of their own
(dengue's sfRNA rule) so an element can be scored on both; where records overlap
(the 124-nt SL5 inside the 145-nt one, the FSE inside two ribosome-complex mRNAs) only the
best by identity then length is kept. The FASTA is the genome only; each record gets a
NCBI_GTF-style transcript/exon pair placed on the genome at its aligned coordinates. The 5' and
3' UTRs are added as <accession>_5UTR/_3UTR, bounded by Wuhan-Hu-1's ORF1ab start and ORF10 stop
(its RefSeq UTR bounds) found in each genome, since the GISAID-only genomes have no annotation.
"""
import gzip, json, os, sys, time, urllib.request
from Bio import Align

STRAINS = [
    # datasets, label, NCBI accession or GISAID id, GEO .map name for GISAID-only genomes
    ("rnastruct00025,rnastruct00028,rnastruct00030", "Leiden0002", "MT510999.1", None),
    ("rnastruct00031", "USAWA1", "MT576563.1", None),
    ("rnastruct00072", "WT", "NC_045512.2", None),
    ("rnastruct00101", "Alpha", "EPI_ISL_754083", "ALPHA"),
    ("rnastruct00102", "Beta", "EPI_ISL_1173248", "BETA"),
    ("rnastruct00103", "Delta", "EPI_ISL_2621925", "DELTA"),
    ("rnastruct00104", "OmicronXBB", "EPI_ISL_14917728", "XBB"),
]
GEO_MAP = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE279nnn/GSE279203/suppl/GSE279203_{}.map.gz"
# Illumina read-2 adapter left on the 5' end of the Delta consensus
ADAPTER = "GGAGTTCAGACGTGTGCTCTTCCGATCT"
TAXID, MIN_LEN, MIN_IDENT = 2697049, 60, 0.99
# Wuhan-Hu-1 (NC_045512.2, also the WT genome) ORF1ab first base and ORF10 last base; K-mer anchors
REF, ORF1AB_START, ORF10_END, K = "NC_045512.2", 266, 29674, 25

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


def fetch_geo(label, epi, name):
    rows = gzip.decompress(get(GEO_MAP.format(name))).decode().splitlines()
    seq = "".join(row.split("\t")[3] for row in rows).upper().replace("U", "T").removeprefix(ADAPTER)
    return f"{epi} Severe acute respiratory syndrome coronavirus 2 {label}, GISAID {epi} consensus from GEO GSE279203_{name}.map", seq


def rnacentral_records():
    d = json.loads(get(f"https://www.ebi.ac.uk/ebisearch/ws/rest/rnacentral?query=TAXONOMY:{TAXID}&fields=id,length,description&size=100&format=json"))
    for e in d["entries"]:
        if int(e["fields"]["length"][0]) < MIN_LEN:
            continue
        time.sleep(0.5)
        # rnacentral.org's cache serves its HTML home page to some uncompressed API requests
        req = urllib.request.Request(f"https://rnacentral.org/api/v1/rna/{e['id'].split('_')[0]}?format=json",
                                     headers={"Accept-Encoding": "gzip"})
        seq = json.loads(gzip.decompress(get(req)))["sequence"]
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


genomes = fetch_genomes([acc for _, _, acc, geo in STRAINS if not geo])
genomes.update({acc: fetch_geo(label, acc, geo) for _, label, acc, geo in STRAINS if geo})
candidates = list(rnacentral_records())
os.makedirs("work", exist_ok=True)
rows = ["datasets\tstrain\taccession\trecord\trecord_len\tidentity\tgenome_start\tgenome_end\tdescription"]
for datasets, label, acc, _ in STRAINS:
    header, genome = genomes[acc]
    gtf = [(acc, acc, 1, len(genome))]
    rows.append(f"{datasets}\t{label}\t{acc}\t{acc}\t{len(genome)}\t1.0000\t1\t{len(genome)}\tgenome")
    hits = []
    for urs, seq, desc in candidates:
        a = aligner.align(genome, seq)[0]
        ident = sum(x == y for x, y in zip(*a)) / len(seq)
        if ident >= MIN_IDENT:
            hits.append((ident, len(seq), urs, desc, a.aligned[0][0][0] + 1, a.aligned[0][-1][1]))
    kept = []
    for ident, n, urs, desc, start, end in sorted(hits, key=lambda h: (-h[0], -h[1])):
        if any(start <= k[5] and k[4] <= end for k in kept):
            continue
        kept.append((ident, n, urs, desc, start, end))
    for ident, n, urs, desc, start, end in sorted(kept, key=lambda k: k[4]):
        gtf.append((urs, acc, start, end))
        rows.append(f"{datasets}\t{label}\t{acc}\t{urs}\t{n}\t{ident:.4f}\t{start}\t{end}\t{desc}")
    ref = genomes[REF][1]
    cds_start = genome.index(ref[ORF1AB_START - 1:ORF1AB_START - 1 + K]) + 1
    cds_end = genome.index(ref[ORF10_END - K:ORF10_END]) + K
    gtf += [(f"{acc}_5UTR", acc, 1, cds_start - 1), (f"{acc}_3UTR", acc, cds_end + 1, len(genome))]
    write_fasta(f"work/SARS2_{label}.fa", [(header, genome)])
    write_gtf(f"work/SARS2_{label}.gtf", gtf)
    print("\n".join(rows[-(len(gtf) - 2):]), file=sys.stderr)
open("MANIFEST-genome.tsv", "w").write("\n".join(rows) + "\n")
