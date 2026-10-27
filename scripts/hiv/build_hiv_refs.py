"""Build the HIV-1 NL4-3 reference (rnastruct00119): NCBI genome + every RNAcentral record that sits in it.

RNAcentral records of >= 60 nt that align to the genome at >= 99 % identity become transcripts of
their own; overlapping records keep only the best by identity, then length, then ID. The FASTA is the
genome only; the GTF has the genome, each record at its aligned span, and the 5' and 3' UTRs cut at
the first CDS start and last CDS end, as <accession>_5UTR/_3UTR. Records are read from work/records.fa (fetched once if absent).
"""
import gzip, json, os, re, sys, time, urllib.request
from Bio import Align

STRAINS = [
    # datasets, label, NCBI accession, last proviral base (AF324493.2 continues into cellular DNA and pUC vector)
    ("rnastruct00119", "NL43", "AF324493.2", 9709),
]
TAXID, MIN_LEN, MIN_IDENT, MIN_UTR, PREFIX = 11676, 60, 0.99, 20, "HIV1"
CACHE = "work/records.fa"
os.makedirs("work", exist_ok=True)

aligner = Align.PairwiseAligner(mode="local", match_score=1, mismatch_score=-1,
                                open_gap_score=-2, extend_gap_score=-0.5)


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


def fetch_genomes(accs):
    text = get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={','.join(accs)}&rettype=fasta&retmode=text").decode()
    out = {}
    for block in text.strip().split(">")[1:]:
        lines = block.splitlines()
        out[lines[0].split()[0]] = (lines[0], "".join(lines[1:]).upper())
    return out


def fetch_cds(acc):
    gb = get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={acc}&rettype=gb&retmode=text").decode()
    # first CDS start to last CDS end over the plus-strand CDS features (not a plasmid's complement(bla))
    locs = re.findall(r"^ {5}CDS +(\S[^\n]*(?:\n {21}[^/ ][^\n]*)*)", gb, re.M)
    nums = [int(n) for loc in locs if not loc.startswith("complement") for n in re.findall(r"\d+", loc)]
    return min(nums), max(nums)


def read_records(path):
    out, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name, desc = line[1:].rstrip("\n").split(" ", 1); out[name] = ["", desc]
        else:
            out[name][0] += line.strip()
    return {k: tuple(v) for k, v in out.items()}


def rnacentral_records():
    """{urs_taxid: (seq, description)} of >= MIN_LEN nt, cached in records.fa."""
    if os.path.exists(CACHE):
        return read_records(CACHE)
    d = json.loads(get(f"https://www.ebi.ac.uk/ebisearch/ws/rest/rnacentral?query=TAXONOMY:{TAXID}&fields=id,length,description&size=100&format=json"))
    out = {}
    for e in d["entries"]:
        if int(e["fields"]["length"][0]) < MIN_LEN:
            continue
        time.sleep(1)
        # rnacentral.org's cache serves its HTML home page to some uncompressed API requests
        req = urllib.request.Request(f"https://rnacentral.org/api/v1/rna/{e['id'].split('_')[0]}?format=json",
                                     headers={"Accept-Encoding": "gzip"})
        raw = get(req)
        seq = json.loads(gzip.decompress(raw) if raw[:2] == b"\x1f\x8b" else raw)["sequence"]
        out[e["id"]] = (seq.upper().replace("U", "T"), e["fields"]["description"][0])
    with open(CACHE, "w") as fh:
        for rid, (seq, desc) in out.items():
            fh.write(f">{rid} {desc}\n{seq}\n")
    return out


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


genomes = fetch_genomes([acc for _, _, acc, _ in STRAINS])
candidates = [(urs, seq, desc) for urs, (seq, desc) in rnacentral_records().items() if len(seq) >= MIN_LEN]
rows = ["datasets\tstrain\taccession\trecord\trecord_len\tidentity\tgenome_start\tgenome_end\tdescription"]
for datasets, label, acc, provirus_end in STRAINS:
    header, genome = genomes[acc]
    genome = genome[:provirus_end]
    gtf, first = [(acc, acc, 1, len(genome))], len(rows)
    rows.append(f"{datasets}\t{label}\t{acc}\t{acc}\t{len(genome)}\t1.0000\t1\t{len(genome)}\tgenome")
    hits = []
    for urs, seq, desc in candidates:
        a = aligner.align(genome, seq)[0]
        ident = sum(x == y for x, y in zip(*a)) / len(seq)
        if ident >= MIN_IDENT:
            hits.append((ident, len(seq), urs, seq, desc, a.aligned[0][0][0] + 1, a.aligned[0][-1][1]))
    kept = []
    for ident, n, urs, seq, desc, start, end in sorted(hits, key=lambda h: (-h[0], -h[1], h[2])):
        if any(start <= k[6] and k[5] <= end for k in kept):
            continue
        kept.append((ident, n, urs, seq, desc, start, end))
    for ident, n, urs, seq, desc, start, end in sorted(kept, key=lambda k: (k[5], k[2])):
        gtf.append((urs, acc, start, end))
        rows.append(f"{datasets}\t{label}\t{acc}\t{urs}\t{n}\t{ident:.4f}\t{start}\t{end}\t{desc}")
    cds_start, cds_end = fetch_cds(acc)
    # a UTR shorter than MIN_UTR is not worth a transcript of its own
    gtf += [(f"{acc}_{n}", acc, s, e) for n, s, e in (("5UTR", 1, cds_start - 1), ("3UTR", cds_end + 1, len(genome))) if e - s + 1 >= MIN_UTR]
    write_fasta(f"work/{PREFIX}_{label}.fa", [(header, genome)])
    write_gtf(f"work/{PREFIX}_{label}.gtf", gtf)
    print("\n".join(rows[first:]), file=sys.stderr)
    print(f"CDS {cds_start}..{cds_end}", file=sys.stderr)
open("MANIFEST-genome.tsv", "w").write("\n".join(rows) + "\n")
