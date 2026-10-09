"""Build the E. coli reference: K-12 MG1655 genome + every RNAcentral record that sits in it.

Same rule as the viral builds: RNAcentral K-12 records (taxids 511145 MG1655 and 83333
K-12) of >= 60 nt that align to NC_000913.3 at >= 99 % identity are placed on it as their own
records; where records overlap on the same strand only the best by identity then length
is kept. The FASTA is the genome; the GTF places each kept record on it as an NCBI_GTF-style
transcript/exon pair. The genome itself is not a GTF transcript: a 4.6 Mb rf-fold never finishes.
"""
import gzip, json, os, sys, time, urllib.request
from Bio import Align

ACC, TAXIDS, MIN_LEN, MIN_IDENT, SEED = "NC_000913.3", (511145, 83333), 60, 0.99, 24
CACHE = "work/records.fa"
os.makedirs("work", exist_ok=True)

aligner = Align.PairwiseAligner(mode="local", match_score=1, mismatch_score=-1,
                                open_gap_score=-2, extend_gap_score=-0.5)


def get(url, tries=6):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url) as r:
                return r.read()
        except urllib.error.HTTPError as e:
            if e.code == 429 or e.code >= 500:
                time.sleep(3 * (i + 1)); continue
            raise
    raise RuntimeError(f"gave up on {url}")


def read_records(path):
    out, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name, desc = line[1:].rstrip("\n").split(" ", 1); out[name] = ["", desc]
        else:
            out[name][0] += line.strip()
    return {k: tuple(v) for k, v in out.items()}


def rnacentral_records():
    """{urs_taxid: (seq, description)}, one entry per URS, cached in records.fa."""
    if os.path.exists(CACHE):
        return read_records(CACHE)
    entries, seen = [], set()
    for taxid in TAXIDS:
        for start in range(0, 10000, 100):
            d = json.loads(get(f"https://www.ebi.ac.uk/ebisearch/ws/rest/rnacentral?query=TAXONOMY:{taxid}"
                               f"&fields=length,description&size=100&start={start}&format=json"))
            for e in d["entries"]:
                urs = e["id"].split("_")[0]
                if urs not in seen and int(e["fields"]["length"][0]) >= MIN_LEN:
                    seen.add(urs); entries.append((e["id"], e["fields"]["description"][0]))
            if start + 100 >= d["hitCount"]:
                break
    # appended per record so a failed run resumes; renamed to records.fa once complete
    partial = CACHE + ".partial"
    out = read_records(partial) if os.path.exists(partial) else {}
    fh = open(partial, "a")
    for rid, desc in entries:
        if rid in out:
            continue
        # rnacentral.org throttles by serving its HTML home page instead of JSON after ~60
        # quick requests; pace at 1/s and back off, alternating gzip and plain requests
        for i in range(10):
            time.sleep(1 if i == 0 else 10 * i)
            req = urllib.request.Request(f"https://rnacentral.org/api/v1/rna/{rid.split('_')[0]}?format=json",
                                         headers={"Accept-Encoding": "gzip"} if i % 2 == 0 else {})
            raw = get(req)
            try:
                seq = json.loads(gzip.decompress(raw) if raw[:2] == b"\x1f\x8b" else raw)["sequence"]
                break
            except ValueError:
                continue
        else:
            raise RuntimeError(f"no JSON for {rid}")
        out[rid] = (seq.upper().replace("U", "T"), desc)
        fh.write(f">{rid} {desc}\n{out[rid][0]}\n"); fh.flush()
    fh.close()
    os.rename(partial, CACHE)
    return out


def revcomp(s):
    return s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


def locate(genome, seq):
    """Best (identity, strand, start, end) of seq on genome: an exact seed, then local alignment around it."""
    best, tried = (0.0, None, 0, 0), set()
    for strand, q in (("+", seq), ("-", revcomp(seq))):
        for off in (0, len(q) // 4, len(q) // 2, 3 * len(q) // 4, len(q) - SEED):
            # every occurrence: rRNA operons, tRNAs and sRNA repeats are multi-copy
            p = genome.find(q[off:off + SEED])
            while p >= 0:
                if (strand, p - off) not in tried:
                    tried.add((strand, p - off))
                    lo = max(0, p - off - 50)
                    a = aligner.align(genome[lo:p - off + len(q) + 50], q)[0]
                    ident = sum(x == y for x, y in zip(*a)) / len(q)
                    if ident > best[0]:
                        best = (ident, strand, lo + a.aligned[0][0][0] + 1, lo + a.aligned[0][-1][1])
                p = genome.find(q[off:off + SEED], p + 1)
    return best


text = get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={ACC}&rettype=fasta&retmode=text").decode()
header, genome = text.splitlines()[0][1:], "".join(text.splitlines()[1:]).upper()
candidates = rnacentral_records()
hits = []
for rid, (seq, desc) in candidates.items():
    ident, strand, start, end = locate(genome, seq)
    if ident >= MIN_IDENT:
        hits.append((ident, len(seq), rid, seq, desc, strand, start, end))
kept = []
for h in sorted(hits, key=lambda h: (-h[0], -h[1])):
    if any(h[5] == k[5] and h[6] <= k[7] and k[6] <= h[7] for k in kept):
        continue
    kept.append(h)
kept.sort(key=lambda k: k[6])

with open("work/ECOLI_MG1655.fa", "w") as fh:
    fh.write(f">{header}\n" + "\n".join(genome[i:i + 60] for i in range(0, len(genome), 60)) + "\n")
with open("work/ECOLI_MG1655.gtf", "w") as fh:
    for _, _, rid, _, _, strand, start, end in kept:
        attrs = f'gene_id "{rid}"; transcript_id "{rid}";'
        for feat in ("transcript", "exon"):
            fh.write(f"{ACC}\tncbi\t{feat}\t{start}\t{end}\t.\t{strand}\t.\t{attrs}\n")
with open("MANIFEST-genome.tsv", "w") as fh:
    fh.write("record\trecord_len\tidentity\tstrand\tgenome_start\tgenome_end\tdescription\n")
    fh.write(f"{ACC}\t{len(genome)}\t1.0000\t+\t1\t{len(genome)}\tgenome\n")
    for ident, n, rid, _, desc, strand, start, end in kept:
        fh.write(f"{rid}\t{n}\t{ident:.4f}\t{strand}\t{start}\t{end}\t{desc}\n")
print(f"{len(candidates)} RNAcentral records >= {MIN_LEN} nt, {len(hits)} at >= {MIN_IDENT:.0%} identity, "
      f"{len(kept)} kept after overlap", file=sys.stderr)
