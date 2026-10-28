#!/usr/bin/env python3
"""Build release<N>/<org>/rfeval/rrna_common.{db,windows.tsv}, the method of the archived build_rrna_common.py.

The cytoplasmic rRNAs are found by cmsearch --cut_ga among the transcripts a run is keyed on
(the filtered GTF's rRNA records for eukaryotes, the transcriptome FASTA's records for bacteria
and archaea), one locus per molecule within 10 kb of the SSU, SS_cons projected without
pseudoknots, and only structures of >= 0.10 bp/nt kept.

    python3 scripts/rfeval_rrna_common.py <org> <release>    # from the repository root
"""
import glob, gzip, json, os, subprocess, sys, time, urllib.error, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from project_ss_cons import pair_columns, project, read_pfam_sto

IMG = "rnacentral/r2dt:latest"
CMS = os.path.join(HERE, "rrna_cm")
EUK = {"RF01960": "18S_rRNA", "RF00002": "5.8S_rRNA", "RF02543": "28S_rRNA", "RF00001": "5S_rRNA"}
BAC = {"RF00177": "16S_rRNA", "RF02541": "23S_rRNA", "RF00001": "5S_rRNA"}
ARC = {"RF01959": "16S_rRNA", "RF02540": "23S_rRNA", "RF00001": "5S_rRNA"}
# 5S genes outside the 45S rDNA unit are searched genome-wide ("unlinked")
ORGANISMS = {
    "arabidopsis_thaliana": {"species": "arabidopsis_thaliana", "families": {**EUK, "RF02543": "25S_rRNA"}, "unlinked": {"RF00001"}},
    "danio_rerio": {"species": "danio_rerio", "families": EUK, "unlinked": {"RF00001"}},
    "plasmodium_falciparum": {"species": "plasmodium_falciparum", "families": EUK, "unlinked": {"RF00001"}},
    "triticum_turgidum": {"species": "triticum_turgidum", "families": {**EUK, "RF02543": "26S_rRNA"}, "unlinked": {"RF00001"}},
    "oryza_sativa": {"species": "oryza_sativa", "families": {**EUK, "RF02543": "25S_rRNA"}, "unlinked": {"RF00001"}},
    "bacillus_subtilis": {"families": BAC},
    "listeria_monocytogenes": {"families": BAC},
    "pseudomonas_savastanoi": {"families": BAC},
    "methanosarcina_acetivorans": {"families": ARC},
}
MIN_BP_PER_NT, UNIT = 0.10, 10000


def post_json(url, body, tries=4):
    req = urllib.request.Request(url, json.dumps(body).encode(), {"Content-Type": "application/json", "Accept": "application/json"})
    for i in range(tries):
        try:
            return json.load(urllib.request.urlopen(req, timeout=300))
        except (urllib.error.URLError, TimeoutError):
            if i == tries - 1:
                raise
            time.sleep(5 * (i + 1))


def region_seqs(species, regions):
    """{region: sequence} from Ensembl REST; a batch that keeps failing is split to isolate the bad region."""
    try:
        return {r["query"]: r["seq"].upper() for r in post_json(f"https://rest.ensembl.org/sequence/region/{species}", {"regions": regions})}
    except urllib.error.HTTPError:
        if len(regions) == 1:
            raise
        mid = len(regions) // 2
        return {**region_seqs(species, regions[:mid]), **region_seqs(species, regions[mid:])}


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]; seqs[name] = []
        elif name:
            seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def gtf_rrna_exons(path):
    """rRNA-biotype transcripts >= 100 nt: {transcript_id: [(chrom, start, end, strand), ...]} in transcript order."""
    exons = {}
    for line in gzip.open(path, "rt"):
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon" or 'transcript_biotype "rRNA"' not in f[8]:
            continue
        tx = f[8].split('transcript_id "', 1)[1].split('"', 1)[0]
        exons.setdefault(tx, []).append((f[0], int(f[3]), int(f[4]), f[6]))
    for ex in exons.values():
        ex.sort(key=lambda e: e[1], reverse=ex[0][3] == "-")
    return {tx: ex for tx, ex in exons.items() if sum(e[2] - e[1] + 1 for e in ex) >= 100}


def candidates(org, orgdir, fa):
    """Write the candidate sequences; return {transcript: (chrom, start, end)} genome loci."""
    if "species" not in org:
        man = [l.rstrip("\n").split("\t") for l in open(os.path.join(orgdir, "MANIFEST-genome.tsv"))]
        col = {h: i for i, h in enumerate(man[0])}
        rows = [r for r in man[1:] if r[col["description"]] != "genome"]
        chrom = lambda r: r[col["replicon"]] if "replicon" in col else "chr"
        seqs = read_fasta(glob.glob(os.path.join(orgdir, "reference", "*.transcriptome.fa"))[0])
        with open(fa, "w") as fh:
            fh.writelines(f">{r[0]}\n{seqs[r[0]]}\n" for r in rows)
        return {r[0]: (chrom(r), int(r[col["genome_start"]]), int(r[col["genome_end"]])) for r in rows}
    exons = gtf_rrna_exons(glob.glob(os.path.join(orgdir, "reference", "*.filtered.gtf.gz"))[0])
    # Ensembl REST spells the same assembly the pipeline downloads
    regions = sorted({f"{c}:{s}..{e}:{1 if st == '+' else -1}" for ex in exons.values() for c, s, e, st in ex})
    got = {}
    for i in range(0, len(regions), 50):
        got.update(region_seqs(org["species"], regions[i:i + 50]))
    with open(fa, "w") as fh:
        for tx, ex in sorted(exons.items()):
            fh.write(f">{tx}\n" + "".join(got[f"{c}:{s}..{e}:{1 if st == '+' else -1}"] for c, s, e, st in ex) + "\n")
    return {tx: (ex[0][0], min(x[1] for x in ex), max(x[2] for x in ex)) for tx, ex in exons.items()}


def unpair_orphans(dot):
    """A pair whose partner fell outside the trimmed span is no longer a pair."""
    stack, fixed = [], list(dot)
    for i, ch in enumerate(dot):
        if ch == "(":
            stack.append(i)
        elif ch == ")":
            if stack: stack.pop()
            else: fixed[i] = "."
    for i in stack:
        fixed[i] = "."
    return "".join(fixed)


def main():
    name, release = sys.argv[1], sys.argv[2]
    org, orgdir = ORGANISMS[name], os.path.join(f"release{release}", name)
    work, out = os.path.abspath(os.path.join(orgdir, "work", "rrna_common")), os.path.join(orgdir, "rfeval")
    os.makedirs(work, exist_ok=True); os.makedirs(out, exist_ok=True)
    docker = lambda *cmd: subprocess.run(["docker", "run", "--rm", "-u", f"{os.getuid()}:{os.getgid()}", "-v", f"{work}:/w",
                                          "-v", f"{CMS}:/cm:ro", "-w", "/w", IMG, *cmd], check=True, capture_output=True)
    fa = os.path.join(work, "candidates.fa")
    loci = candidates(org, orgdir, fa)
    seqs = read_fasta(fa)
    print(f"{name}: {len(seqs)} candidate sequences")

    found = {}
    for acc, molecule in org["families"].items():
        # --toponly: only sense-strand hits are measured on a transcript
        tbl = os.path.join(work, f"{acc}.tbl")
        if not os.path.exists(tbl):
            docker("cmsearch", "--cpu", "8", "--cut_ga", "--nohmmonly", "--toponly", "--tblout", f"{acc}.tbl", f"/cm/{acc}.cm", "candidates.fa")
        hits = [l.split() for l in open(tbl) if not l.startswith("#")]
        found[acc] = sorted(((float(h[14]), len(seqs[h[0]]), h[0], int(h[7]), int(h[8])) for h in hits), key=lambda h: (-h[0], h[1], h[2]))
        for h in found[acc][:5]:
            print(f"  {molecule:10s} {h[2]:32s} {h[3]}-{h[4]} {h[0]:.1f} bits (transcript {h[1]} nt, {loci[h[2]][0]}:{loci[h[2]][1]})")
    near = lambda tx, a: (l := loci[tx])[0] == a[0] and max(l[1], a[1]) - min(l[2], a[2]) < UNIT
    ssu, *rest = org["families"]
    linked = [acc for acc in rest if acc not in org.get("unlinked", ())]
    # identical rRNA copies share one URS, so anchor on the SSU whose unit holds the most linked molecules, then by score
    anchor = max((loci[h[2]] for h in found[ssu]), key=lambda a: sum(any(near(h[2], a) for h in found[acc]) for acc in linked), default=None)

    db, windows, manifest = [], [], []
    for acc, molecule in org["families"].items():
        hits = found[acc]
        if acc != ssu and acc not in org.get("unlinked", ()):
            # the rest must come from the SSU's own rDNA unit/operon, not a pseudogene or an organelle
            hits = [h for h in hits if near(h[2], anchor)]
        elif acc == ssu:
            hits = [h for h in hits if loci[h[2]] == anchor]
        if not hits:
            print(f"  {molecule}: no GA hit in the rDNA unit"); continue
        # identical copies tie on score; the shortest transcript is the mature molecule itself
        score, _, tx, start, end = min(hits, key=lambda h: (-h[0], h[1], h[2]))
        open(os.path.join(work, f"{acc}.fa"), "w").write(f">{tx}\n{seqs[tx][start - 1:end]}\n")
        docker("cmalign", "--cpu", "8", "--outformat", "Pfam", "-o", f"{acc}.sto", f"/cm/{acc}.cm", f"{acc}.fa")
        aligned, ss_cons = read_pfam_sto(os.path.join(work, f"{acc}.sto"))
        seq, dot, lo, hi = project(aligned[tx], pair_columns(ss_cons))
        seq, dot = seq[lo:hi + 1], unpair_orphans(dot[lo:hi + 1])
        sid, bp = f"{tx}__{molecule}", dot.count("(")
        print(f"  -> {sid} {start + lo}-{start + hi}: {bp} bp over {len(seq)} nt ({bp / len(seq):.2f} bp/nt)")
        if bp / len(seq) < MIN_BP_PER_NT:
            print(f"  {molecule}: below {MIN_BP_PER_NT} bp/nt, dropped"); continue
        db.append(f">{sid}\n{seq}\n{dot}\n")
        windows.append(f"{tx}\t{start + lo}\t{start + hi}\t{sid}\n")
        manifest.append(f"{sid}\t{acc}\t{tx}\t{start + lo}\t{start + hi}\t{loci[tx][0]}:{loci[tx][1]}-{loci[tx][2]}\t{len(seq)}\t{bp}\t{score:.1f}\n")

    open(os.path.join(out, "rrna_common.db"), "w").writelines(db)
    open(os.path.join(out, "rrna_common.windows.tsv"), "w").writelines(["# ref_seq_id\tstart\tend\tstructure_id\n"] + windows)
    open(os.path.join(orgdir, "MANIFEST-rfeval.tsv"), "w").writelines(
        ["structure\trfam\tref_seq_id\tstart\tend\tlocus\tlength\tbase_pairs\tbits\n"] + manifest)


if __name__ == "__main__":
    main()
