#!/usr/bin/env python3
"""Build release<N>/<org>/rfeval/<set>.{db,windows.tsv} for viral references from Rfam (skill section E).

Each genome or segment in reference/*.transcriptome.fa is scanned with cmscan --cut_ga (clan
competition on, sense strand only) against the full Rfam, or against the organism's own family list
where an earlier set fixed one (dengue, influenza: their new strains are appended to that set). Each
hit is cmaligned back to its model with the seed structure mapped in, the SS_cons projected onto it
and trimmed to the model-covered span; pseudoknots are dropped unless the family is the pseudoknot,
hairpins with loops under 3 nt are opened, and structures under 0.10 bp/nt are dropped. Every
element is windowed on its genome (_ncbi) and on any RNAcentral record that contains it (_rnac).

    python3 scripts/rfeval_viral_elements.py <org> <release> <Rfam dir>    # from the repository root
"""
import glob, os, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from project_ss_cons import pair_columns, project, read_pfam_sto

IMG = "rnacentral/r2dt:latest"
MIN_BP_PER_NT, MIN_LOOP = 0.10, 3
ORGANISMS = {
    "yellow_fever_virus": {}, "west_nile_virus": {}, "murine_norovirus": {}, "porcine_epidemic_diarrhea_virus": {}, "human_bocavirus": {}, "bluetongue_virus": {},
    # the provirus runs LTR to LTR; the RNA genome is R to R (AF324493.2 454-9626)
    "hiv": {"rna_span": {"AF324493.2": (454, 9626)}},
    # family -> element labels in genome order, as scripts/dengue_virus/build_rfeval_refs.py
    "dengue_virus": {"set": "dengue_elements", "stems": ["DENV1_WP75", "DENV2_16681", "DENV4_H241", "DENV4_TVP360"],
               "elements": {"RF02340": ["SLA"], "RF03546": ["5UTR"], "RF00617": ["cHP"], "RF00525": ["DB1", "DB2"], "RF00185": ["3SL"]}},
    "influenza_a_virus": {"set": "iav_elements", "stems": ["IAV_PR8"], "elements": {"RF01099": ["NS_PK"]}, "suffix": False,
                  "keep_pk": {"RF01099"}},
}


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]; seqs[name] = []
        elif name:
            seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def pk_columns(ss_cons):
    """Pseudoknot pairs from WUSS letters: an uppercase run opens, its lowercase run closes."""
    opens, pairs = {}, {}
    for col, ch in enumerate(ss_cons):
        if ch.isupper():
            opens.setdefault(ch, []).append(col)
        elif ch.islower():
            i = opens[ch.upper()].pop()
            pairs[i], pairs[col] = col, i
    return pairs


def clip(dot):
    """Unpair any bracket whose partner was cut off with the flank."""
    out, stacks, close = list(dot), {"(": [], "[": []}, {")": "(", "]": "["}
    for i, c in enumerate(dot):
        if c in stacks:
            stacks[c].append(i)
        elif c in close:
            if stacks[close[c]]:
                stacks[close[c]].pop()
            else:
                out[i] = "."
    for st in stacks.values():
        for i in st:
            out[i] = "."
    return "".join(out)


def open_short_loops(dot):
    """Unpair the closing pair of any hairpin whose loop is under MIN_LOOP nt, until none is left."""
    out, changed = list(dot), True
    while changed:
        changed, stack = False, []
        for i, c in enumerate(out):
            if c == "(":
                stack.append(i)
            elif c == ")":
                j = stack.pop()
                if all(x == "." for x in out[j + 1:i]) and i - j - 1 < MIN_LOOP:
                    out[j] = out[i] = "."; changed = True
    return "".join(out)


def read_tbl(path):
    """cmscan --fmt 2 rows: family accession and name, sequence, span, E-value, included."""
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.split()
        yield dict(name=f[1], fam=f[2], seq=f[3], start=int(f[9]), end=int(f[10]), strand=f[11], evalue=float(f[17]), inc=f[18] == "!")


def gtf_records(path):
    """{record: (contig, start, end)} for the RNAcentral records the GTF places on a genome."""
    out = {}
    for line in open(path):
        f = line.split("\t")
        if len(f) > 8 and f[2] == "transcript":
            tid = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
            if tid.startswith("URS"):
                out[tid] = (f[0], int(f[3]), int(f[4]))
    return out


def main():
    name, release, rfam = sys.argv[1], sys.argv[2], os.path.abspath(sys.argv[3])
    org, orgdir = ORGANISMS[name], os.path.join(f"release{release}", name)
    work, out = os.path.abspath(os.path.join(orgdir, "work", "elements")), os.path.join(orgdir, "rfeval")
    os.makedirs(work, exist_ok=True)
    docker = lambda *cmd: subprocess.run(["docker", "run", "--rm", "-u", f"{os.getuid()}:{os.getgid()}", "-v", f"{work}:/w",
                                                "-v", f"{rfam}:/rfam:ro", "-w", "/w", IMG, *cmd], check=True, capture_output=True)
    elements, suffix, keep_pk = org.get("elements"), org.get("suffix", True), org.get("keep_pk", set())
    cm = "/rfam/Rfam.cm"
    if elements:
        open(os.path.join(work, "set.keys"), "w").write("\n".join(elements) + "\n")
        open(os.path.join(work, "set.cm"), "wb").write(docker("cmfetch", "-f", "/rfam/Rfam.cm", "set.keys").stdout)
        docker("cmpress", "-F", "set.cm")
        cm = "set.cm"
    stems = org.get("stems") or sorted(os.path.basename(p).removesuffix(".transcriptome.fa") for p in glob.glob(os.path.join(orgdir, "reference", "*.transcriptome.fa")))

    db, windows, manifest = [], [], []
    for stem in stems:
        seqs = read_fasta(os.path.join(orgdir, "reference", f"{stem}.transcriptome.fa"))
        records = gtf_records(os.path.join(orgdir, "work", f"{stem}.gtf"))
        genomes = {k: v for k, v in seqs.items() if not k.startswith("URS") and not re.search(r"_[35]UTR$", k)}
        with open(os.path.join(work, f"{stem}.fa"), "w") as fh:
            fh.writelines(f">{k}\n{v}\n" for k, v in genomes.items())
        # clan competition only applies to the full Rfam; a curated family set has no competing clan members
        clans = [] if elements else ["--clanin", "/rfam/Rfam.clanin", "--oskip"]
        docker("cmscan", "--cpu", "8", "--cut_ga", "--nohmmonly", "--toponly", "--fmt", "2", *clans,
               "--tblout", f"{stem}.tbl", "-o", f"{stem}.log", cm, f"{stem}.fa")
        span = org.get("rna_span", {})
        hits = [h for h in read_tbl(os.path.join(work, f"{stem}.tbl")) if h["inc"] and h["strand"] == "+"
                and (h["seq"] not in span or span[h["seq"]][0] <= h["start"] and h["end"] <= span[h["seq"]][1])]
        if elements:
            # one hit per element label; families with two copies take the two best by E
            chosen = []
            for fam, labels in elements.items():
                best = sorted((h for h in hits if h["fam"] == fam), key=lambda h: h["evalue"])[:len(labels)]
                chosen += [dict(h, element=l) for l, h in zip(labels, sorted(best, key=lambda h: (h["seq"], h["start"])))]
        else:
            counts, seen, chosen = {}, {}, []
            for h in hits:
                counts[(h["seq"], h["name"])] = counts.get((h["seq"], h["name"]), 0) + 1
            for h in sorted(hits, key=lambda h: (h["seq"], h["start"])):
                k = (h["seq"], h["name"]); seen[k] = seen.get(k, 0) + 1
                seg = f"{h['seq'].split('.')[0]}_" if len(genomes) > 1 else ""
                chosen.append(dict(h, element=seg + h["name"] + (str(seen[k]) if counts[k] > 1 else "")))
        for h in sorted(chosen, key=lambda h: (h["seq"], h["start"])):
            sid, fam, acc = f"{stem}_{h['element']}", h["fam"], h["seq"]
            open(os.path.join(work, f"{sid}.fa"), "w").write(f">{sid}\n{genomes[acc][h['start'] - 1:h['end']]}\n")
            open(os.path.join(work, f"{fam}.cm"), "wb").write(docker("cmfetch", "/rfam/Rfam.cm", fam).stdout)
            open(os.path.join(work, f"{fam}.seed.sto"), "wb").write(docker("esl-afetch", "/rfam/Rfam.seed", fam).stdout)
            docker("cmalign", "--outformat", "Pfam", "--mapali", f"{fam}.seed.sto", "--mapstr", "-o", f"{sid}.sto", f"{fam}.cm", f"{sid}.fa")
            aligned, ss_cons = read_pfam_sto(os.path.join(work, f"{sid}.sto"))
            nested = "".join("." if c.isalpha() else c for c in ss_cons)
            seq, dot, lo, hi = project(aligned[sid], pair_columns(nested))
            if fam in keep_pk:
                _, pk, _, _ = project(aligned[sid], pk_columns(ss_cons))
                dot = "".join(p.translate(str.maketrans("()", "[]")) if p != "." else d for d, p in zip(dot, pk))
            seq, dot = seq[lo:hi + 1], open_short_loops(clip(dot[lo:hi + 1]))
            gstart, gend, bp = h["start"] + lo, h["start"] + hi, dot.count("(") + dot.count("[")
            if bp / len(seq) < MIN_BP_PER_NT:
                print(f"[{stem}] {h['element']} {fam} {acc}:{gstart}-{gend}: {bp} bp over {len(seq)} nt, dropped"); continue
            targets = [("ncbi", acc, gstart, gend)] + [("rnac", r, gstart - s + 1, gend - s + 1) for r, (c, s, e) in records.items()
                                                       if c == acc and s <= gstart and gend <= e]
            for source, ref, start, end in targets:
                oid = f"{sid}_{source}" if suffix else sid
                if not suffix and source != "ncbi":
                    continue
                db.append(f">{oid}\n{seq}\n{dot}\n")
                windows.append(f"{ref}\t{start}\t{end}\t{oid}\n")
                manifest.append(f"{oid}\t{stem}\t{fam}\t{h['name']}\t{acc}\t{gstart}\t{gend}\t{ref}\t{start}\t{end}\t{len(seq)}\t{dot.count('(')}\t{dot.count('[')}\t{h['evalue']:.2g}\n")
            print(f"[{stem}] {h['element']:16s} {fam} {acc}:{gstart}-{gend} {len(seq)} nt {dot.count('(')} bp + {dot.count('[')} pk "
                  f"E={h['evalue']:.2g} on {len(targets)} record(s)")

    if not db:
        print(f"{name}: no Rfam GA hit kept; no rf-eval set written"); return
    os.makedirs(out, exist_ok=True)
    setname = org.get("set", f"{name}_elements")
    dbp, winp = os.path.join(out, f"{setname}.db"), os.path.join(out, f"{setname}.windows.tsv")
    # an existing set keeps its other strains' entries; this run's stems are replaced
    old_db = open(dbp).read().split(">")[1:] if org.get("set") and os.path.exists(dbp) else []
    old_win = [l for l in open(winp)][1:] if org.get("set") and os.path.exists(winp) else []
    mine = lambda sid: any(sid.startswith(f"{s}_") for s in stems)
    db = [f">{b}" for b in old_db if not mine(b.split("\n")[0])] + db
    windows = [l for l in old_win if not mine(l.split("\t")[3].strip())] + windows
    open(dbp, "w").writelines(db)
    open(winp, "w").writelines(["# ref_seq_id\tstart\tend\tstructure_id\n"] + windows)
    man = os.path.join(orgdir, "MANIFEST-rfeval-added.tsv" if org.get("set") else "MANIFEST-rfeval.tsv")
    open(man, "w").writelines(["structure\tstem\trfam\trfam_id\tgenome\tgenome_start\tgenome_end\tref_seq_id\tstart\tend\t"
                               "length\tbase_pairs\tpk_pairs\tEvalue\n"] + manifest)


if __name__ == "__main__":
    main()
