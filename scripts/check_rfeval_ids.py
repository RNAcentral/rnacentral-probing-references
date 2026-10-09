#!/usr/bin/env python3
"""Check an rf-eval windows file against a new release's GTF and re-point renumbered IDs (see README, step 3)."""
import argparse, gzip, os, re, sys

TID = re.compile(r'transcript_id "([^"]+)"')


def urs(tid):
    """URS…_<taxid>.<n> -> URS…_<taxid>; IDs without a copy number are returned unchanged."""
    return tid.rsplit(".", 1)[0] if "." in tid else tid


def read_windows(path):
    rows = []
    for line in open(path, encoding="utf-8"):
        if line.startswith("#") or not line.strip():
            continue
        rows.append(line.rstrip("\n").split("\t"))
    return rows


def gtf_exons(path, keep):
    """{transcript_id: sorted [(chrom, strand, start, end)]} for transcripts where keep(tid) is true."""
    exons = {}
    with (gzip.open(path, "rt", encoding="utf-8") if path.endswith(".gz") else open(path, encoding="utf-8")) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            m = TID.search(f[8])
            if m and keep(m.group(1)):
                exons.setdefault(m.group(1), []).append((f[0], f[6], int(f[3]), int(f[4])))
    return {t: sorted(e) for t, e in exons.items()}


def overlaps(a, b):
    return any(c1 == c2 and s1 == s2 and x1 <= y2 and x2 <= y1
               for c1, s1, x1, y1 in a for c2, s2, x2, y2 in b)


def locus(ex):
    return f"{ex[0][0]}:{ex[0][2]}-{ex[-1][3]}({ex[0][1]})"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("windows", help="rf-eval windows manifest, e.g. release27/human/rfeval/rrna_common.windows.tsv")
    ap.add_argument("old_gtf", help="GTF the windows were built on")
    ap.add_argument("new_gtf", help="GTF of the new release")
    ap.add_argument("--outdir", help="write re-pointed <name>.windows.tsv and <name>.db here")
    args = ap.parse_args()

    rows = read_windows(args.windows)
    wanted = sorted({r[0] for r in rows})
    bases = {urs(t) for t in wanted}
    old = gtf_exons(args.old_gtf, lambda t: t in wanted)
    new = gtf_exons(args.new_gtf, lambda t: urs(t) in bases)

    repoint, unresolved = {}, 0
    for tid in wanted:
        if tid not in old:
            sys.exit(f"{tid}: not in the old GTF; is it the GTF the windows were built on?")
        if tid in new:
            status = "ok" if new[tid] == old[tid] else f"moved {locus(old[tid])} -> {locus(new[tid])}"
            print(f"{tid}\t{status}")
            continue
        same_urs = sorted(t for t in new if urs(t) == urs(tid))
        here = [t for t in same_urs if new[t] == old[tid]] or [t for t in same_urs if overlaps(new[t], old[tid])]
        if here:
            repoint[tid] = here[0]
            extra = f" (also {', '.join(here[1:])})" if len(here) > 1 else ""
            print(f"{tid}\trepoint -> {here[0]} at {locus(new[here[0]])}{extra}")
        elif same_urs:
            unresolved += 1
            print(f"{tid}\tcandidates at other loci: " + ", ".join(f"{t} {locus(new[t])}" for t in same_urs))
        else:
            unresolved += 1
            print(f"{tid}\tmissing: no {urs(tid)} in the new GTF; rebuild this structure")

    if args.outdir:
        stem = os.path.basename(args.windows).removesuffix(".windows.tsv")
        db = os.path.join(os.path.dirname(args.windows), stem + ".db")

        os.makedirs(args.outdir, exist_ok=True)
        for src in (args.windows, db):
            with open(src, encoding="utf-8") as fh:
                text = fh.read()
            for a, b in repoint.items():  # not followed by a digit, so .12 never matches inside .123
                text = re.sub(re.escape(a) + r"(?!\d)", b, text)
            with open(os.path.join(args.outdir, os.path.basename(src)), "w", encoding="utf-8") as out:
                out.write(text)
        print(f"wrote {stem}.windows.tsv and {stem}.db to {args.outdir}")
    sys.exit(1 if unresolved else 0)


if __name__ == "__main__":
    main()
