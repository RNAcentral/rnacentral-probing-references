#!/usr/bin/env python3
"""Project an Rfam SS_cons onto each cmalign-aligned sequence as a dot-bracket.

cmalign places every query residue in either a model (match) column or an insert
column, so a consensus base pair transfers to the query whenever both of its
columns hold a residue in that query. Pairs with a gapped partner are dropped --
the query has no partner to pair with, so it is called unpaired.

Each structure is trimmed to the model-covered span: everything from the first
to the last residue sitting in a match column. Flanks the model never covered
would otherwise be scored as confidently-unpaired ground truth when the model
simply says nothing about them (yeast U2 is 1175 nt against a ~190 nt model).
The span is written to a windows manifest so the XML can be sliced to match.

WUSS pseudoknot annotation (A-Z/a-z) is dropped: dot-bracket cannot express it
and the R2DT arm is nested-only too, so the two arms stay comparable.
"""
import argparse, sys

OPEN_TO_CLOSE = {"<": ">", "(": ")", "[": "]", "{": "}"}
CLOSE_TO_OPEN = {v: k for k, v in OPEN_TO_CLOSE.items()}
GAP = ".-_~"


def read_pfam_sto(path):
    """Parse a Pfam-format (one line per sequence) Stockholm into (seqs, ss_cons)."""
    seqs, ss_cons = {}, None
    for line in open(path, encoding="utf-8"):
        line = line.rstrip("\n")
        if line.startswith("#=GC SS_cons"):
            ss_cons = line.split(None, 2)[2]
        elif line.startswith("#") or line.startswith("//") or not line.strip():
            continue
        else:
            name, aln = line.split(None, 1)
            if name in seqs:
                sys.exit(f"{path}: {name} twice; expected Pfam (single-block) format")
            seqs[name] = aln.strip()
    if ss_cons is None:
        sys.exit(f"{path}: no #=GC SS_cons line")
    return seqs, ss_cons


def pair_columns(ss_cons):
    """Map alignment column -> partner column, for the nested WUSS brackets only."""
    stacks, pairs = {k: [] for k in OPEN_TO_CLOSE}, {}
    for col, ch in enumerate(ss_cons):
        if ch in OPEN_TO_CLOSE:
            stacks[ch].append(col)
        elif ch in CLOSE_TO_OPEN:
            opener = CLOSE_TO_OPEN[ch]
            if not stacks[opener]:
                sys.exit(f"unbalanced SS_cons at column {col}")
            i = stacks[opener].pop()
            pairs[i], pairs[col] = col, i
    unclosed = [k for k, v in stacks.items() if v]
    if unclosed:
        sys.exit(f"unclosed SS_cons brackets: {unclosed}")
    return pairs


def project(aln, col_pairs):
    """Return (sequence, dot-bracket, first_match_pos, last_match_pos), 0-based positions.

    In Pfam-format cmalign output a residue in a match column is uppercase and a
    residue in an insert column is lowercase, which is how the model-covered span
    is recovered.
    """
    seq, is_match = [], []
    col_to_pos = {}
    for col, ch in enumerate(aln):
        if ch in GAP:
            continue
        col_to_pos[col] = len(seq)
        is_match.append(ch.isupper())
        seq.append(ch.upper().replace("T", "U"))
    db = ["."] * len(seq)
    for col, partner in col_pairs.items():
        if col >= partner:
            continue  # each pair handled once, from its opener
        if col in col_to_pos and partner in col_to_pos:
            db[col_to_pos[col]] = "("
            db[col_to_pos[partner]] = ")"
    matched = [i for i, m in enumerate(is_match) if m]
    if not matched:
        return "".join(seq), "".join(db), None, None
    return "".join(seq), "".join(db), matched[0], matched[-1]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sto", required=True, nargs="+", help="cmalign --outformat Pfam alignments")
    ap.add_argument("--output", required=True, help="output .db")
    ap.add_argument("--windows", required=True, help="output windows manifest")
    ap.add_argument("--manifest", required=True, help="output TSV of id/family/span/bp")
    ap.add_argument("--min-bp-density", type=float, default=0.10,
                    help="drop structures below this bp-per-nt (default: %(default)s)")
    args = ap.parse_args()

    entries, windows, rows, dropped = {}, [], [], []
    for sto in args.sto:
        family = sto.rsplit("/", 1)[-1].split(".")[0]
        seqs, ss_cons = read_pfam_sto(sto)
        col_pairs = pair_columns(ss_cons)
        for name, aln in sorted(seqs.items()):
            seq, db, lo, hi = project(aln, col_pairs)
            if lo is None:
                dropped.append(f"{name}: no match columns in {family}")
                continue
            sub_seq, sub_db = seq[lo:hi + 1], db[lo:hi + 1]
            # a pair whose partner fell outside the trimmed span is no longer a pair
            stack, out = [], list(sub_db)
            for i, ch in enumerate(sub_db):
                if ch == "(":
                    stack.append(i)
                elif ch == ")":
                    if stack:
                        stack.pop()
                    else:
                        out[i] = "."
            for i in stack:
                out[i] = "."
            sub_db = "".join(out)
            npairs = sub_db.count("(")
            density = npairs / float(len(sub_seq)) if sub_seq else 0.0
            kind = "full" if (lo, hi) == (0, len(seq) - 1) else "partial"
            if density < args.min_bp_density:
                dropped.append(f"{name}: {npairs} bp over {len(sub_seq)} nt = "
                               f"{density:.3f} bp/nt, below --min-bp-density ({family})")
                rows.append((name, family, "dropped_low_bp", lo + 1, hi + 1,
                             len(sub_seq), len(seq), npairs))
                continue
            if name in entries:
                dropped.append(f"{name}: aligned to more than one family; keeping the first")
                continue
            entries[name] = (sub_seq, sub_db)
            if kind == "partial":
                windows.append((name, lo + 1, hi + 1, name))
            rows.append((name, family, kind, lo + 1, hi + 1, len(sub_seq), len(seq), npairs))

    with open(args.output, "w", encoding="utf-8") as fh:
        for name in sorted(entries):
            seq, db = entries[name]
            fh.write(f">{name}\n{seq}\n{db}\n")
    with open(args.windows, "w", encoding="utf-8") as fh:
        fh.write("# ref_seq_id\tstart\tend\tstructure_id -- model-covered span, partial hits only\n")
        for row in sorted(windows):
            fh.write("\t".join(str(v) for v in row) + "\n")
    with open(args.manifest, "w", encoding="utf-8") as fh:
        fh.write("structure\trfam\tkind\tstart\tend\tstructure_len\ttranscript_len\tbase_pairs\n")
        for row in sorted(rows):
            fh.write("\t".join(str(v) for v in row) + "\n")

    for msg in dropped:
        print(f"[project] drop {msg}", file=sys.stderr)
    nfull = sum(1 for r in rows if r[2] == "full")
    print(f"[project] {len(entries)} structures ({nfull} full, {len(windows)} partial)")


if __name__ == "__main__":
    main()
