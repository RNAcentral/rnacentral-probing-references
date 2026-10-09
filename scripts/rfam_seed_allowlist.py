#!/usr/bin/env python3
"""Build RFAM_SEED_ALLOWLIST entries for filter_rnacentral_gff3.py.

Lists the transcripts in a RNAcentral genome-coordinate GFF3 that the
authoritative-DB pre-filter would drop but whose sequence exactly matches a
sequence in an Rfam seed alignment. Output is grouped by biotype in the
format of the allowlist block in the filter script.

Usage: python rfam_seed_allowlist.py Rfam.seed.gz input.gff3[.gz] rnacentral_<species>.fasta[.gz]
"""
import gzip
import sys
from collections import defaultdict

from filter_rnacentral_gff3 import _open, _parse_attrs, has_authoritative_db


def seed_sequences(path: str) -> dict[str, str]:
    """Ungapped seed sequence (DNA alphabet) -> Rfam family accession."""
    seqs: dict[str, str] = {}
    family = ""
    with gzip.open(path, "rt", encoding="latin-1") as fh:
        for line in fh:
            if line.startswith("#=GF AC"):
                family = line.split()[2]
            elif line and line[0] not in "#/\n":
                _, seq = line.split(None, 1)
                seq = seq.strip().replace(".", "").replace("-", "").upper().replace("U", "T")
                seqs.setdefault(seq, family)
    return seqs


def main() -> None:
    seed_path, gff3_path, fasta_path = sys.argv[1:4]
    seeds = seed_sequences(seed_path)
    print(f"{len(seeds):,} distinct seed sequences", file=sys.stderr)

    candidates: dict[str, str] = {}
    with _open(gff3_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "transcript":
                continue
            attrs = _parse_attrs(parts[8])
            biotype = attrs.get("type", "")
            if biotype in ("circRNA", "piRNA"):
                continue
            if not has_authoritative_db(attrs.get("databases", ""), biotype):
                candidates[attrs.get("Name", "")] = biotype
    print(f"{len(candidates):,} URS IDs fail the pre-filter", file=sys.stderr)

    hits: dict[str, list[tuple[str, str]]] = defaultdict(list)
    name, chunks = "", []

    def flush() -> None:
        if name in candidates:
            seq = "".join(chunks).upper().replace("U", "T")
            if seq in seeds:
                hits[candidates[name]].append((name, seeds[seq]))

    with _open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                flush()
                name, chunks = line[1:].split(None, 1)[0], []
            else:
                chunks.append(line.strip())
        flush()

    total = sum(len(v) for v in hits.values())
    print(f"{total:,} match an Rfam seed sequence", file=sys.stderr)
    for biotype, entries in sorted(hits.items(), key=lambda kv: -len(kv[1])):
        families = sorted({f for _, f in entries})
        print(f"    # {biotype} ({len(entries)} entries; {', '.join(families)})")
        ids = sorted(n for n, _ in entries)
        for i in range(0, len(ids), 3):
            print("    " + " ".join(f'"{n}",' for n in ids[i:i + 3]))


if __name__ == "__main__":
    main()
