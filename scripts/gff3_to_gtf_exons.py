#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert a GFF3 file into a minimal exon-only GTF."
    )
    parser.add_argument("--input", required=True, help="Input GFF3 or GFF3.GZ path")
    parser.add_argument("--output", required=True, help="Output GTF path")
    return parser.parse_args()


def open_text(path: Path, mode: str):
    if path.suffix == ".gz":
        return gzip.open(path, mode, encoding="utf-8")
    return path.open(mode, encoding="utf-8")


def parse_attributes(field: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in field.strip().split(";"):
        if not item or "=" not in item:
            continue
        key, value = item.split("=", 1)
        attrs[key] = value
    return attrs


def quote_gtf(value: str) -> str:
    return value.replace("\\", "\\\\").replace('"', '\\"')


def exon_parent_ids(attrs: dict[str, str]) -> list[str]:
    parent_field = attrs.get("Parent", "").strip()
    if not parent_field:
        return []
    return [parent.strip() for parent in parent_field.split(",") if parent.strip()]


def _unique_sorted(csv: str) -> str:
    if not csv:
        return ""
    return ",".join(sorted({d.strip() for d in csv.split(",") if d.strip()}))


def _chrom_sort_key(chrom: str, start: int, end: int, feature_rank: int) -> tuple:
    """Sort key: numeric chromosomes first (1..22), then MT, X, Y alphabetically."""
    try:
        return (0, int(chrom), start, end, feature_rank)
    except ValueError:
        return (1, chrom, start, end, feature_rank)


def first_pass_collect_transcripts(input_path: Path) -> dict[str, dict[str, str]]:
    """Collect gene_id, biotype, databases, source_db, and coordinates keyed by raw transcript ID."""
    transcripts: dict[str, dict[str, str]] = {}
    with open_text(input_path, "rt") as handle:
        for raw_line in handle:
            if not raw_line or raw_line.startswith("#"):
                continue
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            if fields[2] != "transcript":
                continue
            attrs = parse_attributes(fields[8])
            transcript_id = attrs.get("ID")
            if not transcript_id:
                continue
            raw_gene_id = attrs.get("Parent", "").split(",", 1)[0].strip() or transcript_id
            transcripts[transcript_id] = {
                "gene_id":   raw_gene_id,
                "biotype":   attrs.get("type", ""),
                "databases": _unique_sorted(attrs.get("databases", "")),
                "source_db": _unique_sorted(attrs.get("providing_databases", "")),
                "chrom":     fields[0],
                "start":     int(fields[3]),
                "end":       int(fields[4]),
            }
    return transcripts


def _build_gtf_attrs(meta: dict, raw_tid: str) -> str:
    gene_id   = meta.get("gene_id", raw_tid)
    biotype   = meta.get("biotype", "")
    databases = meta.get("databases", "")
    source_db = meta.get("source_db", "")
    attrs = f'gene_id "{quote_gtf(gene_id)}"; transcript_id "{quote_gtf(raw_tid)}";'
    if biotype:
        attrs += f' transcript_biotype "{quote_gtf(biotype)}";'
    if databases:
        attrs += f' databases "{quote_gtf(databases)}";'
    if source_db:
        attrs += f' source_db "{quote_gtf(source_db)}";'
    return attrs


def convert_gff3_to_gtf(input_path: Path, output_path: Path) -> tuple[int, int]:
    tx_meta = first_pass_collect_transcripts(input_path)

    # The RNAcentral GFF3 places noncoding_exon lines BEFORE their parent
    # transcript line, so we buffer by transcript and flush in the right order.
    tx_gtf:   dict[str, str]        = {}             # raw_tid → transcript GTF line
    exon_gtf: dict[str, list[str]]  = {}             # raw_tid → exon GTF lines
    tx_seen:  list[str]             = []             # insertion-order transcript IDs

    with open_text(input_path, "rt") as reader:
        for raw_line in reader:
            if not raw_line or raw_line.startswith("#"):
                continue
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            feature_type = fields[2]

            if feature_type == "transcript":
                attrs = parse_attributes(fields[8])
                raw_tid = attrs.get("ID", "")
                if not raw_tid:
                    continue
                meta = tx_meta.get(raw_tid, {})
                out = fields[:]
                out[2] = "transcript"
                out[8] = _build_gtf_attrs(meta, raw_tid)
                tx_gtf[raw_tid] = "\t".join(out) + "\n"
                if raw_tid not in exon_gtf:
                    exon_gtf[raw_tid] = []
                    tx_seen.append(raw_tid)

            elif feature_type in {"exon", "noncoding_exon"}:
                attrs = parse_attributes(fields[8])
                for raw_tid in exon_parent_ids(attrs):
                    meta = tx_meta.get(raw_tid, {})
                    out = fields[:]
                    out[2] = "exon"
                    out[8] = _build_gtf_attrs(meta, raw_tid)
                    if raw_tid not in exon_gtf:
                        exon_gtf[raw_tid] = []
                        tx_seen.append(raw_tid)
                    exon_gtf[raw_tid].append("\t".join(out) + "\n")

    # Sort transcripts by genomic coordinate; sort exons within each transcript too.
    def _tx_key(raw_tid: str) -> tuple:
        m = tx_meta.get(raw_tid, {})
        return _chrom_sort_key(
            m.get("chrom", "ZZ"), m.get("start", 0), m.get("end", 0), 0
        )

    def _exon_key(line: str) -> tuple:
        f = line.split("\t")
        return _chrom_sort_key(f[0], int(f[3]), int(f[4]), 1)

    def _line_key(line: str) -> tuple:
        f = line.split("\t")
        return _chrom_sort_key(f[0], int(f[3]), int(f[4]), 0 if f[2] == "transcript" else 1)

    all_lines: list[str] = []
    for raw_tid in tx_seen:
        if raw_tid in tx_gtf:
            all_lines.append(tx_gtf[raw_tid])
        all_lines.extend(exon_gtf.get(raw_tid, []))

    written = 0
    with output_path.open("wt", encoding="utf-8") as writer:
        for line in sorted(all_lines, key=_line_key):
            writer.write(line)
            written += 1

    skipped = len(tx_meta) - len(tx_gtf)
    return written, skipped


def main() -> int:
    args = parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output)
    written, skipped = convert_gff3_to_gtf(input_path, output_path)
    print(f"Wrote {written} exon record(s); skipped {skipped} malformed/parentless line(s).")
    return 0 if written > 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
