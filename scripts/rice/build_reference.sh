#!/usr/bin/env bash
# Filter RNAcentral's rice GFF3 for one release into release<N>/rice/reference/.
# Usage, from the repository root: scripts/rice/build_reference.sh <release number>
# --no-prefilter: this species' main sources are not on the per-biotype trusted-database lists.
set -euo pipefail
release=$1
sp=oryza_sativa.IRGSP-1.0
out=release$release/rice
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
mkdir -p "$out/reference"
curl -fsSL -o "$tmp/$sp.gff3.gz" "https://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/$release.0/genome_coordinates/gff3/$sp.gff3.gz"
python3 scripts/filter_rnacentral_gff3.py --no-prefilter --pirna-max-loci 5 "$tmp/$sp.gff3.gz" "$tmp/filtered.gff3" 2> "$out/filter.log"
python3 scripts/gff3_to_gtf_exons.py --input "$tmp/filtered.gff3" --output "$out/reference/$sp.filtered.gtf" >> "$out/filter.log"
gzip -9 -n -f "$out/reference/$sp.filtered.gtf"
cat "$out/filter.log"
