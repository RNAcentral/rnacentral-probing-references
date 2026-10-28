# West Nile virus NY99 (rnastruct00136) — reference and Rfam structures for rf-eval

## Reference

`WNV_NY99.transcriptome.fa`: the NY99 infectious-clone genome AF404756.1 (11 029 nt)
and `AF404756.1_5UTR` (96 nt) / `_3UTR` (631 nt), either side of the polyprotein CDS
(`build_wnv_refs.py`, `MANIFEST-genome.tsv`). RNAcentral records of ≥ 60 nt at ≥ 99 %
identity would become transcripts of their own; none sits in this genome. Where
records overlap only the best by identity, then length, then ID is kept.

## rf-eval reference

`scripts/rfeval_viral_elements.py west_nile_virus 28 <Rfam dir>`: `cmscan --cut_ga` over
the full Rfam 15.1 with clan competition, each hit `cmalign`ed back to its model with the
seed structure mapped in, the `SS_cons` projected and trimmed to the model-covered
span, pseudoknots dropped, and only structures of at least 0.10 bp/nt kept.
Windows are on the genome record (`_ncbi`).

| element | Rfam | AF404756.1 | nt | bp | E |
|---|---|---|---|---|---|
| Flavivirus-5UTR | RF03546 | 1–151 | 151 | 46 | 1.4e-30 |
| flavi_FSE | RF01768 | 3558–3619 | 62 | 11 | 9.2e-20 |
| Flavivirus_SLIV1 | RF01415 | 10501–10581 | 81 | 18 | 1.3e-07 |
| Flavivirus_SLIV2 | RF01415 | 10661–10735 | 75 | 19 | 3.1e-15 |
| Flavivirus_DB1 | RF00525 | 10760–10832 | 73 | 19 | 3.2e-15 |
| Flavivirus_DB2 | RF00525 | 10837–10906 | 70 | 18 | 1.7e-15 |
| Flavi_CRE | RF00185 | 10936–11029 | 94 | 24 | 4e-12 |

7 entries, `_ncbi` only (`WNV_NY99_<element>_ncbi`).

Rebuild the reference: `cd release<N>/west_nile_virus && python3 ../../scripts/west_nile_virus/build_wnv_refs.py` (Biopython, network) writes the genome FASTA/GTF to `work/`; then, from the root,
`python3 scripts/make_transcriptome_fa.py release<N>/west_nile_virus` cuts `reference/<stem>.transcriptome.fa` from them.
rf-eval rebuild: `scripts/rfeval_viral_elements.py west_nile_virus <N> <Rfam dir>` (Docker, R2DT image for Infernal).

## Files

- `reference/WNV_NY99.transcriptome.fa` — the genome plus its 5′/3′ UTR; `--fasta` with
  `--transcriptome`. `MANIFEST-genome.tsv` lists them.
- `rfeval/west_nile_elements.db` — `--rfeval_reference`.
- `rfeval/west_nile_elements.windows.tsv` — `ref_seq_id start end structure_id`, 1-based
  inclusive. `--rfeval_windows`.
- `MANIFEST-rfeval.tsv` — element, family, coordinates, pair counts, E-value.

## Running

```
nextflow run . -profile <…> -resume \
  --transcriptome --fasta release28/west_nile_virus/reference/WNV_NY99.transcriptome.fa \
  --rfeval_reference release28/west_nile_virus/rfeval/west_nile_elements.db \
  --rfeval_windows   release28/west_nile_virus/rfeval/west_nile_elements.windows.tsv
```
