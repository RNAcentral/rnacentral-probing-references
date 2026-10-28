# Murine norovirus MNV-1 CW3 (rnastruct00097/00098) — reference and Rfam structure for rf-eval

## Reference

`MNV_CW3.transcriptome.fa`: the genome EF014462.1 (7 382 nt) and `EF014462.1_3UTR`
(75 nt), after the last CDS end (ORF1–3; `build_mnv_refs.py`, `MANIFEST-genome.tsv`).
The 5′UTR is 5 nt, under the 20 nt minimum, so there is no `_5UTR` transcript.
RNAcentral records of ≥ 60 nt at ≥ 99 % identity would become transcripts of their
own; none sits in this genome.

## rf-eval reference

`scripts/rfeval_viral_elements.py murine_norovirus 28 <Rfam dir>`: `cmscan --cut_ga` over
the full Rfam 15.1 with clan competition, each hit `cmalign`ed back to its model with the
seed structure mapped in, the `SS_cons` projected and trimmed to the model-covered
span, pseudoknots dropped, and only structures of at least 0.10 bp/nt kept.

| element | Rfam | EF014462.1 | nt | bp | E |
|---|---|---|---|---|---|
| MNV_3UTR | RF02532 | 7161–7382 | 222 | 74 | 1.1e-63 |

1 entry, `_ncbi` only (`MNV_CW3_MNV_3UTR_ncbi`). The element starts before the
`_3UTR` transcript (75 nt), so it is scored on the genome record.

Rebuild the reference: `cd release<N>/murine_norovirus && python3 ../../scripts/murine_norovirus/build_mnv_refs.py` (Biopython, network) writes the genome FASTA/GTF to `work/`; then, from the root,
`python3 scripts/make_transcriptome_fa.py release<N>/murine_norovirus` cuts `reference/<stem>.transcriptome.fa` from them.
rf-eval rebuild: `scripts/rfeval_viral_elements.py murine_norovirus <N> <Rfam dir>` (Docker, R2DT image for Infernal).

## Files

- `reference/MNV_CW3.transcriptome.fa` — the genome plus its 3′UTR; `--fasta` with
  `--transcriptome`. `MANIFEST-genome.tsv` lists them.
- `rfeval/norovirus_elements.db` — `--rfeval_reference`.
- `rfeval/norovirus_elements.windows.tsv` — `ref_seq_id start end structure_id`, 1-based
  inclusive. `--rfeval_windows`.
- `MANIFEST-rfeval.tsv` — element, family, coordinates, pair count, E-value.

## Running

```
nextflow run . -profile <…> -resume \
  --transcriptome --fasta release28/murine_norovirus/reference/MNV_CW3.transcriptome.fa \
  --rfeval_reference release28/murine_norovirus/rfeval/norovirus_elements.db \
  --rfeval_windows   release28/murine_norovirus/rfeval/norovirus_elements.windows.tsv
```
