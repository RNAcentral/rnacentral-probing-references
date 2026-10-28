# Porcine epidemic diarrhea virus AJ1102 (rnastruct00074) — reference and Rfam structures for rf-eval

## Reference

`PEDV_AJ1102.transcriptome.fa`: the genome MK584552.1 (28 044 nt) and
`MK584552.1_5UTR` (292 nt) / `_3UTR` (343 nt), cut at the first CDS start and last CDS
end (`build_pedv_refs.py`, `MANIFEST-genome.tsv`). RNAcentral records of ≥ 60 nt at
≥ 99 % identity would become transcripts of their own; none sits in this genome.

## rf-eval reference

`scripts/rfeval_viral_elements.py porcine_epidemic_diarrhea_virus 28 <Rfam dir>`: `cmscan --cut_ga` over
the full Rfam 15.1 with clan competition, each hit `cmalign`ed back to its model with the
seed structure mapped in, the `SS_cons` projected and trimmed to the model-covered
span, pseudoknots dropped, and only structures of at least 0.10 bp/nt kept.
Windows are on the genome record (`_ncbi`).

| element | Rfam | MK584552.1 | nt | bp | E |
|---|---|---|---|---|---|
| aCoV-5UTR | RF03116 | 1–320 | 320 | 88 | 2.3e-67 |
| Corona_FSE | RF00507 | 12617–12691 | 75 | 17 | 4.8e-16 |
| aCoV-3UTR | RF03121 | 27735–28035 | 301 | 47 | 6.2e-52 |

3 entries, `_ncbi` only (`PEDV_AJ1102_<element>_ncbi`).

Rebuild the reference: `cd release<N>/porcine_epidemic_diarrhea_virus && python3 ../../scripts/porcine_epidemic_diarrhea_virus/build_pedv_refs.py` (Biopython, network) writes the genome FASTA/GTF to `work/`; then, from the root,
`python3 scripts/make_transcriptome_fa.py release<N>/porcine_epidemic_diarrhea_virus` cuts `reference/<stem>.transcriptome.fa` from them.
rf-eval rebuild: `scripts/rfeval_viral_elements.py porcine_epidemic_diarrhea_virus <N> <Rfam dir>` (Docker, R2DT image for Infernal).

## Files

- `reference/PEDV_AJ1102.transcriptome.fa` — the genome plus its 5′/3′ UTR; `--fasta` with
  `--transcriptome`. `MANIFEST-genome.tsv` lists them.
- `rfeval/pedv_elements.db` — `--rfeval_reference`.
- `rfeval/pedv_elements.windows.tsv` — `ref_seq_id start end structure_id`, 1-based
  inclusive. `--rfeval_windows`.
- `MANIFEST-rfeval.tsv` — element, family, coordinates, pair counts, E-value.

## Running

```
nextflow run . -profile <…> -resume \
  --transcriptome --fasta release28/porcine_epidemic_diarrhea_virus/reference/PEDV_AJ1102.transcriptome.fa \
  --rfeval_reference release28/porcine_epidemic_diarrhea_virus/rfeval/pedv_elements.db \
  --rfeval_windows   release28/porcine_epidemic_diarrhea_virus/rfeval/pedv_elements.windows.tsv
```
