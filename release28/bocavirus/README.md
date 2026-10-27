# Human bocavirus 1 Salvador1 (rnastruct00118) — reference

## Reference

`HBOV1_Salvador1.transcriptome.fa`: the infectious-clone genome JQ923422.1 (5 543 nt)
and `JQ923422.1_5UTR` (348 nt) / `_3UTR` (376 nt), cut at the first CDS start and last
CDS end (`build_hbov_refs.py`, `MANIFEST-genome.tsv`). RNAcentral records of ≥ 60 nt at
≥ 99 % identity would become transcripts of their own; none sits in this genome.

## rf-eval reference

None. No Rfam family reaches GA on the genome, so there is nothing to score.

Rebuild the reference: `cd release<N>/bocavirus && python3 ../../scripts/bocavirus/build_hbov_refs.py` (Biopython, network) writes the genome FASTA/GTF to `work/`; then, from the root,
`python3 scripts/make_transcriptome_fa.py release<N>/bocavirus` cuts `reference/<stem>.transcriptome.fa` from them.

## Files

- `reference/HBOV1_Salvador1.transcriptome.fa` — the genome plus its 5′/3′ UTR; `--fasta` with
  `--transcriptome`. `MANIFEST-genome.tsv` lists them.

## Running

```
nextflow run . -profile <…> -resume \
  --transcriptome --fasta release28/bocavirus/reference/HBOV1_Salvador1.transcriptome.fa
```
