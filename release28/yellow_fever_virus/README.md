# Yellow fever virus (rnastruct00075 17D, 00105 Asibi, 00106 Dakar) — references and Rfam structures for rf-eval

## References

One FASTA per strain (`build_yf_refs.py`, `MANIFEST-genome.tsv`):

| strain | datasets | genome | note |
|---|---|---|---|
| 17D | 00075 | PX660172.1, 10 862 nt | the 17D-204 full-length infectious clone; the metadata says not to use NC_002031 / X03700 |
| Asibi | 00105 | AY640589.1, 10 862 nt | complete genome |
| Dakar | 00106 | MN106242.1, 10 784 nt | isolate Dakar HD 1279; its 3′UTR is about 77 nt shorter (431 nt against 508) |

RNAcentral records of ≥ 60 nt at ≥ 99 % identity to the genome would become
transcripts of their own (the dengue sfRNA rule); none of the three strains has
one, so each FASTA holds the genome and the `<accession>_5UTR`/`_3UTR` transcripts
either side of the polyprotein CDS (118 and 508 nt for 17D and Asibi, 117 and
431 nt for Dakar). Where records overlap only the best by identity, then length,
then ID is kept.

## rf-eval reference

`scripts/rfeval_viral_elements.py yellow_fever_virus 28 <Rfam dir>`: `cmscan --cut_ga` over
the full Rfam 15.1 with clan competition, each hit `cmalign`ed back to its model with the
seed structure mapped in, the `SS_cons` projected and trimmed to the model-covered
span, pseudoknots dropped, and only structures of at least 0.10 bp/nt kept.
Windows are on the genome record (`_ncbi`).

| element | Rfam | 17D | Asibi | nt | bp |
|---|---|---|---|---|---|
| 3SL | RF00185 Flavi_CRE | 10762–10860 | 10762–10860 | 99 | 24 |

2 entries. Nothing reaches GA on Dakar, so there is nothing to score on MN106242.1.

Rebuild the reference: `cd release<N>/yellow_fever_virus && python3 ../../scripts/yellow_fever_virus/build_yf_refs.py` (Biopython, network) writes the genome FASTA/GTF to `work/`; then, from the root,
`python3 scripts/make_transcriptome_fa.py release<N>/yellow_fever_virus` cuts `reference/<stem>.transcriptome.fa` from them.
rf-eval rebuild: `scripts/rfeval_viral_elements.py yellow_fever_virus <N> <Rfam dir>` (Docker, R2DT image for Infernal).

## Files

- `reference/YFV_<strain>.transcriptome.fa` — the genome plus its 5′/3′ UTR; `--fasta` with
  `--transcriptome`. `MANIFEST-genome.tsv` lists them.
- `rfeval/yellow_fever_elements.db` — 17D + Asibi. `--rfeval_reference`.
- `rfeval/yellow_fever_elements.windows.tsv` — `ref_seq_id start end structure_id`, 1-based
  inclusive. `--rfeval_windows`.
- `MANIFEST-rfeval.tsv` — element, family, coordinates, pair counts, E-value.

## Running

```
nextflow run . -profile <…> -resume \
  --transcriptome --fasta release28/yellow_fever_virus/reference/YFV_17D.transcriptome.fa \
  --rfeval_reference release28/yellow_fever_virus/rfeval/yellow_fever_elements.db \
  --rfeval_windows   release28/yellow_fever_virus/rfeval/yellow_fever_elements.windows.tsv
```

Swap `17D` for `Asibi` for rnastruct00105. For rnastruct00106 there is nothing to
score on MN106242.1; omit the two `rfeval_*` params.
