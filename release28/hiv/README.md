# HIV-1 NL4-3 (rnastruct00119) — reference and Rfam structures for rf-eval

## Reference

`HIV1_NL43.transcriptome.fa`: the NL4-3 provirus, AF324493.2 bases 1–9709 (5′ LTR to
3′ LTR), plus two RNAcentral records and `AF324493.2_5UTR` (1–789) / `_3UTR`
(9408–9709), either side of gag start to nef end (`build_hiv_refs.py`,
`MANIFEST-genome.tsv`). AF324493.2 is the whole pNL4-3 plasmid (14 825 nt): after the
3′ LTR it runs into cellular flanking DNA and pUC vector, which is not viral RNA and is
cut off.

RNAcentral records of ≥ 60 nt at ≥ 99 % identity become transcripts of their own;
where records overlap only the best by identity, then length, then ID is kept:

| record | nt | identity | AF324493.2 | description |
|---|---|---|---|---|
| `URS00006F9B91_11676` | 61 | 100 % | 454–514 | hiv1-mir-TAR precursor |
| `URS0000D56C30_11676` | 101 | 100 % | 577–677 | viral RNA genome fragment (PDB 6WB0, chain C) |

The probed virus is VSV-G-pseudotyped and delta Env; the reference is the full NL4-3
provirus, so env reads will be sparse.

## rf-eval reference

`scripts/rfeval_viral_elements.py hiv 28 <Rfam dir>`: `cmscan --cut_ga` over the full
Rfam 15.1 with clan competition, each hit `cmalign`ed back to its model with the seed
structure mapped in, the `SS_cons` projected and trimmed to the model-covered span,
pseudoknots dropped, and only structures of at least 0.10 bp/nt kept. Hits are kept only
inside the RNA genome, R to R (454–9626): the provirus's 5′ U3 and 3′ U5 are not in the
viral RNA, so the second PBS hit (9654–9709) is left out. HIV_GSL3 (854–937) reaches GA
but has 8 bp over 84 nt and is dropped.

| element | Rfam | AF324493.2 | nt | bp | E | windows |
|---|---|---|---|---|---|---|
| mir-TAR1 | RF00250 | 454–514 | 61 | 24 | 6.4e-21 | `_ncbi`, `_rnac` |
| HIV_PBS | RF00375 | 579–677 | 99 | 32 | 2.2e-27 | `_ncbi`, `_rnac` |
| HIV-1_DIS | RF00175 | 694–733 | 40 | 11 | 1.6e-11 | `_ncbi` |
| HIV-1_SD | RF01380 | 736–754 | 19 | 7 | 3.4e-05 | `_ncbi` |
| HIV-1_SL3 | RF01381 | 763–785 | 23 | 5 | 5.8e-09 | `_ncbi` |
| HIV-1_SL4 | RF01382 | 791–810 | 20 | 5 | 0.00062 | `_ncbi` |
| HIV_FE | RF00480 | 2085–2136 | 52 | 18 | 1.8e-12 | `_ncbi` |
| HIV_POL-1_SL | RF01418 | 2466–2578 | 113 | 34 | 4.9e-30 | `_ncbi` |
| RRE | RF00036 | 7709–8045 | 337 | 116 | 2.7e-125 | `_ncbi` |
| mir-TAR2 | RF00250 | 9529–9589 | 61 | 24 | 6.4e-21 | `_ncbi` |

12 entries: every element on the genome record (`HIV1_NL43_<element>_ncbi`), and the
two inside an RNAcentral record also on that record (`_rnac`). TAR appears twice
because R is repeated in both LTRs.

Rebuild the reference: `cd release<N>/hiv && python3 ../../scripts/hiv/build_hiv_refs.py` (Biopython, network) writes the genome FASTA/GTF to `work/`; then, from the root,
`python3 scripts/make_transcriptome_fa.py release<N>/hiv` cuts `reference/<stem>.transcriptome.fa` from them.
rf-eval rebuild: `scripts/rfeval_viral_elements.py hiv <N> <Rfam dir>` (Docker, R2DT image for Infernal).

## Files

- `reference/HIV1_NL43.transcriptome.fa` — the provirus, the two records and the 5′/3′
  UTR; `--fasta` with `--transcriptome`. `MANIFEST-genome.tsv` lists them.
- `rfeval/hiv_elements.db` — `--rfeval_reference`.
- `rfeval/hiv_elements.windows.tsv` — `ref_seq_id start end structure_id`, 1-based
  inclusive. `--rfeval_windows`.
- `MANIFEST-rfeval.tsv` — element, family, coordinates, pair counts, E-value.

## Running

```
nextflow run . -profile <…> -resume \
  --transcriptome --fasta release28/hiv/reference/HIV1_NL43.transcriptome.fa \
  --rfeval_reference release28/hiv/rfeval/hiv_elements.db \
  --rfeval_windows   release28/hiv/rfeval/hiv_elements.windows.tsv
```
