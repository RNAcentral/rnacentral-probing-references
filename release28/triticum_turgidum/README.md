# Durum wheat (T. turgidum Svevo.v1) — release 28

## Reference

- **GTF:** `reference/triticum_turgidum.Svevo.v1.filtered.gtf.gz`, the RNAcentral release 28
  ncRNAs on Svevo.v1 after filtering (below).
- **FASTA:** not stored. The pipeline downloads the Ensembl Plants Svevo.v1 genome itself
  (`organism: triticum_turgidum`). Chromosome names (`1A` … `7B`) match that FASTA.

Transcript IDs are `URS…_4567.<n>`, one ID per genomic copy. The metadata's taxid
4567 is run on the Svevo assembly.

## How the GTF was filtered

Input is `triticum_turgidum.Svevo.v1.gff3.gz` from
`https://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/28.0/genome_coordinates/gff3/`.

```
scripts/triticum_turgidum/build_reference.sh 28   # from the repository root
```

It runs `scripts/filter_rnacentral_gff3.py --no-prefilter --pirna-max-loci 5`, then
`scripts/gff3_to_gtf_exons.py`, and writes `filter.log` beside this README.

`--no-prefilter` skips the trusted-database check, because wheat's main sources
(ENA, SILVA, Rfam) are not on the per-RNA-type lists, which were written for human;
a trial run with the check on removed 1,784 of the 2,193 transcripts (81 %).
Main chromosomes are still enforced: the unplaced `Un` sequence is excluded. Wheat
has no circRNAs or piRNAs, so step 0 and the piRNA cap do nothing.

| step | what it does | durum_wheat |
|---|---|---|
| 0 | drop circRNAs | none; 2,193 transcripts |
| 1 | same sequence listed several times at one locus: keep the best-ranked copy | 1,645 kept |
| 2 | overlapping records of one RNA type on one strand: keep one per RNAcentral predicted gene (or one per group) | **1,032 kept** |

Biotypes kept: rRNA 577, tRNA 447, misc_RNA 7, snoRNA 1.

## rf-eval references

- `rfeval/rrna_common.*`: 18S, 26S and 5S (below). No 5.8S reached GA inside the rDNA unit.

## Common rRNA set (`rfeval/rrna_common.*`)

Built by `scripts/rfeval_rrna_common.py triticum_turgidum 28`: `cmsearch --cut_ga` on the
GTF's rRNA records (Ensembl REST sequences), one locus per molecule within 10 kb of
the SSU, anchored on the SSU whose unit holds the most linked molecules; 5S is
searched genome-wide. Structures are the Rfam `SS_cons` projected with `cmalign`,
without pseudoknots, and only structures of at least 0.10 bp/nt are kept.
`MANIFEST-rfeval.tsv` has the bit scores and loci.

| structure | transcript (filtered GTF) | window | locus | bp |
|---|---|---|---|---|
| `URS00007DC593_4567.200__18S_rRNA` | `URS00007DC593_4567.200` | 8–1485 | 1B:29877319-29878810 | 388 |
| `URS00007BFDE6_4567.206__26S_rRNA` | `URS00007BFDE6_4567.206` | 144–2879 | 1B:29882213-29885098 | 643 |
| `URS000035C26F_4567.443__5S_rRNA` | `URS000035C26F_4567.443` | 1–119 | 1B:37493377-37493495 | 34 |

Svevo.v1 does not carry complete rDNA units: the best 18S and 26S records are
partial (1,478 of about 1,810 nt scored; 632 bits, against 1,785 for Arabidopsis's
full-length 18S; 2,736 nt of 26S), and no 5.8S sits between them at GA. Treat the
set as a partial yardstick.

The names only match runs on
`release28/triticum_turgidum/reference/triticum_turgidum.Svevo.v1.filtered.gtf.gz`.

```
--gtf              release28/triticum_turgidum/reference/triticum_turgidum.Svevo.v1.filtered.gtf.gz \
--rfeval_reference release28/triticum_turgidum/rfeval/rrna_common.db \
--rfeval_windows   release28/triticum_turgidum/rfeval/rrna_common.windows.tsv
```
