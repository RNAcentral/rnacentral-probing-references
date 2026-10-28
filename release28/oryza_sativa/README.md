# Rice (O. sativa IRGSP-1.0) — release 28

## Reference

- **GTF:** `reference/oryza_sativa.IRGSP-1.0.filtered.gtf.gz`, the RNAcentral release 28
  ncRNAs on IRGSP-1.0 after filtering (below).
- **FASTA:** not stored. The pipeline downloads the Ensembl Plants IRGSP-1.0 genome
  itself (`organism: oryza_sativa`).

Transcript IDs are `URS…_39947.<n>`, one ID per genomic copy.

## How the GTF was filtered

Input is `oryza_sativa.IRGSP-1.0.gff3.gz` from
`https://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/28.0/genome_coordinates/gff3/`.

```
scripts/oryza_sativa/build_reference.sh 28   # from the repository root
```

It runs `scripts/filter_rnacentral_gff3.py --no-prefilter --pirna-max-loci 5`, then
`scripts/gff3_to_gtf_exons.py`, and writes `filter.log` beside this README.

`--no-prefilter` skips the trusted-database check, because rice's main sources
(Ensembl Plants, Rfam) are not on the per-RNA-type lists, which were written for
human. Main chromosomes are still enforced (numbered plus Mt/Pt). Rice has no
circRNAs or piRNAs, so step 0 and the piRNA cap do nothing.

| step | what it does | rice |
|---|---|---|
| 0 | drop circRNAs | none; 19,108 transcripts |
| 1 | same sequence listed several times at one locus: keep the best-ranked copy | 17,753 kept |
| 2 | overlapping records of one RNA type on one strand: keep one per RNAcentral predicted gene (or one per group) | **15,093 kept** |

Biotypes kept: pre_miRNA 8,749, ncRNA 2,337, miRNA 1,374, tRNA 892, snoRNA 773,
rRNA 594, snRNA 220, and smaller groups. The ncRNA group is labelled sRNA in the
release 27 GTF.

## rf-eval references

- `rfeval/rrna_common.*`: 18S, 5.8S, 25S and 5S (below). New in release 28; release 27
  had no rice rf-eval set.

## Common rRNA set (`rfeval/rrna_common.*`)

Built by `scripts/rfeval_rrna_common.py oryza_sativa 28`: `cmsearch --cut_ga` on the
GTF's rRNA records (Ensembl REST sequences), one locus per molecule within 10 kb of
the SSU, anchored on the SSU whose unit holds the most linked molecules; 5S is
searched genome-wide. Structures are the Rfam `SS_cons` projected with `cmalign`,
without pseudoknots, and only structures of at least 0.10 bp/nt are kept.
`MANIFEST-rfeval.tsv` has the bit scores and loci.

| structure | transcript (filtered GTF) | window | locus | bp |
|---|---|---|---|---|
| `URS000035D754_39947.17217__18S_rRNA` | `URS000035D754_39947.17217` | 1–1811 | 9:3626-5436 | 447 |
| `URS00006A9AE4_39947.17356__5.8S_rRNA` | `URS00006A9AE4_39947.17356` | 1–156 | 9:5637-5792 | 25 |
| `URS0000BFF4C1_39947.17383__25S_rRNA` | `URS0000BFF4C1_39947.17383` | 1–3383 | 9:6026-9408 | 872 |
| `URS0000C4C995_39947.4391__5S_rRNA` | `URS0000C4C995_39947.4391` | 1–120 | 11:12100779-12100898 | 34 |

The 18S, 5.8S and 25S are one full-length unit at the start of chromosome 9 (the
nucleolus organiser); 5S comes from the chromosome 11 5S cluster.

The names only match runs on
`release28/oryza_sativa/reference/oryza_sativa.IRGSP-1.0.filtered.gtf.gz`.

```
--gtf              release28/oryza_sativa/reference/oryza_sativa.IRGSP-1.0.filtered.gtf.gz \
--rfeval_reference release28/oryza_sativa/rfeval/rrna_common.db \
--rfeval_windows   release28/oryza_sativa/rfeval/rrna_common.windows.tsv
```
