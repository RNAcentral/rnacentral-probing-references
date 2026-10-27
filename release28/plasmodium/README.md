# Plasmodium (P. falciparum 3D7, GCA000002765v3) — release 28

## Reference

- **GTF:** `reference/plasmodium_falciparum.GCA000002765v3.filtered.gtf.gz`, the RNAcentral
  release 28 ncRNAs on the 3D7 assembly GCA000002765v3 after filtering (below).
- **FASTA:** not stored. The pipeline downloads the Ensembl Protists genome itself
  (`organism: plasmodium_falciparum`). Chromosome names (`Pf3D7_01_v3` … `Pf3D7_14_v3`,
  `Pf3D7_API_v3`, `Pf3D7_MIT_v3`) match that FASTA.

Transcript IDs are `URS…_36329.<n>`, one ID per genomic copy. The metadata's taxid
36329 is run on the 3D7 assembly.

RNAcentral also publishes `plasmodium_falciparum.ASM276v2.gff3.gz`. It is not used:
its chromosome names do not match the Ensembl FASTA.

## How the GTF was filtered

Input is `plasmodium_falciparum.GCA000002765v3.gff3.gz` from
`https://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/28.0/genome_coordinates/gff3/`.

```
scripts/plasmodium/build_reference.sh 28   # from the repository root
```

It runs `scripts/filter_rnacentral_gff3.py --no-prefilter --pirna-max-loci 5`, then
`scripts/gff3_to_gtf_exons.py`, and writes `filter.log` beside this README.

`--no-prefilter` skips the trusted-database check, because the Ensembl
Protists-backed tRNA, rRNA and snoRNA records are not on the per-RNA-type lists,
which were written for human, and would be removed (a trial run with the check on
dropped 232 of the 1,641). Main chromosomes are still enforced. P. falciparum has
no circRNAs or piRNAs, so step 0 and the piRNA cap do nothing.

| step | what it does | plasmodium |
|---|---|---|
| 0 | drop circRNAs | none; 1,641 transcripts |
| 1 | same sequence listed several times at one locus: keep the best-ranked copy | 1,599 kept |
| 2 | overlapping records of one RNA type on one strand: keep one per RNAcentral predicted gene (or one per group) | **1,371 kept** |

Biotypes kept: ncRNA 793, antisense_RNA 365, tRNA 80, rRNA 73, snoRNA 53, snRNA 5,
and single SRP and RNase P records.

## rf-eval references

- `rfeval/rrna_common.*`: 18S, 5.8S, 28S and 5S (below).

## Common rRNA set (`rfeval/rrna_common.*`)

Built by `scripts/rfeval_rrna_common.py plasmodium 28`: `cmsearch --cut_ga` on the
GTF's rRNA records (Ensembl REST sequences), one locus per molecule within 10 kb of
the SSU, anchored on the SSU whose unit holds the most linked molecules; 5S is
searched genome-wide. Structures are the Rfam `SS_cons` projected with `cmalign`,
without pseudoknots, and only structures of at least 0.10 bp/nt are kept.
`MANIFEST-rfeval.tsv` has the bit scores and loci.

The unit is the A-type rDNA unit on chr5 (`Pf3D7_05_v3`); the 5S is from chr14.

| structure | transcript (filtered GTF) | window | locus | bp |
|---|---|---|---|---|
| `URS0000868510_36329.431__18S_rRNA` | `URS0000868510_36329.431` | 1–1607 | Pf3D7_05_v3:1289601-1291690 | 406 |
| `URS0000868556_36329.438__5.8S_rRNA` | `URS0000868556_36329.438` | 5–150 | Pf3D7_05_v3:1292053-1292207 | 24 |
| `URS0000290BFA_36329.441__28S_rRNA` | `URS0000290BFA_36329.441` | 4–3790 | Pf3D7_05_v3:1292410-1296199 | 847 |
| `URS000080E122_36329.1373__5S_rRNA` | `URS000080E122_36329.1373` | 1–118 | Pf3D7_14_v3:779377-779494 | 34 |

The names only match runs on
`release28/plasmodium/reference/plasmodium_falciparum.GCA000002765v3.filtered.gtf.gz`.

```
--gtf              release28/plasmodium/reference/plasmodium_falciparum.GCA000002765v3.filtered.gtf.gz \
--rfeval_reference release28/plasmodium/rfeval/rrna_common.db \
--rfeval_windows   release28/plasmodium/rfeval/rrna_common.windows.tsv
```
