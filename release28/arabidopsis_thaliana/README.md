# Arabidopsis (A. thaliana TAIR10) — release 28

## Reference

- **GTF:** `reference/arabidopsis_thaliana.TAIR10.filtered.gtf.gz`, the RNAcentral release 28
  ncRNAs on TAIR10 after filtering (below).
- **FASTA:** not stored. The pipeline downloads the Ensembl Plants TAIR10 genome itself
  (`organism: arabidopsis_thaliana`). Chromosome names (1–5, Mt, Pt) match that FASTA.

Transcript IDs are `URS…_3702.<n>`, one ID per genomic copy.

## How the GTF was filtered

Input is `arabidopsis_thaliana.TAIR10.gff3.gz` from
`https://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/28.0/genome_coordinates/gff3/`.

```
scripts/arabidopsis_thaliana/build_reference.sh 28   # from the repository root
```

It runs `scripts/filter_rnacentral_gff3.py --pirna-max-loci 5`, then
`scripts/gff3_to_gtf_exons.py`, and writes `filter.log` beside this README.
The pre-filter is on. Arabidopsis has no circRNAs or piRNAs in the GFF3, so step 0
and the piRNA cap do nothing.

| step | what it does | arabidopsis |
|---|---|---|
| 0 | drop circRNAs | none; 485,452 transcripts |
| pre-filter | keep main chromosomes only; each record needs a trusted database for its RNA type (e.g. miRBase for miRNA), or a sequence in the Rfam seed alignments | 110,485 dropped (264 rescued by the Rfam list) |
| 1 | same sequence listed several times at one locus (once per database): keep the best-ranked copy; copies at other loci stay | 237,076 kept |
| 2 | overlapping records of one RNA type on one strand: keep one per RNAcentral predicted gene (or one per group); records with more than 5 databases always stay | **131,036 kept** |

The pre-filter removes mostly ENA-only records. The Rfam list is a new block of 217
IDs, made from release 28 (Rfam 15.1).

Biotypes kept: siRNA 95,915, ncRNA 26,570, lncRNA 3,508, antisense_RNA 2,763,
tRNA 727, miRNA 425, snoRNA 372, pre_miRNA 338, rRNA 229, snRNA 93, and smaller groups.

## rf-eval references

- `rfeval/rrna_common.*`: 18S, 5.8S, 25S and 5S (below).

## Common rRNA set (`rfeval/rrna_common.*`)

Built by `scripts/rfeval_rrna_common.py arabidopsis_thaliana 28`: `cmsearch --cut_ga` on the
GTF's rRNA records (Ensembl REST sequences), one locus per molecule within 10 kb of
the SSU, anchored on the SSU whose unit holds the most linked molecules; 5S is
searched genome-wide. Structures are the Rfam `SS_cons` projected with `cmalign`,
without pseudoknots, and only structures of at least 0.10 bp/nt are kept.
`MANIFEST-rfeval.tsv` has the bit scores and loci.

| structure | transcript (filtered GTF) | window | locus | bp |
|---|---|---|---|---|
| `URS000059A818_3702.257375__18S_rRNA` | `URS000059A818_3702.257375` | 1–1808 | 3:14197677-14199484 | 447 |
| `URS00004C0A7E_3702.264163__5.8S_rRNA` | `URS00004C0A7E_3702.264163` | 6–161 | 3:14199753-14199916 | 25 |
| `URS00039548E4_3702.264909__25S_rRNA` | `URS00039548E4_3702.264909` | 2–3383 | 3:14200105-14203491 | 872 |
| `URS0000EFF782_3702.473488__5S_rRNA` | `URS0000EFF782_3702.473488` | 1–119 | 5:12820987-12828915 | 34 |

The names only match runs on
`release28/arabidopsis_thaliana/reference/arabidopsis_thaliana.TAIR10.filtered.gtf.gz`.

```
--gtf              release28/arabidopsis_thaliana/reference/arabidopsis_thaliana.TAIR10.filtered.gtf.gz \
--rfeval_reference release28/arabidopsis_thaliana/rfeval/rrna_common.db \
--rfeval_windows   release28/arabidopsis_thaliana/rfeval/rrna_common.windows.tsv
```
