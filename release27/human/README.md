# Human (GRCh38) — release 27

## Reference

- **GTF:** `reference/homo_sapiens.GRCh38.filtered.gtf.gz`, the RNAcentral release 27
  ncRNAs on GRCh38 after filtering (below). It is gzipped, like every GTF here, and
  the pipeline reads `.gtf.gz` directly.
- **FASTA:** not stored. The pipeline downloads the Ensembl GRCh38 genome itself
  (`organism: homo_sapiens`).

Transcript IDs are `URS…_9606.<n>`, one ID per genomic copy. To link back to
RNAcentral, split on the last `.`.

## How the GTF was filtered

Input is `homo_sapiens.GRCh38.gff3.gz` from
`https://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/27.0/genome_coordinates/gff3/`.

```
scripts/human/build_reference.sh 27   # from the repository root
```

It runs `scripts/filter_rnacentral_gff3.py --pirna-max-loci 5`, then
`scripts/gff3_to_gtf_exons.py`, and writes `filter.log` beside this README.

| step | what it does | human |
|---|---|---|
| 0 | drop circRNAs (they need junction-aware mapping and folding) | 1,496,713 dropped; 889,995 left |
| piRNA cap | drop every entry of a piRNA sequence with more than 5 entries (repeat-derived, reads multimap) | 108,676 dropped |
| pre-filter | keep main chromosomes only; each record needs a trusted database for its RNA type (e.g. miRBase for miRNA), or a sequence in the Rfam seed alignments | 82,972 dropped (663 rescued by the Rfam list) |
| 1 | same sequence listed several times at one locus (once per database): keep the best-ranked copy; copies at other loci stay | 582,701 kept |
| 2 | overlapping records of one RNA type on one strand: keep one per RNAcentral predicted gene (or one per group); records with more than 5 databases always stay | **268,971 kept** |

Biotypes kept: lncRNA 234,869, piRNA 22,068, miRNA 3,106, snoRNA 2,106, pre_miRNA 1,963,
snRNA 1,886, and smaller groups.

The piRNA cap was added for human on 2026-10-08, and it removed only piRNAs.
Runs before that date (e.g. rnastruct000038/40) used the uncapped GTF, which had
105,323 piRNAs. Every other line is identical, so the rf-eval IDs below are still
in this GTF.

## rf-eval references

- `rfeval/rrna_common.*`: 28S, 18S and 5.8S, scored in every human run (below).

## Common rRNA set (`rfeval/rrna_common.*`)

The three nuclear rRNAs scored in every human run, taken unchanged from the best
sets of the rnastruct000038 and rnastruct000040 runs (byte-identical in both): the
best-covered rDNA copy of each molecule in those runs' RDATs;
only the structure names are changed: `<URS>__<molecule>` instead of `<URS>__<Rfam id>`.

| structure | transcript (filtered GTF) | window |
|---|---|---|
| `URS0000ABD879_9606.1467066__5.8S_rRNA` | `URS0000ABD879_9606.1467066` | 6573–6725 |
| `URS0000ABD879_9606.1467066__18S_rRNA` | `URS0000ABD879_9606.1467066` | 3632–5500 |
| `URS0000ABD87F_9606.1468181__28S_rRNA` | `URS0000ABD87F_9606.1468181` | 7927–12989 |

Both transcripts are in `release27/human/reference/homo_sapiens.GRCh38.filtered.gtf.gz`
(chr21 rDNA). The names only match runs on that GTF.

```
--rfeval_reference release27/human/rfeval/rrna_common.db \
--rfeval_windows   release27/human/rfeval/rrna_common.windows.tsv
```
