# Human (GRCh38) — release 28

## Reference

- **GTF:** `reference/homo_sapiens.GRCh38.filtered.gtf.gz`, the RNAcentral release 28
  ncRNAs on GRCh38 after filtering (below). It is gzipped, like every GTF here, and
  the pipeline reads `.gtf.gz` directly.
- **FASTA:** not stored. The pipeline downloads the Ensembl GRCh38 genome itself
  (`organism: homo_sapiens`).
- **Small-RNA transcriptome:** `reference/homo_sapiens.small_rna.transcriptome.fa`, for
  size-selected libraries only (below). Its records are named `URS…_9606`, one per sequence.

Transcript IDs are `URS…_9606.<n>`, one ID per genomic copy. To link back to
RNAcentral, split on the last `.`.

## How the GTF was filtered

Input is `homo_sapiens.GRCh38.gff3.gz` from
`https://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/28.0/genome_coordinates/gff3/`.

```
scripts/homo_sapiens/build_reference.sh 28   # from the repository root
```

It runs `scripts/filter_rnacentral_gff3.py --pirna-max-loci 5`, then
`scripts/gff3_to_gtf_exons.py`, and writes `filter.log` beside this README.

| step | what it does | human |
|---|---|---|
| 0 | drop circRNAs (they need junction-aware mapping and folding) | 1,395,740 dropped; 1,631,638 left |
| piRNA cap | drop every entry of a piRNA sequence with more than 5 entries (repeat-derived, reads multimap) | 743,975 dropped |
| pre-filter | keep main chromosomes only; each record needs a trusted database for its RNA type (e.g. miRBase for miRNA), or a sequence in the Rfam seed alignments | 123,446 dropped (658 rescued by the Rfam list) |
| 1 | same sequence listed several times at one locus (once per database): keep the best-ranked copy; copies at other loci stay | 720,107 kept |
| 2 | overlapping records of one RNA type on one strand: keep one per RNAcentral predicted gene (or one per group); records with more than 5 databases always stay | **310,225 kept** |

Biotypes kept: lncRNA 240,085, piRNA 44,115, ncRNA 14,549, miRNA 3,111, snoRNA 2,107,
pre_miRNA 1,963, snRNA 1,891, and smaller groups.

The Rfam list was regenerated from release 28 (599 IDs, 598 of them also in the
release 27 list; Rfam 15.1). The IDs that left now pass the database check on
their own as ncRNA backed by Rfam.

Release 27 kept 268,971 transcripts. The growth is mainly RNAcentral content:
about 22,000 PirBase-only piRNAs now pass (184 in release 27; 93 % of them are
24–32 nt) and about 13,000 ENA/GeneCards ncRNA records. They are kept on
purpose, because PirBase is a curated set. Expect weak structure signal on the
piRNAs (short, PIWI-bound, mostly germline, size-selected libraries). They do not
take reads from host transcripts (STAR `AllBestScore`, no primary-only, MAPQ 0).

## rf-eval references

- `rfeval/rrna_common.*`: 28S, 18S and 5.8S, scored in every genome-route human run (below).
- `rfeval/small_rna.*`: every hit of 14 small-RNA families for runs on the small-RNA transcriptome (below).

## Common rRNA set (`rfeval/rrna_common.*`)

The three nuclear rRNAs scored in every human run, the same sets and windows as
release 27 (taken unchanged from the best sets of the rnastruct000038 and
rnastruct000040 runs): the best-covered rDNA copy of each molecule in those runs'
RDATs; only the structure names are changed: `<URS>__<molecule>` instead of
`<URS>__<Rfam id>`. Release 28 renumbered the copies, so `scripts/check_rfeval_ids.py`
re-pointed the windows to the new IDs, at the same loci.

| structure | transcript (filtered GTF) | window |
|---|---|---|
| `URS0000ABD879_9606.1785618__5.8S_rRNA` | `URS0000ABD879_9606.1785618` | 6573–6725 |
| `URS0000ABD879_9606.1785618__18S_rRNA` | `URS0000ABD879_9606.1785618` | 3632–5500 |
| `URS0000ABD87F_9606.1800971__28S_rRNA` | `URS0000ABD87F_9606.1800971` | 7927–12989 |

Both transcripts are in `release28/homo_sapiens/reference/homo_sapiens.GRCh38.filtered.gtf.gz`
(chr21 rDNA). The names only match runs on that GTF.

```
--gtf              release28/homo_sapiens/reference/homo_sapiens.GRCh38.filtered.gtf.gz \
--rfeval_reference release28/homo_sapiens/rfeval/rrna_common.db \
--rfeval_windows   release28/homo_sapiens/rfeval/rrna_common.windows.tsv
```

## Small-RNA transcriptome (`reference/homo_sapiens.small_rna.transcriptome.fa`)

For size-selected small-RNA libraries (metadata `rna_type: sRNA`: rnastruct00040, and
rnastruct00065 for its rerun). On the genome route, rnastruct00040's tRNA structures were
keyed on lncRNA records of 0.9–4.2 kb that span the tRNA genes (one a nuclear copy of
mitochondrial DNA on chr1), and every in vivo rf-eval structure was at chance. The
authors mapped to mature small-RNA sequences under ~200 nt (miRBase, GENCODE snoRNAs and
snRNAs, GtRNAdb tRNAs, RefSeq vault and Y RNAs, 5S rRNA), which this follows.

```
python3 scripts/build_small_rna_transcriptome.py homo_sapiens 28 [rnacentral_species_specific_ids.fasta.gz]
python3 scripts/rfeval_small_rna.py homo_sapiens 28
```

Built from this release's filtered GTF, with the same locus rule as the filter:

| step | records |
|---|---|
| GTF records in the kept classes (tRNA, snRNA, scaRNA, snoRNA, pre-miRNA, Y RNA, vault, SRP, RNase P/MRP, telomerase, and the listed rRNAs) | 7,300 |
| best record per locus (overlapping, lengths within 2×; database rank, then databases, then length) | 7,234 |
| one record per URS | 6,860 |
| longer than 500 nt, dropped (the listed rRNAs excepted) | −7 |
| within 2 substitutions of a better-ranked record, merged into it | −188 |
| **records** | **6,665** |

By class: 2,022 snoRNA, 1,886 pre-miRNA, 1,770 snRNA (7SK, URS0000AAF09C, is one), 630 Y RNA,
252 tRNA (19 mitochondrial), 58 SRP, 27 scaRNA, 12 vault, one each of RNase P, RNase MRP and
telomerase RNA, and five rRNAs: the two chr21 rDNA units `rrna_common` sits on
(URS0000ABD879, URS0000ABD87F), 5S (URS00000F9D45), mt 16S (URS000080E357) and mt 12S
(URS000044DFF6, class `ncRNA` in the GTF). Mature miRNAs, piRNAs, lncRNAs and the standalone
18S/5.8S/28S records (they would duplicate the units) are left out.
`MANIFEST-small_rna.tsv` lists each record's class, genomic copies, first locus and the
records merged into it.

- **Near-identical copies are merged.** 112 tRNAs were merged into 55 representatives. In a
  tRNA and rRNA test of this route, URS0000282AB2 scored below chance in vivo (AUROC 0.36):
  two isodecoders differ from it at 2 bases each, 3 of them paired, so their reads landed
  on it and the differences counted as mutations at paired bases. They are now merged into it.
- **No CCA tail.** Records stay identical to their URS sequence; Bowtie2 runs in local mode,
  so the 3 tail bases of a read are soft-clipped.
- **5.8S** is identical in both rDNA units, so its reads split between them.

## Small-RNA rf-eval set (`rfeval/small_rna.*`)

Every family below is scored on every record it hits, not on one representative, so a run is
judged on all the structures it covers: rf-eval skips records without a profile, and the class
medians are what to compare. Each record's best `cmsearch --cut_ga` hit is aligned with
`cmalign` and its `SS_cons` projected without pseudoknots; structures under 0.10 bp/nt are
dropped. 18S is the exception: its hits are the identical rDNA units, so only the best is
kept (`MANIFEST-rfeval-small_rna.tsv`).

| family | Rfam | records hit | structures |
|---|---|---|---|
| tRNA | RF00005 | 236 | 236 |
| 5S rRNA | RF00001 | 1 | 1 |
| 18S rRNA | RF01960 | 3 | 1 (URS0000ABD879, 3632–5500) |
| U1 / U2 / U4 / U5 snRNA | RF00003 / RF00004 / RF00015 / RF00020 | 90 / 150 / 99 / 31 | 90 / 125 / 96 / 31 |
| U6 snRNA | RF00026 | 1,224 | 0 (`SS_cons` under 0.10 bp/nt) |
| 7SL / 7SK | RF00017 / RF00100 | 54 / 137 | 54 / 136 |
| RNase P / RNase MRP | RF00009 / RF00030 | 1 / 1 | 1 / 1 |
| vault / Y RNA | RF00006 / RF00019 | 9 / 596 | 9 / 596 |

1,377 structures in all. Most snRNA, 7SK and Y RNA records are pseudogene copies with few or no
reads; they cost nothing when unprofiled, but a profiled pseudogene can carry reads from its
parent gene, so look at per-structure scores before reading much into a low class median. The
earlier set (one structure per family) scored only 5S and 18S in rnastruct00040, because the
chosen tRNA, U1 and 7SL records had no profile.

**Why this route.** The tRNA and rRNA test of rnastruct00040 (commit `dd2042a`) gave, as AUROC
(standard deviations above the random baseline):

| structure | in vivo | in vitro |
|---|---|---|
| tRNA URS000030BAD5 | 0.82 (+3.5) | 0.83 (+3.9) |
| mt tRNA URS00004206E2 | 0.76 (+3.3) | 0.80 (+3.7) |
| tRNA URS000042F13F | 0.66 (+1.8) | 0.79 (+3.8) |
| tRNA URS0000282AB2 | 0.36 (−1.9) | 0.54 (+0.6) |
| 18S | 0.72 (+2.8) | 0.65 (baseline not computed) |
| 5.8S | 0.53 (+0.3) | 0.67 (+2.4) |
| 28S | 0.66 (+4.4) | 0.63 (+4.0) |

5.8S at chance in vivo is expected: in the ribosome it pairs with 28S and is bound by
proteins, which its standalone structure does not include.

Run with `transcriptome: true`, `fasta` and `bowtie_all: false` (with `-a`, reads also
land on near-identical copies, where the differences count as mutations). The `rnacentral`
QC gate stays off until the mapped % on this reference has been checked.

```
--transcriptome --fasta release28/homo_sapiens/reference/homo_sapiens.small_rna.transcriptome.fa \
--bowtie_all false \
--rfeval_reference release28/homo_sapiens/rfeval/small_rna.db \
--rfeval_windows   release28/homo_sapiens/rfeval/small_rna.windows.tsv
```
