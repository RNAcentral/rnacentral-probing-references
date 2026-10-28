# Human (GRCh38) — release 28

## Reference

- **GTF:** `reference/homo_sapiens.GRCh38.filtered.gtf.gz`, the RNAcentral release 28
  ncRNAs on GRCh38 after filtering (below). It is gzipped, like every GTF here, and
  the pipeline reads `.gtf.gz` directly.
- **FASTA:** not stored. The pipeline downloads the Ensembl GRCh38 genome itself
  (`organism: homo_sapiens`).

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

- `rfeval/rrna_common.*`: 28S, 18S and 5.8S, scored in every human run (below).

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

## Small-RNA test transcriptome (`homo_sapiens.trna_rrna.*`)

A test of the transcriptome route for small-RNA libraries, made for rnastruct00040
(HEK293T icSHAPE-MaP). On the genome route its tRNA profiles disagreed with the
authors' and every in vivo rf-eval reference was at chance; the authors mapped to
mature small-RNA sequences instead. If this run fixes that, the reference grows to
the other small-RNA classes.

```
python3 scripts/homo_sapiens/build_trna_rrna_test.py 28   # from the repository root
```

- **tRNAs:** every tRNA in the filtered GTF, one record per URS (556 genomic
  copies, 367 sequences, 19 of them mitochondrial; the filter keeps 19 of the 22). Sequences come from RNAcentral's
  GtRNAdb export, the rest from rnacentral.org; `work/` caches both.
- **No CCA tail.** Records stay identical to their URS sequence; Bowtie2 runs in
  local mode, so the 3 tail bases of a read are soft-clipped.
- **rRNA:** the two chr21 rDNA units the rf-eval windows sit on (URS0000ABD879,
  URS0000ABD87F), 5S (URS00000F9D45), mt 16S (URS000080E357) and mt 12S
  (URS000044DFF6, the Ensembl/HGNC MT-RNR1 record, 954 nt). The standalone
  18S/5.8S/28S records are left out because they would duplicate the units. 12S is
  not in the filtered GTF, so it is added by URS; without it, 12S reads would be
  pushed onto other records by local alignment. 5.8S is identical in both units, so
  its reads split between them and its window gets about half the depth.
- **rf-eval** (`rfeval/trna_rrna.*`): the three `rrna_common` structures and four
  tRNA structures taken from RNAcentral for the rnastruct00040 run, not from its
  probing data (`scripts/homo_sapiens/rnastruct00040_trna.db`). They are re-keyed
  onto these records, not rebuilt. In release 27 they sat on lncRNA records of
  0.9–4.2 kb that span the tRNA genes, one of them a nuclear copy of mitochondrial
  DNA on chr1; here they are on the tRNAs themselves. The fourth tRNA is
  mitochondrial (URS00004206E2).

Run with `transcriptome: true`, `fasta` and `bowtie_all: false`: with `-a`,
Bowtie2 also reports worse alignments, and reads would land on near-identical tRNA
genes, where the differences count as mutations. `rnacentral: false` turns the QC
gate off for the test: reads from the small RNAs not in this reference (snRNAs,
snoRNAs, 7SK, 7SL, Y RNAs) leave the mapped % below the gate's 60 %. The checks
still show in MultiQC. See `params/rnastruct00040-params.yaml`.
