# Mouse (GRCm39) — release 28

## Reference

- **GTF:** `reference/mus_musculus.GRCm39.filtered.gtf.gz`, the RNAcentral release 28
  ncRNAs on GRCm39 after filtering (below).
- **FASTA:** not stored. The pipeline downloads the Ensembl GRCm39 genome itself
  (`organism: mus_musculus`).

Transcript IDs are `URS…_10090.<n>`, one ID per genomic copy. To link back to
RNAcentral, split on the last `.`.

## How the GTF was filtered

Input is `mus_musculus.GRCm39.gff3.gz` from
`https://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/28.0/genome_coordinates/gff3/`.

```
scripts/mouse/build_reference.sh 28   # from the repository root
```

It runs `scripts/filter_rnacentral_gff3.py --pirna-max-loci 5`, then
`scripts/gff3_to_gtf_exons.py`, and writes `filter.log` beside this README.

| step | what it does | mouse |
|---|---|---|
| 0 | drop circRNAs (they need junction-aware mapping and folding) | 414,929 dropped; 8,307,015 left |
| piRNA cap | drop every entry of a piRNA sequence with more than 5 entries (repeat-derived, reads multimap) | 7,767,589 dropped |
| pre-filter | keep main chromosomes only; each record needs a trusted database for its RNA type (e.g. miRBase for miRNA), or a sequence in the Rfam seed alignments | 180,257 dropped (442 rescued by the Rfam list) |
| 1 | same sequence listed several times at one locus (once per database): keep the best-ranked copy; copies at other loci stay | 314,899 kept |
| 2 | overlapping records of one RNA type on one strand: keep one per RNAcentral predicted gene (or one per group); records with more than 5 databases always stay | **111,716 kept** |

Biotypes kept: lncRNA 63,827, piRNA 38,061, miRNA 3,704, snRNA 1,719, snoRNA 1,498,
pre_miRNA 1,474, tRNA 473, and smaller groups.

The Rfam list was regenerated from release 28 (386 IDs, 384 of them also in the
release 27 list; Rfam 15.1).

The piRNA cap matters most in mouse: it drops 7.77M of the 8.31M entries left after
the circRNA step. Release 27 kept 91,700 transcripts; the growth is mainly
RNAcentral content (piRNAs 24,218 in release 27, 38,061 now). piRNAs are kept on
purpose, because PirBase is a curated set. Expect weak structure signal on them,
and they do not take reads from host transcripts (STAR `AllBestScore`, no
primary-only, MAPQ 0).

## rf-eval references

- `rfeval/rrna_common.*`: 18S and 5S only. GRCm39 does not assemble the rDNA arrays, so
  28S and 5.8S cannot be placed on the genome (below). It is a
  placeholder: rebuild it from the first mouse run's RDATs by picking
  the best-covered locus per molecule, as was done for human.

## Common rRNA set (`rfeval/rrna_common.*`)

18S and 5S only. GRCm39 does not assemble the rDNA arrays: its one copy
(chr17:40153888–40159679, Ensembl `Rn18s-rs5`) is 5′ ETS, 18S and the start of
ITS1, with no gap after it. RNAcentral has the 28S (`URS000004B853_10090`) and
5.8S (`URS0000436019_10090`) but cannot place them on the genome, so no run on
this GTF measures them. The only other LSU/5.8S GA hits are mt 16S and a
5.8S-like chr18 pseudogene (Gm54867); both are excluded by requiring the rDNA
unit. 5S is the best-scoring copy of the chr8 Rn5s cluster (identical repeats).
Structures are the Rfam `SS_cons` projected with `cmalign` onto `cmsearch --cut_ga` hits.
The sets and windows are those of release 27; release 28 renumbered the copies, so
`scripts/check_rfeval_ids.py` re-pointed them to the new IDs, at the same loci.

| structure | transcript (filtered GTF) | window |
|---|---|---|
| `URS00006DD47C_10090.3638184__18S_rRNA` | `URS00006DD47C_10090.3638184` | 1–1849 |
| `URS000064027A_10090.7591433__5S_rRNA` | `URS000064027A_10090.7591433` | 1–119 |

The names only match runs on
`release28/mouse/reference/mus_musculus.GRCm39.filtered.gtf.gz`.
Built by `scripts/build_rrna_common.py mouse`, archived as run for release 27; it expects its inputs beside it, so adjust paths before rerunning (Docker, network).

```
--gtf              release28/mouse/reference/mus_musculus.GRCm39.filtered.gtf.gz \
--rfeval_reference release28/mouse/rfeval/rrna_common.db \
--rfeval_windows   release28/mouse/rfeval/rrna_common.windows.tsv
```
