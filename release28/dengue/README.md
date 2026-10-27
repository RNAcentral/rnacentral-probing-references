# Dengue (release28 strains: four EDEN, WP75, 16681, H241, TVP/360) — Rfam reference structures for rf-eval

## References

One FASTA per strain (`build_dengue_refs.py`, `MANIFEST-genome.tsv`):

| strain | datasets | genome | RNAcentral sfRNA record (length, identity, genome span) |
|---|---|---|---|
| DENV1 EDEN2402 | 00046 | EU081230.1, 10 735 nt | URS00019FFBCB_11053 (416 nt, 100 %, 10320–10735) |
| DENV2 EDEN3295 | 00047 | EU081177.1, 10 723 nt | URS000221D037_11060 (425 nt, 99.29 %, 10299–10723) |
| DENV3 EDEN863 | 00048 | EU081190.1, 10 706 nt | URS0003C2DF58_11069 (418 nt, 99.76 %, 10290–10706) |
| DENV4 EDEN2270 | 00049 | GQ398256.1, 10 653 nt | URS00025332D5_11070 (359 nt, 99.72 %, 10295–10653) |
| DENV1 WP75 | 00113 | AY145121.1, 10 735 nt | URS0002211232_11053 (416 nt, 99.04 %, 10320–10735) |
| DENV2 16681 | 00114 | NC_001474.2, 10 723 nt | URS00030298E1_11060 (425 nt, 100 %, 10299–10723) |
| DENV4 H241 | 00115 | KR011349.2, 10 664 nt | URS0003BD2A38_11070 (357 nt, 99.72 %, 10307–10664) |
| DENV4 TVP/360 | 00116 | KU513442.1, 10 649 nt | URS000252FE21_11070 (358 nt, 100 %, 10292–10649) |

The first four are the release 27 strains; the last four are new. Each FASTA holds
the genome, the sfRNA record and the `<accession>_5UTR`/`_3UTR` transcripts. The
DV4-H241 metadata names "DV4-H24, probably H241" and gives no GenBank accession.
KR011349.2 is the FDA/WHO reference-reagent H241 genome (Anez et al. 2016).

Against release 27: EDEN2402, EDEN863 and EDEN2270 are byte-identical. EDEN3295's
sfRNA record changed from `URS000287775D_11060` to `URS000221D037_11060`: a tie at
0.9929 identity, same coordinates (10299–10723), broken by record ID. Its three
`_rnac` windows were re-pointed to the new record, and both records match the
genome over every window.

## rf-eval reference

Known-structure references for the eight `DENV*_*.fa` references, built the
same way as `tmp/RF-eval/Zika`: Infernal `cmsearch` with the Rfam flavivirus
CMs, the hit span cut out, `cmalign`ed back to its model, and the family
`SS_cons` projected onto it (`project_ss_cons.py`, copied from
`tmp/RF-eval/Human/rfam`). Only hits above the family GA threshold are kept
(`cmsearch --cut_ga`):
the xrRNA models (RF01415 SL-IV, RF03547 Flavi_xrRNA) top out at 29 bits on
dengue against a GA of 54.5/37, so SL-I/SL-II are deliberately absent rather
than scored against an unreliable reference. Hairpins the projection left
with a loop under 3 nt have their innermost pairs unpaired.

The four release 27 strains' structures are carried over; the four new strains
were added with the same families and labels by `scripts/rfeval_viral_elements.py`
(`dengue` set, `cmscan --cut_ga`, 0.10 bp/nt minimum) and are listed in
`MANIFEST-rfeval-added.tsv`.

Rebuild the reference: `cd release<N>/dengue && python3 ../../scripts/dengue/build_dengue_refs.py` (Biopython, network) writes the genome FASTA/GTF to `work/`; then, from the root,
`python3 scripts/make_transcriptome_fa.py release<N>/dengue` cuts `reference/<stem>.transcriptome.fa` from them.
rf-eval rebuild: `scripts/dengue/build_rfeval_refs.py`, archived as run for release 27; it expects its inputs beside it, so adjust paths before rerunning (Docker, R2DT image for Infernal). Reads `MANIFEST-genome.tsv` for the serotype/accession/sfRNA mapping.

## Elements (per strain, see MANIFEST.tsv and MANIFEST-rfeval-added.tsv)

| element | Rfam | scored on | note |
|---|---|---|---|
| SLA  | RF02340 DENV_SLA        | genome record | nested inside the 5UTR window — pick one when reporting |
| 5UTR | RF03546 Flavivirus-5UTR | genome record | SLA + SLB; EDEN3295, EDEN2270, WP75 and 16681 only — the other four score under the GA of 80 bits |
| cHP  | RF00617                 | genome record | capsid-coding hairpin |
| DB1, DB2 | RF00525 Flavivirus_DB | genome + sfRNA records | dumbbells |
| 3SL  | RF00185 Flavi_CRE       | genome + sfRNA records | 3′ stem-loop / CRE |

Structure IDs carry the record they are scored on: `_ncbi` is the genome
accession, `_rnac` the RNAcentral sfRNA transcript (`URS…_<taxid>`), matching
the two records in each `DENV*.gtf`. The 3′UTR elements appear twice, once
per record, so the same structure can be compared across the two alignments in
the rf-eval metrics table. The 5′ elements lie outside the sfRNA record and
exist as `_ncbi` only. The release 27 strains' IDs are `DENV<n>_<element>`; the new
strains' are `DENV<n>_<strain>_<element>`.

## Files

- `reference/DENV<n>_<strain>.transcriptome.fa` — the genome plus the sfRNA record and 5′/3′ UTR;
  `--fasta` with `--transcriptome`. `MANIFEST-genome.tsv` lists them.
- `rfeval/dengue_elements.db` — 68 structures (34 from release 27, 34 for the new
  strains), one entry per window. `--rfeval_reference`.
- `rfeval/dengue_elements.windows.tsv` — `ref_seq_id start end structure_id`, 1-based
  inclusive in the record's own coordinates. `--rfeval_windows`.
- `MANIFEST.tsv` — release 27 strains: element, family, genome and record coordinates, bp count,
  cmsearch E-value.
- `MANIFEST-rfeval-added.tsv` — the same for the four new strains.
- The Rfam models are in `scripts/dengue/cm/`; search tables, cut spans and alignments
  are not kept.

## Running

```
nextflow run . -profile <…> -resume \
  --transcriptome --fasta release28/dengue/reference/DENV1_EDEN2402.transcriptome.fa \
  --rfeval_reference release28/dengue/rfeval/dengue_elements.db \
  --rfeval_windows   release28/dengue/rfeval/dengue_elements.windows.tsv
```

Swap the FASTA for the strain of the dataset (`DENV2_EDEN3295`, `DENV3_EDEN863`,
`DENV4_EDEN2270`, `DENV1_WP75`, `DENV2_16681`, `DENV4_H241`, `DENV4_TVP360`);
windows on other strains' genomes are skipped with a warning.
