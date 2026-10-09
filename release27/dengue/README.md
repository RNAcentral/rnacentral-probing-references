# Dengue (release27 EDEN strains) — Rfam reference structures for rf-eval

Known-structure references for the four `DENV*_*.fa` references, built the
same way as `tmp/RF-eval/Zika`: Infernal `cmsearch` with the Rfam flavivirus
CMs, the hit span cut out, `cmalign`ed back to its model, and the family
`SS_cons` projected onto it (`project_ss_cons.py`, copied from
`tmp/RF-eval/Human/rfam`). Only hits above the family GA threshold are kept
(`cmsearch --cut_ga`):
the xrRNA models (RF01415 SL-IV, RF03547 Flavi_xrRNA) top out at 29 bits on
dengue against a GA of 54.5/37, so SL-I/SL-II are deliberately absent rather
than scored against an unreliable reference. Hairpins the projection left
with a loop under 3 nt have their innermost pairs unpaired.

Rebuild the reference: `cd release<N>/dengue && python3 ../../scripts/dengue/build_dengue_refs.py` (Biopython, network) writes the genome FASTA/GTF to `work/`; then, from the root,
`python3 scripts/make_transcriptome_fa.py release<N>/dengue` cuts `reference/<stem>.transcriptome.fa` from them.
rf-eval rebuild: `scripts/dengue/build_rfeval_refs.py`, archived as run for release 27; it expects its inputs beside it, so adjust paths before rerunning (Docker, R2DT image for Infernal). Reads `MANIFEST-genome.tsv` for the serotype/accession/sfRNA mapping.

## Elements (per serotype, see MANIFEST.tsv)

| element | Rfam | scored on | note |
|---|---|---|---|
| SLA  | RF02340 DENV_SLA        | genome record | nested inside the 5UTR window — pick one when reporting |
| 5UTR | RF03546 Flavivirus-5UTR | genome record | SLA + SLB; DENV2 and DENV4 only — DENV1/DENV3 score just under the GA of 80 bits |
| cHP  | RF00617                 | genome record | capsid-coding hairpin |
| DB1, DB2 | RF00525 Flavivirus_DB | genome + sfRNA records | dumbbells |
| 3SL  | RF00185 Flavi_CRE       | genome + sfRNA records | 3′ stem-loop / CRE |

Structure IDs carry the record they are scored on: `_ncbi` is the genome
accession, `_rnac` the RNAcentral sfRNA transcript (`URS…_<taxid>`), matching
the two records in each `DENV*.gtf`. The 3′UTR elements appear twice, once
per record, so the same structure can be compared across the two alignments in
the rf-eval metrics table. The 5′ elements lie outside the sfRNA record and
exist as `_ncbi` only.

## Files

- `reference/DENV<n>_<strain>.transcriptome.fa` — the genome plus the sfRNA record and 5′/3′ UTR;
  `--fasta` with `--transcriptome`. `MANIFEST-genome.tsv` lists them.
- `rfeval/dengue_elements.db` — 34 structures (per serotype: 2–3 `_ncbi` 5′
  elements, and the 3 3′ elements on both records), one entry per window.
  `--rfeval_reference`.
- `rfeval/dengue_elements.windows.tsv` — `ref_seq_id start end structure_id`, 1-based
  inclusive in the record's own coordinates. `--rfeval_windows`.
- `MANIFEST.tsv` — element, family, genome and record coordinates, bp count,
  cmsearch E-value.
- The Rfam models are in `scripts/dengue/cm/`; search tables, cut spans and alignments
  are not kept.

## Running

```
nextflow run . -profile <…> -resume \
  --rfeval_reference release27/dengue/rfeval/dengue_elements.db \
  --rfeval_windows   release27/dengue/rfeval/dengue_elements.windows.tsv
```
