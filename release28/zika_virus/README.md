# Zika virus (rnastruct00009/00051 MR766, 00052 BeH815744, 00053/00117 H/PF/2013) — references and Rfam structures for rf-eval

## References

One FASTA/GTF per strain (`build_zika_refs.py`, `MANIFEST-genome.tsv`), the
accessions the probing-metadata YAMLs name and `conf/viral_genomes.config`
resolves:

| strain | datasets | genome | note |
|---|---|---|---|
| MR766 | 00009, 00051 | AY632535.2, 10 794 nt | complete genome |
| BeH815744 | 00052 | KU365780.1, 10 662 nt | **CDS only**: starts 14 nt into the 5′UTR, ends ~130 nt before the genome 3′ end |
| HPF2013 | 00053, 00117 | KJ776791.2, 10 807 nt | complete genome |

RNAcentral records are appended under the 99 % rule used for SARS-CoV-2. Release 28
adds one for HPF2013: `URS0004101ADB_64320` (sfRNA1, 413 nt, 100 %, 10395–10807 of
KJ776791.2). MR766 and BeH815744 have none: the Zika partial sfRNAs in RNAcentral
(ENA, 103–124 nt) come from another Asian-lineage isolate (≤ 97.5 % to HPF2013,
≤ 93 % to MR766) and the PDB xrRNA constructs (5TPY, 7U4A) are engineered
(≤ 90 %). The FASTA is the genome only. Besides the genome, the GTF has
`<accession>_5UTR`/`_3UTR` transcripts either side of the GenBank CDS, as for dengue
(BeH815744's are partial: 93 and 297 nt). The other two strains' FASTAs are
byte-identical to release 27; HPF2013's now also holds the sfRNA record.

## rf-eval reference

The seven flavivirus families searched at GA (`cmsearch --cut_ga`). Three
hit the complete genomes; 5UTR (RF03546, GA 80) and cHP (RF00617, GA 40)
score under their thresholds on Zika, and the xrRNA families (RF03547,
RF01415) do as on dengue. BeH815744 has no hit at all: SLA needs the first
14 nt the record lacks, and the 3′UTR elements lie past its end. Structures
are cmaligned back with the seed mapped in and projected as `()` nested /
`[]` pseudoknot pairs (none in these three families).

| element | Rfam | MR766 | HPF2013 | nt | bp |
|---|---|---|---|---|---|
| SLA | RF02340 DENV_SLA | 1–73 | 1–73 | 73 | 26 |
| DB | RF00525 Flavivirus_DB | 10605–10678 | 10618–10691 | 74 | 19 |
| 3SL | RF00185 Flavi_CRE | 10697–10794 | 10711–10807 | 98 / 97 | 23 / 24 |

6 entries, `_ncbi` only. The windows are keyed on the genome, so they carry over from
release 27 unchanged; the new HPF2013 sfRNA record has no windows of its own.

Rebuild the reference: `cd release<N>/zika_virus && python3 ../../scripts/zika_virus/build_zika_refs.py` (Biopython, network) writes the genome FASTA/GTF to `work/`; then, from the root,
`python3 scripts/make_transcriptome_fa.py release<N>/zika_virus` cuts `reference/<stem>.transcriptome.fa` from them.
rf-eval rebuild: `scripts/zika_virus/build_rfeval_refs.py`, archived as run for release 27; it expects its inputs beside it, so adjust paths before rerunning (Docker, R2DT image for Infernal).

## Files

- `reference/ZIKV_<strain>.transcriptome.fa` — the genome plus its 5′/3′ UTR (HPF2013 also the sfRNA record); `--fasta` with
  `--transcriptome`. `MANIFEST-genome.tsv` lists them.
- `rfeval/zikv_elements.db` — MR766 + HPF2013. `--rfeval_reference`.
- `rfeval/zikv_elements.windows.tsv` — `ref_seq_id start end structure_id`, 1-based
  inclusive. `--rfeval_windows`.
- `MANIFEST.tsv` — element, family, coordinates, pair counts, E-value.
- The Rfam models are in `scripts/zika_virus/cm/`; search tables, cut spans and alignments
  are not kept.

## Running

```
nextflow run . -profile <…> -resume \
  --transcriptome --fasta release28/zika_virus/reference/ZIKV_MR766.transcriptome.fa \
  --rfeval_reference release28/zika_virus/rfeval/zikv_elements.db \
  --rfeval_windows   release28/zika_virus/rfeval/zikv_elements.windows.tsv
```

Swap `MR766` for `HPF2013` for rnastruct00053 and rnastruct00117. For rnastruct00052 there is
nothing to score on KU365780.1; omit the two `rfeval_*` params.
