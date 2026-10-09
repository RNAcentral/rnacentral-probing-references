# Rotavirus A strain RF (rnastruct00024) — reference and Rfam structure for rf-eval

## Reference

`RVA_RF.fa` / `.gtf`: the eleven segment records, one transcript per segment,
the same KF729687–KF729697 accessions `conf/viral_genomes.config` resolves for
this strain, pinned locally (`build_rva_refs.py`, `MANIFEST.tsv`). RNAcentral
has no strain-RF transcript — two AU-repeat oligos from a polymerase structure,
and VP4 genes of other strains at 68–77 % identity — so as for influenza there
is no second record and no `_rnac` comparison.

## rf-eval reference

Rfam has a single rotavirus family, RF00501 Rota_CRE, the 3′ cis-acting
replication element of the NSP2 mRNA. `cmsearch --cut_ga` against the eleven
segments gives one hit above GA (40 bits): segment 8 (NSP2) at 992–1059, the
last 68 nt of the segment, 96.4 bits. The next best hit is 10 bits.

The structure is cmaligned back to the model with the seed alignment mapped in
(`--mapali --mapstr`) and projected; the family has no pseudoknot, so it is
17 nested pairs in 68 nt.

| structure | segment | KF729694.1 | nt | bp | Rfam | E |
|---|---|---|---|---|---|---|
| RVA_RF_NSP2_CRE | 8 (NSP2) | 992–1059 | 68 | 17 | RF00501 | 8e-25 |

Rebuild the reference: `cd release<N>/rotavirus && python3 ../../scripts/rotavirus/build_rva_refs.py` (Biopython, network) writes the genome FASTA/GTF to `work/`; then, from the root,
`python3 scripts/make_transcriptome_fa.py release<N>/rotavirus` cuts `reference/<stem>.transcriptome.fa` from them.
rf-eval rebuild: `scripts/rotavirus/build_rfeval_refs.py`, archived as run for release 27; it expects its inputs beside it, so adjust paths before rerunning (Docker, R2DT image for Infernal).

## Files

- `reference/RVA_RF.transcriptome.fa` — the eleven segments; `--fasta` with `--transcriptome`.
  `MANIFEST.tsv` lists them.
- `rfeval/rva_elements.db` — `--rfeval_reference`.
- `rfeval/rva_elements.windows.tsv` — `ref_seq_id start end structure_id`, 1-based
  inclusive on the segment. `--rfeval_windows`.
- `MANIFEST-rfeval.tsv` — coordinates, pair count, E-value.
- The Rfam models are in `scripts/rotavirus/cm/`; search tables, cut spans and alignments
  are not kept.

## Running

```
nextflow run . -profile <…> -resume \
  --transcriptome --fasta release27/rotavirus/reference/RVA_RF.transcriptome.fa \
  --rfeval_reference release27/rotavirus/rfeval/rva_elements.db \
  --rfeval_windows   release27/rotavirus/rfeval/rva_elements.windows.tsv
```
