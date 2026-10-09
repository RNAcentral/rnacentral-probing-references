# Influenza A/WSN/1933(H1N1) (rnastruct00013) — reference and Rfam structure for rf-eval

## Reference

`IAV_WSN1933.fa` / `.gtf`: the eight segment records, one transcript per
segment, the same CY034132–CY034139 accessions `conf/viral_genomes.config`
resolves for this strain, pinned locally (`build_iav_refs.py`,
`MANIFEST.tsv`). RNAcentral has no full-length influenza A transcript — only
PDB promoter fragments of 12–70 nt — so unlike dengue there is no second
record and no `_rnac` comparison.

## rf-eval reference

Rfam has a single influenza A family, RF01099 PK-IAV, a 48-nt pseudoknot in
the NS segment (Gultyaev 2007; Moss 2010; Pseudobase). `cmsearch --cut_ga`
against the eight segments gives one hit above GA (47 bits): segment 8 at
498–545, 77.5 bits. The next best hit, in NA, is 15 bits and is not kept.

The structure is cmaligned back to the model with the seed alignment mapped
in (`--mapali --mapstr`) so the pseudoknot annotation survives, then
projected: 5 nested pairs as `()`, 6 pseudoknot pairs as `[]`. The dengue and
human builds drop pseudoknots to stay comparable with R2DT; there is no R2DT
arm for influenza and for this family the pseudoknot *is* the structure, so it
is kept. The pipeline scores it (`rfeval_keep_pseudoknots = true`).

| structure | segment | CY034136.1 | nt | bp | pk bp | Rfam | E |
|---|---|---|---|---|---|---|---|
| IAV_WSN1933_NS_PK | 8 (NS) | 498–545 | 48 | 5 | 6 | RF01099 | 6.8e-20 |

Rebuild the reference: `cd release<N>/influenza && python3 ../../scripts/influenza/build_iav_refs.py` (Biopython, network) writes the genome FASTA/GTF to `work/`; then, from the root,
`python3 scripts/make_transcriptome_fa.py release<N>/influenza` cuts `reference/<stem>.transcriptome.fa` from them.
rf-eval rebuild: `scripts/influenza/build_rfeval_refs.py`, archived as run for release 27; it expects its inputs beside it, so adjust paths before rerunning (Docker, R2DT image for Infernal).

## Files

- `reference/IAV_WSN1933.transcriptome.fa` — the eight segments; `--fasta` with `--transcriptome`.
  `MANIFEST.tsv` lists them.
- `rfeval/iav_elements.db` — `--rfeval_reference`.
- `rfeval/iav_elements.windows.tsv` — `ref_seq_id start end structure_id`, 1-based
  inclusive on the segment. `--rfeval_windows`.
- `MANIFEST-rfeval.tsv` — coordinates, pair counts, E-value.
- The Rfam models are in `scripts/influenza/cm/`; search tables, cut spans and alignments
  are not kept.

## Running

```
nextflow run . -profile <…> -resume \
  --transcriptome --fasta release27/influenza/reference/IAV_WSN1933.transcriptome.fa \
  --rfeval_reference release27/influenza/rfeval/iav_elements.db \
  --rfeval_windows   release27/influenza/rfeval/iav_elements.windows.tsv
```
