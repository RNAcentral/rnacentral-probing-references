# Influenza A/WSN/1933(H1N1) (rnastruct00013) and A/Puerto Rico/8/1934(H1N1) (rnastruct00014/00134/00135) — references and Rfam structures for rf-eval

## References

| strain | datasets | stem | segments |
|---|---|---|---|
| A/WSN/1933 | 00013 | `IAV_WSN1933` | CY034132–CY034139 |
| A/Puerto Rico/8/1934 | 00014, 00134, 00135 | `IAV_PR8` | RefSeq NC_002016–NC_002023 |

`IAV_WSN1933.fa` / `IAV_PR8.fa` and their `.gtf`: the eight segment records of each strain, one transcript
per segment, the same accessions `conf/viral_genomes.config` resolves for these
strains, pinned locally (`build_iav_refs.py`, `MANIFEST.tsv`). RNAcentral has no
full-length influenza A transcript — only PDB promoter fragments of 12–70 nt — so
unlike dengue there is no second record and no `_rnac` comparison. WSN1933 is
byte-identical to release 27; PR8 is new.

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
is kept. The pipeline scores it, but runs need `rfeval_keep_pseudoknots = true`.

| structure | segment | accession | span | nt | bp | pk bp | Rfam | E |
|---|---|---|---|---|---|---|---|---|
| IAV_WSN1933_NS_PK | 8 (NS) | CY034136.1 | 498–545 | 48 | 5 | 6 | RF01099 | 6.8e-20 |
| IAV_PR8_NS_PK | 8 (NS) | NC_002020.1 | 523–570 | 48 | 5 | 6 | RF01099 | 1.1e-20 |

The PR8 structure was added in release 28 by `scripts/rfeval_viral_elements.py`
(`influenza` set, same family and labels) and is listed in `MANIFEST-rfeval-added.tsv`.

Rebuild the reference: `cd release<N>/influenza_a_virus && python3 ../../scripts/influenza_a_virus/build_iav_refs.py` (Biopython, network) writes the genome FASTA/GTF to `work/`; then, from the root,
`python3 scripts/make_transcriptome_fa.py release<N>/influenza_a_virus` cuts `reference/<stem>.transcriptome.fa` from them.
rf-eval rebuild: `scripts/influenza_a_virus/build_rfeval_refs.py`, archived as run for release 27; it expects its inputs beside it, so adjust paths before rerunning (Docker, R2DT image for Infernal).

## Files

- `reference/IAV_<strain>.transcriptome.fa` — the eight segments of one strain; `--fasta` with
  `--transcriptome`. `MANIFEST.tsv` lists them.
- `rfeval/iav_elements.db` — `--rfeval_reference`.
- `rfeval/iav_elements.windows.tsv` — `ref_seq_id start end structure_id`, 1-based
  inclusive on the segment. `--rfeval_windows`.
- `MANIFEST-rfeval.tsv` — coordinates, pair counts, E-value (WSN1933);
  `MANIFEST-rfeval-added.tsv` — the same for PR8.
- The Rfam models are in `scripts/influenza_a_virus/cm/`; search tables, cut spans and alignments
  are not kept.

## Running

```
nextflow run . -profile <…> -resume \
  --transcriptome --fasta release28/influenza_a_virus/reference/IAV_WSN1933.transcriptome.fa \
  --rfeval_reference release28/influenza_a_virus/rfeval/iav_elements.db \
  --rfeval_windows   release28/influenza_a_virus/rfeval/iav_elements.windows.tsv \
  --rfeval_keep_pseudoknots true
```

Swap `IAV_WSN1933` for `IAV_PR8` for rnastruct00014, 00134 and 00135; the window
on the other strain's segment is skipped with a warning.
