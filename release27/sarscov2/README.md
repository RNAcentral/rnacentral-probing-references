# SARS-CoV-2 (rnastruct00025/00028/00030 Leiden-0002, rnastruct00031 USA-WA1/2020, rnastruct00072/00101–00104 GSE279203 variants) — references and Rfam structures for rf-eval

## References

One FASTA/GTF per strain (`build_sars_refs.py`, `MANIFEST-genome.tsv`):

| strain | datasets | genome | RNAcentral records in the GTF |
|---|---|---|---|
| Leiden0002 | 00025 (SHAPE), 00028 (DMS), 00030 (DMS) | MT510999.1, 29 890 nt | URS0002915722 (SL5, 145 nt, 99.3 %), URS00021ED9B3 (FSE, 88 nt, 100 %) |
| USAWA1 | 00031 (SHAPE) | MT576563.1, 29 850 nt | URS0002915722 (SL5, 100 %), URS00021ED9B3 (FSE, 100 %) |
| WT | 00072 (SHAPE) | NC_045512.2, 29 903 nt | SL5 (100 %), FSE (100 %) |
| Alpha | 00101 (SHAPE) | EPI_ISL_754083, 29 869 nt | SL5 (99.3 %), FSE (100 %) |
| Beta | 00102 (SHAPE) | EPI_ISL_1173248, 29 789 nt | FSE (100 %) |
| Delta | 00103 (SHAPE) | EPI_ISL_2621925, 29 840 nt | FSE (100 %) |
| OmicronXBB | 00104 (SHAPE) | EPI_ISL_14917728, 29 847 nt | SL5 (99.3 %), FSE (100 %) |

Leiden0002 and USAWA1 use the same accessions `conf/viral_genomes.config`
resolves, pinned locally. GSE279203's five variants were mapped by the authors
to GISAID assemblies. Its WT is Wuhan-Hu-1 base for base, so NC_045512.2 is
used; the other four exist only in GISAID and are read from the per-base
sequence in the authors' ShapeMapper `.map` files on GEO. Those genomes are the
GISAID consensus as deposited: Alpha keeps 4 N and 1 R (masked by the pipeline
before rf-count), Beta ends 92 nt short of the Wuhan-Hu-1 3′ end, and Delta's
leading 28-nt Illumina read-2 adapter is trimmed off. XBB lacks the s2m
(26-nt deletion). RNAcentral's SARS-CoV-2 records are PDB constructs and ENA
miRNAs; those of ≥ 60 nt at ≥ 99 % identity to the genome become transcripts
of their own (dengue's sfRNA rule) so an element can be scored on both. The
FASTA is the genome only; the GTF places each record on it at its aligned span.
Where records overlap only the best by identity then length is kept: the
124-nt SL5 (8UYS) inside the 145-nt SL5 (8QO5), and two 217/220-nt
ribosome-complex mRNAs (7O7Z/7O80) around the exact 88-nt FSE (6XRZ). The
SL5 record carries the USA-WA1 base at one position where Leiden differs; cut
from the genome, the transcript carries Leiden's. `<accession>_5UTR`/`_3UTR`
transcripts sit either side of Wuhan-Hu-1's ORF1ab start and ORF10 stop (its
RefSeq UTR bounds), found by exact 25-mer in each genome since the GISAID ones
have no annotation; they match the GenBank CDS bounds of Leiden0002 and USAWA1.

## rf-eval reference

Rfam's six betacoronavirus families searched at GA (`cmsearch --cut_ga`);
five hit Leiden0002, USAWA1, WT and Alpha; the packaging signal RF00182 (an
MHV element) hits none. Beta has no 3UTR (its genome stops short of the model's
3′ end), Delta no s2m (misses GA), XBB no s2m (deleted). The table gives the
Leiden0002/USAWA1 spans; the other strains' are in `MANIFEST.tsv`.
Structures are cmaligned back with the seed mapped in (`--mapali --mapstr`)
so pseudoknots survive, projected as `()` nested / `[]` pseudoknot pairs,
and scored with `-kp`. Structure IDs carry the record: `_ncbi` the genome,
`_rnac` the RNAcentral record.

| element | Rfam | span (both strains) | nt | bp | pk | on |
|---|---|---|---|---|---|---|
| 5UTR | RF03117 bCoV-5UTR | 11–287 | 277 | 37 | – | genome |
| 5UTR_SL5 | RF03117, clipped to the SL5 record | 137–281 | 145 | 16 | – | genome + SL5 record |
| FSE | RF00507 Corona_FSE | 13456–13533 | 78 | 20 | 5 | genome + FSE record |
| 3UTR | RF03122 bCoV-3UTR | 29506–29857 / 29850 | 352 / 345 | 78 / 73 | – | genome |
| pk3 | RF00165 Corona_pk3 | 29590–29649 | 60 | 9 | 6 | genome |
| s2m | RF00164 s2m | 29714–29756 | 43 | 12 | – | genome |

`5UTR_SL5` is the 5UTR consensus cut to the SL5 record's span with any pair
crossing the boundary unpaired — Rfam's conserved-only model is sparse there
(16 bp; the cryo-EM structure behind the record has many more). pk3 and s2m
lie inside the 3UTR window, as SLA does inside the dengue 5UTR: pick one level
when reporting. WT, Alpha, Beta and XBB hit the whole 5′ end (1–300, 47 bp)
where Leiden0002/USAWA1 start at 11. 49 entries in total: 8 each for
Leiden0002, USAWA1, WT and Alpha, 5 for Beta and Delta (no SL5 record), 7 for XBB.

Rebuild the reference: `cd release<N>/sarscov2 && python3 ../../scripts/sarscov2/build_sars_refs.py` (Biopython, network) writes the genome FASTA/GTF to `work/`; then, from the root,
`python3 scripts/make_transcriptome_fa.py release<N>/sarscov2` cuts `reference/<stem>.transcriptome.fa` from them.
rf-eval rebuild: `scripts/sarscov2/build_rfeval_refs.py`, archived as run for release 27; it expects its inputs beside it, so adjust paths before rerunning (Docker, R2DT image for Infernal). Each `--max` search over 30 kb takes several minutes, about 35 min for all seven.

## Files

- `reference/SARS2_<strain>.transcriptome.fa` — the genome plus each GTF record (SL5, FSE,
  5′/3′ UTR); `--fasta` with `--transcriptome`. `MANIFEST-genome.tsv` lists the records.
- `rfeval/sars2_elements.db` — all strains. `--rfeval_reference`.
- `rfeval/sars2_elements.windows.tsv` — `ref_seq_id start end structure_id`, 1-based
  inclusive in the record's own coordinates. `--rfeval_windows`.
- `MANIFEST.tsv` — element, family, genome and record coordinates, pair
  counts, E-value.
- The Rfam models are in `scripts/sarscov2/cm/`; search tables, cut spans and alignments
  are not kept.

## Running

```
nextflow run . -profile <…> -resume \
  --transcriptome --fasta release27/sarscov2/reference/SARS2_Leiden0002.transcriptome.fa \
  --rfeval_reference release27/sarscov2/rfeval/sars2_elements.db \
  --rfeval_windows   release27/sarscov2/rfeval/sars2_elements.windows.tsv
```

Swap `Leiden0002` for `USAWA1` (00031), `WT` (00072), `Alpha` (00101),
`Beta` (00102), `Delta` (00103) or `OmicronXBB` (00104); windows on other
strains' genomes are skipped with a warning. The SL5/FSE record windows share
one ID across strains, so each run scores them once per strain that has them.
