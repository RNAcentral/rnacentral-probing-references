# E. coli K-12 MG1655 reference: RNAcentral records placed on the NCBI genome

RNAcentral release 27 has no genome coordinates for E. coli, so this is built the
viral way: every RNAcentral K-12 record (taxids 511145 MG1655 and 83333 K-12) of
≥ 60 nt that sits in `NC_000913.3` at ≥ 99 % identity becomes a transcript. Where
records overlap on the same strand only the best by identity, then length, is kept.
The builder writes the genome and the kept records on it (`work/ECOLI_MG1655.fa`,
`.gtf`). The reference, `reference/ECOLI_MG1655.transcriptome.fa`, is the 328 records
cut from the genome, without the chromosome; it runs with `--transcriptome` and no
`--gtf`. The chromosome is never a transcript: a 4.6 Mb rf-fold never finishes.

| | records |
|---|---|
| RNAcentral K-12 records ≥ 60 nt | 886 |
| ≥ 99 % identity to MG1655 | 770 |
| kept after overlap (`MANIFEST-genome.tsv`) | 328, all 100 % identical |

Kept: Rfam-derived SSU/LSU rRNA records for the operons (the mature 16S/23S
records overlap them and lose on length), 5S, tRNAs, and sRNAs incl. the CssrE and
RtT repeats. RNAcentral stores a sequence once, so identical copies (e.g. rrnE 16S,
four of the eight 5S genes) have no record of their own; their reads map to the identical
record that is kept.

Rejected below 99 %: tRNAs and PDB constructs one or two bases off every genomic
copy (a mismatching record would turn true reads into MaP mutations at that base),
and "MG1655 partial 16S/23S" records that are only 56–80 % identical to MG1655 rRNA.

```
--transcriptome --fasta release27/ecoli/reference/ECOLI_MG1655.transcriptome.fa \
--rfeval_reference release27/ecoli/rfeval/rrna_common.db \
--rfeval_windows   release27/ecoli/rfeval/rrna_common.windows.tsv
```

Rebuild the reference: `cd release<N>/ecoli && python3 ../../scripts/ecoli/build_ecoli_refs.py` (Biopython, network) writes the genome FASTA/GTF to `work/`; then, from the root,
`python3 scripts/make_transcriptome_fa.py release<N>/ecoli` cuts `reference/<stem>.transcriptome.fa` from them. `work/records.fa` caches the
RNAcentral sequences; delete it to refetch. rnacentral.org throttles after ~60 quick
requests by serving its HTML home page, so the fetch runs at 1/s (~20 min).

## Common rRNA set (`rfeval/rrna_common.*`)

The three rRNAs of one *rrn* operon (rrnC, the best-scoring forward-strand copy)
on K-12 MG1655 `NC_000913.3`, found on the genome and then placed on the RNAcentral
records in `reference/ECOLI_MG1655.gtf` that hold the same sequence (rrnC's 5S is identical to
the kept `URS0000049E57_511145`). For runs on `reference/ECOLI_MG1655.fa` + `reference/ECOLI_MG1655.gtf` only.
Structures are the Rfam `SS_cons` projected with `cmalign` onto `cmsearch --cut_ga`
hits, trimmed to the model-covered span.

| structure | transcript | window |
|---|---|---|
| `NC_000913.3__16S_rRNA` | `URS0000801AEA_511145` | 3–1544 |
| `NC_000913.3__23S_rRNA` | `URS0000763AE7_511145` | 10–2913 |
| `NC_000913.3__5S_rRNA` | `URS0000049E57_511145` | 3–118 |

Not for runs on `ecoli_calibration/ecoli_rrna_collab.fa`: the collaborator's 16S/23S
are not MG1655 (6 substitutions in 16S, a 1 nt indel in 23S). Those runs keep
`16S_Ecoli.db` / `23S_Ecoli.db`.
Built by `scripts/build_rrna_common.py ecoli`, archived as run for release 27; it expects its inputs beside it, so adjust paths before rerunning (Docker, network).

```
--rfeval_reference release27/ecoli/rfeval/rrna_common.db \
--rfeval_windows   release27/ecoli/rfeval/rrna_common.windows.tsv
```
