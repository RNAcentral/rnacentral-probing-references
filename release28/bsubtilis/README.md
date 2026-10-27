# B. subtilis 168 reference: RNAcentral records placed on the NCBI genome

RNAcentral release 28 has no genome coordinates for B. subtilis, so this is built the
viral way: every RNAcentral strain 168 record (strain taxid 224308 only) of
≥ 60 nt that sits in `NC_000964.3` at ≥ 99 % identity becomes a transcript. Where
records overlap on the same strand only the best by identity, then length, is kept; ties are broken by record ID.
The builder writes the genome and the kept records on it (`work/BSUB_168.fa`,
`.gtf`). The reference, `reference/BSUB_168.transcriptome.fa`, is the 182 records
cut from the genome, without the chromosome; it runs with `--transcriptome` and no
`--gtf`. The genome is never a transcript: a 4.2 Mb rf-fold never finishes.

| | records |
|---|---|
| RNAcentral strain 168 records ≥ 60 nt | 385 |
| ≥ 99 % identity to the genome | 366 |
| kept after overlap (`MANIFEST-genome.tsv`) | 182 |

Kept: rRNA records for the operons (SSU, LSU and 5S), tRNAs, riboswitches and T-box leaders, sRNAs, and ribosomal-protein leader RNAs. RNAcentral stores a sequence once, so identical copies have no record of their own; their reads map to the identical record that is kept.

Only the strain taxid 224308 is used, not the species taxid.

```
--transcriptome --fasta release28/bsubtilis/reference/BSUB_168.transcriptome.fa \
--rfeval_reference release28/bsubtilis/rfeval/rrna_common.db \
--rfeval_windows   release28/bsubtilis/rfeval/rrna_common.windows.tsv
```

Rebuild the reference: `cd release<N>/bsubtilis && python3 ../../scripts/bsubtilis/build_bsubtilis_refs.py` (Biopython, network) writes the genome FASTA/GTF to `work/`; then, from the root,
`python3 scripts/make_transcriptome_fa.py release<N>/bsubtilis` cuts `reference/<stem>.transcriptome.fa` from them. `work/records.fa` caches the
RNAcentral sequences; delete it to refetch. rnacentral.org throttles after ~60 quick
requests by serving its HTML home page, so a fetch from there runs at 1/s.

## Common rRNA set (`rfeval/rrna_common.*`)

The three rRNAs of one complete 16S/23S/5S operon, built by `scripts/rfeval_rrna_common.py bsubtilis 28`: `cmsearch --cut_ga` with the bacterial SSU/LSU CMs (RF00177, RF02541) and RF00001 for 5S on the records in the transcriptome FASTA, one locus per molecule within 10 kb of the SSU. Structures are the Rfam `SS_cons` projected with `cmalign`, without pseudoknots, and only structures of at least 0.10 bp/nt are kept. They are windowed on the URS records (RNAcentral keeps one record per distinct rRNA sequence, so not every operon has its own records). `MANIFEST-rfeval.tsv` has the bit scores and loci.

| structure | transcript | window | locus | bp |
|---|---|---|---|---|
| `URS0000803274_224308__16S_rRNA` | `URS0000803274_224308` | 3–1552 | chr:946696-948251 | 461 |
| `URS00007FE5A7_224308__23S_rRNA` | `URS00007FE5A7_224308` | 2–2929 | chr:948417-951349 | 844 |
| `URS00004E2288_224308__5S_rRNA` | `URS00004E2288_224308` | 1–116 | chr:951457-951572 | 34 |

For runs on `reference/BSUB_168.transcriptome.fa` only.
