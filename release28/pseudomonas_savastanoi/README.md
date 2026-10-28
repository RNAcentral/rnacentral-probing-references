# P. savastanoi pv. phaseolicola 1448A reference: RNAcentral records placed on the NCBI genome

RNAcentral release 28 has no genome coordinates for P. savastanoi, so this is built the
viral way: every RNAcentral 1448A record (strain taxid 264730 only) of
≥ 60 nt that sits in the chromosome `NC_005773.3` or its plasmids `NC_007274.1` and
`NC_007275.1` at ≥ 99 % identity becomes a transcript, placed on the replicon it matches. Where
records overlap on the same replicon and strand only the best by identity, then length, is kept; ties are broken by record ID.
The builder writes the genome and the kept records on it (`work/PSAV_1448A.fa`,
`.gtf`). The reference, `reference/PSAV_1448A.transcriptome.fa`, is the 79 records
cut from the genome, without the chromosome and plasmids; it runs with `--transcriptome` and no
`--gtf`. No replicon is ever a transcript: a 5.9 Mb rf-fold never finishes.

| | records |
|---|---|
| RNAcentral 1448A records ≥ 60 nt | 82 |
| ≥ 99 % identity to the genome | 82 |
| kept after overlap (`MANIFEST-genome.tsv`) | 79 |

All 79 kept records landed on the chromosome; none are on the plasmids. `MANIFEST-genome.tsv` has a `replicon` column.

Kept: tRNAs, sRNAs, antisense RNAs, riboswitches and the rRNA records. RNAcentral stores a sequence once, so identical copies have no record of their own; their reads map to the identical record that is kept.

Only the strain taxid 264730 is used, not the species taxid.

```
--transcriptome --fasta release28/pseudomonas_savastanoi/reference/PSAV_1448A.transcriptome.fa \
--rfeval_reference release28/pseudomonas_savastanoi/rfeval/rrna_common.db \
--rfeval_windows   release28/pseudomonas_savastanoi/rfeval/rrna_common.windows.tsv
```

Rebuild the reference: `cd release<N>/pseudomonas_savastanoi && python3 ../../scripts/pseudomonas_savastanoi/build_pseudomonas_refs.py` (Biopython, network) writes the genome FASTA/GTF to `work/`; then, from the root,
`python3 scripts/make_transcriptome_fa.py release<N>/pseudomonas_savastanoi` cuts `reference/<stem>.transcriptome.fa` from them. `work/records.fa` caches the
RNAcentral sequences; delete it to refetch. rnacentral.org throttles after ~60 quick
requests by serving its HTML home page, so a fetch from there runs at 1/s.

## Common rRNA set (`rfeval/rrna_common.*`)

The three rRNAs of one complete 16S/23S/5S operon, built by `scripts/rfeval_rrna_common.py pseudomonas_savastanoi 28`: `cmsearch --cut_ga` with the bacterial SSU/LSU CMs (RF00177, RF02541) and RF00001 for 5S on the records in the transcriptome FASTA, one locus per molecule within 10 kb of the SSU. Structures are the Rfam `SS_cons` projected with `cmalign`, without pseudoknots, and only structures of at least 0.10 bp/nt are kept. They are windowed on the URS records (RNAcentral keeps one record per distinct rRNA sequence, so not every operon has its own records). `MANIFEST-rfeval.tsv` has the bit scores and loci.

| structure | transcript | window | locus | bp |
|---|---|---|---|---|
| `URS0000887589_264730__16S_rRNA` | `URS0000887589_264730` | 2–1540 | NC_005773.3:804000-805539 | 455 |
| `URS000013225D_264730__23S_rRNA` | `URS000013225D_264730` | 1–2893 | NC_005773.3:806085-808986 | 830 |
| `URS000048AB4C_264730__5S_rRNA` | `URS000048AB4C_264730` | 3–118 | NC_005773.3:809136-809255 | 34 |

For runs on `reference/PSAV_1448A.transcriptome.fa` only.
