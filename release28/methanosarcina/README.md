# M. acetivorans C2A reference: RNAcentral records placed on the NCBI genome

RNAcentral release 28 has no genome coordinates for M. acetivorans, so this is built the
viral way: every RNAcentral C2A record (strain taxid 188937 only) of
≥ 60 nt that sits in `NC_003552.1` at ≥ 99 % identity becomes a transcript. Where
records overlap on the same strand only the best by identity, then length, is kept; ties are broken by record ID.
The builder writes the genome and the kept records on it (`work/MACE_C2A.fa`,
`.gtf`). The reference, `reference/MACE_C2A.transcriptome.fa`, is the 65 records
cut from the genome, without the chromosome; it runs with `--transcriptome` and no
`--gtf`. The genome is never a transcript: a 5.8 Mb rf-fold never finishes.

| | records |
|---|---|
| RNAcentral C2A records ≥ 60 nt | 82 |
| ≥ 99 % identity to the genome | 81 |
| kept after overlap (`MANIFEST-genome.tsv`) | 65 |

Kept: mostly tRNAs, plus the rRNA records and a few other ncRNAs. RNAcentral stores a sequence once, so identical copies have no record of their own; their reads map to the identical record that is kept.

Only the strain taxid 188937 is used, not the species taxid.

```
--transcriptome --fasta release28/methanosarcina/reference/MACE_C2A.transcriptome.fa \
--rfeval_reference release28/methanosarcina/rfeval/rrna_common.db \
--rfeval_windows   release28/methanosarcina/rfeval/rrna_common.windows.tsv
```

Rebuild the reference: `cd release<N>/methanosarcina && python3 ../../scripts/methanosarcina/build_methanosarcina_refs.py` (Biopython, network) writes the genome FASTA/GTF to `work/`; then, from the root,
`python3 scripts/make_transcriptome_fa.py release<N>/methanosarcina` cuts `reference/<stem>.transcriptome.fa` from them. `work/records.fa` caches the
RNAcentral sequences; delete it to refetch. rnacentral.org throttles after ~60 quick
requests by serving its HTML home page, so a fetch from there runs at 1/s.

## Common rRNA set (`rfeval/rrna_common.*`)

The three rRNAs of one complete 16S/23S/5S operon, built by `scripts/rfeval_rrna_common.py methanosarcina 28`: `cmsearch --cut_ga` with the archaeal SSU/LSU CMs (RF01959, RF02540) and RF00001 for 5S on the records in the transcriptome FASTA, one locus per molecule within 10 kb of the SSU. Structures are the Rfam `SS_cons` projected with `cmalign`, without pseudoknots, and only structures of at least 0.10 bp/nt are kept. They are windowed on the URS records (RNAcentral keeps one record per distinct rRNA sequence, so not every operon has its own records). `MANIFEST-rfeval.tsv` has the bit scores and loci.

| structure | transcript | window | locus | bp |
|---|---|---|---|---|
| `URS0000835833_188937__16S_rRNA` | `URS0000835833_188937` | 1–1478 | chr:1073041-1074518 | 457 |
| `URS0000AB950E_188937__23S_rRNA` | `URS0000AB950E_188937` | 1–2904 | chr:1074708-1077611 | 781 |
| `URS00000DD7E9_188937__5S_rRNA` | `URS00000DD7E9_188937` | 5–127 | chr:1077681-1077815 | 34 |

For runs on `reference/MACE_C2A.transcriptome.fa` only.
