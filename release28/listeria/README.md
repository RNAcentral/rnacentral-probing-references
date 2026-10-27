# L. monocytogenes EGD-e reference: RNAcentral records placed on the NCBI genome

RNAcentral release 28 has no genome coordinates for L. monocytogenes, so this is built the
viral way: every RNAcentral EGD-e record (strain taxid 169963 only) of
≥ 60 nt that sits in `NC_003210.1` at ≥ 99 % identity becomes a transcript. Where
records overlap on the same strand only the best by identity, then length, is kept; ties are broken by record ID.
The builder writes the genome and the kept records on it (`work/LMON_EGDe.fa`,
`.gtf`). The reference, `reference/LMON_EGDe.transcriptome.fa`, is the 125 records
cut from the genome, without the chromosome; it runs with `--transcriptome` and no
`--gtf`. The genome is never a transcript: a 2.9 Mb rf-fold never finishes.

| | records |
|---|---|
| RNAcentral EGD-e records ≥ 60 nt | 158 |
| ≥ 99 % identity to the genome | 157 |
| kept after overlap (`MANIFEST-genome.tsv`) | 125 |

Kept: tRNAs, sRNAs, riboswitches and leader RNAs, and the rRNA records. RNAcentral stores a sequence once, so identical copies have no record of their own; their reads map to the identical record that is kept.

Only the strain taxid 169963 is used, not the species taxid.

```
--transcriptome --fasta release28/listeria/reference/LMON_EGDe.transcriptome.fa \
--rfeval_reference release28/listeria/rfeval/rrna_common.db \
--rfeval_windows   release28/listeria/rfeval/rrna_common.windows.tsv
```

Rebuild the reference: `cd release<N>/listeria && python3 ../../scripts/listeria/build_listeria_refs.py` (Biopython, network) writes the genome FASTA/GTF to `work/`; then, from the root,
`python3 scripts/make_transcriptome_fa.py release<N>/listeria` cuts `reference/<stem>.transcriptome.fa` from them. `work/records.fa` caches the
RNAcentral sequences; delete it to refetch. rnacentral.org throttles after ~60 quick
requests by serving its HTML home page, so a fetch from there runs at 1/s.

## Common rRNA set (`rfeval/rrna_common.*`)

The three rRNAs of one complete 16S/23S/5S operon, built by `scripts/rfeval_rrna_common.py listeria 28`: `cmsearch --cut_ga` with the bacterial SSU/LSU CMs (RF00177, RF02541) and RF00001 for 5S on the records in the transcriptome FASTA, one locus per molecule within 10 kb of the SSU. Structures are the Rfam `SS_cons` projected with `cmalign`, without pseudoknots, and only structures of at least 0.10 bp/nt are kept. They are windowed on the URS records (RNAcentral keeps one record per distinct rRNA sequence, so not every operon has its own records). `MANIFEST-rfeval.tsv` has the bit scores and loci.

| structure | transcript | window | locus | bp |
|---|---|---|---|---|
| `URS00009C7867_169963__16S_rRNA` | `URS00009C7867_169963` | 4–1553 | chr:243487-245047 | 461 |
| `URS00019CA5B9_169963__23S_rRNA` | `URS00019CA5B9_169963` | 3–2934 | chr:245284-248219 | 846 |
| `URS00000EF410_169963__5S_rRNA` | `URS00000EF410_169963` | 1–116 | chr:242274-242389 | 34 |

For runs on `reference/LMON_EGDe.transcriptome.fa` only.
