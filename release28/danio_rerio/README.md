# Zebrafish (D. rerio GRCz11) — release 28

## Reference

- **GTF:** `reference/danio_rerio.GRCz11.filtered.gtf.gz`, the RNAcentral release 28
  ncRNAs on GRCz11 after filtering (below).
- **FASTA:** not stored. The pipeline downloads the Ensembl GRCz11 genome itself
  (`organism: danio_rerio`). Chromosome names (1–25, MT) match that FASTA.

Transcript IDs are `URS…_7955.<n>`, one ID per genomic copy.

## How the GTF was filtered

Input is `danio_rerio.GRCz11.gff3.gz` from
`https://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/28.0/genome_coordinates/gff3/`.

```
scripts/danio_rerio/build_reference.sh 28   # from the repository root
```

It runs `scripts/filter_rnacentral_gff3.py --pirna-max-loci 5`, then
`scripts/gff3_to_gtf_exons.py`, and writes `filter.log` beside this README.
The pre-filter is on. ALT_CTG haplotype sequences are dropped with the other
non-main chromosomes. Zebrafish has no circRNAs in the GFF3, so step 0 does nothing.

| step | what it does | zebrafish |
|---|---|---|
| 0 | drop circRNAs | none; 70,401 transcripts |
| piRNA cap | drop every entry of a piRNA sequence with more than 5 entries (repeat-derived, reads multimap) | 1,590 dropped |
| pre-filter | keep main chromosomes only; each record needs a trusted database for its RNA type (e.g. miRBase for miRNA), or a sequence in the Rfam seed alignments | 49,019 dropped (253 rescued by the Rfam list) |
| 1 | same sequence listed several times at one locus (once per database): keep the best-ranked copy; copies at other loci stay | 16,668 kept |
| 2 | overlapping records of one RNA type on one strand: keep one per RNAcentral predicted gene (or one per group); records with more than 5 databases always stay | **14,415 kept** |

The pre-filter removes mostly Rfam-only rRNA and tRNA genome hits and NONCODE
lncRNAs. The Rfam list is a new block of 208 IDs, made from release 28 (Rfam 15.1).

Biotypes kept: tRNA 7,424, lncRNA 2,840, rRNA 1,981, miRNA 737, snRNA 563,
pre_miRNA 406, snoRNA 255, ncRNA 91, SRP_RNA 78, piRNA 10, and smaller groups.

## rf-eval references

- `rfeval/rrna_common.*`: 18S, 5.8S, 28S and 5S (below).

## Common rRNA set (`rfeval/rrna_common.*`)

Built by `scripts/rfeval_rrna_common.py danio_rerio 28`: `cmsearch --cut_ga` on the
GTF's rRNA records (Ensembl REST sequences), one locus per molecule within 10 kb of
the SSU, anchored on the SSU whose unit holds the most linked molecules; 5S is
searched genome-wide. Structures are the Rfam `SS_cons` projected with `cmalign`,
without pseudoknots, and only structures of at least 0.10 bp/nt are kept.
`MANIFEST-rfeval.tsv` has the bit scores and loci.

| structure | transcript (filtered GTF) | window | locus | bp |
|---|---|---|---|---|
| `URS0000ABD7E7_7955.63204__18S_rRNA` | `URS0000ABD7E7_7955.63204` | 1–1887 | 5:824709-826807 | 447 |
| `URS00006D7FFF_7955.63381__5.8S_rRNA` | `URS00006D7FFF_7955.63381` | 1–154 | 5:836137-836290 | 25 |
| `URS000214A642_7955.63173__28S_rRNA` | `URS000214A642_7955.63173` | 4–3347 | 5:820792-824138 | 652 |
| `URS00006ED473_7955.41686__5S_rRNA` | `URS00006ED473_7955.41686` | 1–119 | 4:34612088-34612206 | 34 |

The 28S record is partial: 3,347 nt, against about 4 kb for the full-length molecule,
so the structure covers only the part the record holds.

The names only match runs on
`release28/danio_rerio/reference/danio_rerio.GRCz11.filtered.gtf.gz`.

```
--gtf              release28/danio_rerio/reference/danio_rerio.GRCz11.filtered.gtf.gz \
--rfeval_reference release28/danio_rerio/rfeval/rrna_common.db \
--rfeval_windows   release28/danio_rerio/rfeval/rrna_common.windows.tsv
```
