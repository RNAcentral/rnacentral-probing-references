# Yeast (S. cerevisiae R64-1-1) — release 28

## Reference

- **GTF:** `reference/saccharomyces_cerevisiae.R64-1-1.filtered.gtf.gz`, the RNAcentral
  release 28 ncRNAs on R64-1-1 after filtering (below).
- **FASTA:** not stored. The pipeline downloads the Ensembl R64-1-1 genome itself
  (`organism: saccharomyces_cerevisiae`).

Transcript IDs are `URS…_559292.<n>`, one ID per genomic copy.

## How the GTF was filtered

Input is `saccharomyces_cerevisiae.R64-1-1.gff3.gz` from
`https://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/28.0/genome_coordinates/gff3/`.

```
scripts/saccharomyces_cerevisiae/build_reference.sh 28   # from the repository root
```

It runs `scripts/filter_rnacentral_gff3.py --no-prefilter --pirna-max-loci 5`, then
`scripts/gff3_to_gtf_exons.py`, and writes `filter.log` beside this README.

`--no-prefilter` skips the trusted-database check, because yeast's main sources
(SGD, Rfam) are not on the per-RNA-type lists, which were written for human.
Main chromosomes are still enforced (roman numerals plus Mito). Yeast has no
circRNAs or piRNAs, so step 0 and the piRNA cap do nothing.

| step | what it does | yeast |
|---|---|---|
| 0 | drop circRNAs | none; 695 transcripts |
| 1 | same sequence listed several times at one locus: keep the best-ranked copy | 645 kept |
| 2 | overlapping records of one RNA type on one strand: keep one per RNAcentral predicted gene (or one per group) | **463 kept** |

Biotypes kept: tRNA 306, snoRNA 80, ncRNA 37, rRNA 21, misc_RNA 7, snRNA 6, and
single telomerase, SRP, RNase P, RNase MRP and lncRNA records plus one `other`.
The ncRNA group is labelled sRNA in the release 27 GTF.

## rf-eval references

- `rfeval/rrna_common.*`: 18S, 5.8S, 25S and 5S from rDNA unit 1 on chrXII (below).

## Common rRNA set (`rfeval/rrna_common.*`)

The four cytoplasmic rRNAs, all from rDNA unit 1 on chrXII (the two R64-1-1
copies are identical; 5S sits in the same repeat). Structures are the Rfam
`SS_cons` projected with `cmalign` onto `cmsearch --cut_ga` hits, trimmed to the
model-covered span. The sets and windows are those of release 27; release 28
renumbered the copies, so `scripts/check_rfeval_ids.py` re-pointed them to the
new IDs, at the same loci.

| structure | transcript (filtered GTF) | window |
|---|---|---|
| `URS00005F2C2D_559292.457__18S_rRNA` | `URS00005F2C2D_559292.457` | 1–1800 |
| `URS00001BF64E_559292.440__5.8S_rRNA` | `URS00001BF64E_559292.440` | 2–155 |
| `URS000061F377_559292.419__25S_rRNA` | `URS000061F377_559292.419` | 4–3396 |
| `URS000055688D_559292.469__5S_rRNA` | `URS000055688D_559292.469` | 1–121 |

The names only match runs on
`release28/saccharomyces_cerevisiae/reference/saccharomyces_cerevisiae.R64-1-1.filtered.gtf.gz`.
Built by `scripts/build_rrna_common.py yeast`, archived as run for release 27; it expects its inputs beside it, so adjust paths before rerunning (Docker, network).

```
--gtf              release28/saccharomyces_cerevisiae/reference/saccharomyces_cerevisiae.R64-1-1.filtered.gtf.gz \
--rfeval_reference release28/saccharomyces_cerevisiae/rfeval/rrna_common.db \
--rfeval_windows   release28/saccharomyces_cerevisiae/rfeval/rrna_common.windows.tsv
```
