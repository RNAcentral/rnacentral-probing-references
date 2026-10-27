# Rice (O. sativa IRGSP-1.0) — release 28

## Reference

- **GTF:** `reference/oryza_sativa.IRGSP-1.0.filtered.gtf.gz`, the RNAcentral release 28
  ncRNAs on IRGSP-1.0 after filtering (below).
- **FASTA:** not stored. The pipeline downloads the Ensembl Plants IRGSP-1.0 genome
  itself (`organism: oryza_sativa`).

Transcript IDs are `URS…_39947.<n>`, one ID per genomic copy.

## How the GTF was filtered

Input is `oryza_sativa.IRGSP-1.0.gff3.gz` from
`https://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/28.0/genome_coordinates/gff3/`.

```
scripts/rice/build_reference.sh 28   # from the repository root
```

It runs `scripts/filter_rnacentral_gff3.py --no-prefilter --pirna-max-loci 5`, then
`scripts/gff3_to_gtf_exons.py`, and writes `filter.log` beside this README.

`--no-prefilter` skips the trusted-database check, because rice's main sources
(Ensembl Plants, Rfam) are not on the per-RNA-type lists, which were written for
human. Main chromosomes are still enforced (numbered plus Mt/Pt). Rice has no
circRNAs or piRNAs, so step 0 and the piRNA cap do nothing.

| step | what it does | rice |
|---|---|---|
| 0 | drop circRNAs | none; 19,108 transcripts |
| 1 | same sequence listed several times at one locus: keep the best-ranked copy | 17,753 kept |
| 2 | overlapping records of one RNA type on one strand: keep one per RNAcentral predicted gene (or one per group) | **15,093 kept** |

Biotypes kept: pre_miRNA 8,749, ncRNA 2,337, miRNA 1,374, tRNA 892, snoRNA 773,
rRNA 594, snRNA 220, and smaller groups. The ncRNA group is labelled sRNA in the
release 27 GTF.

## rf-eval references

None yet.
