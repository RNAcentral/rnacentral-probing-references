# Rice (O. sativa IRGSP-1.0) — release 27

## Reference

- **GTF:** `reference/oryza_sativa.IRGSP-1.0.filtered.gtf.gz`, the RNAcentral release 27
  ncRNAs on IRGSP-1.0 after filtering (below).
- **FASTA:** not stored. The pipeline downloads the Ensembl Plants IRGSP-1.0 genome
  itself (`organism: oryza_sativa`).

## How the GTF was filtered

Input is `oryza_sativa.IRGSP-1.0.gff3.gz` from
`https://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/27.0/genome_coordinates/gff3/`.

```
scripts/rice/build_reference.sh 27   # from the repository root
```

It runs `scripts/filter_rnacentral_gff3.py --no-prefilter --pirna-max-loci 5`, then
`scripts/gff3_to_gtf_exons.py`, and writes `filter.log` beside this README.

`--no-prefilter` skips the trusted-database check, because rice's main sources
(Ensembl Plants, Rfam) are not on the per-RNA-type lists, which were written for
human. Main chromosomes are still enforced (numbered plus Mt/Pt). Rice has no
circRNAs or piRNAs, so step 0 and the piRNA cap do nothing.

| step | what it does | rice |
|---|---|---|
| 0 | drop circRNAs | none; 16,006 transcripts |
| 1 | same sequence listed several times at one locus: keep the best-ranked copy | 15,221 kept |
| 2 | overlapping records of one RNA type on one strand: keep one per RNAcentral predicted gene (or one per group) | **13,869 kept** |

Biotypes kept: pre_miRNA 8,739, sRNA 2,288, tRNA 759, snoRNA 727, miRNA 590,
rRNA 495, snRNA 124, and smaller groups.

## rf-eval references

None yet.
