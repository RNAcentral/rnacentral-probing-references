# RNAcentral release 27

References built from RNAcentral release 27. How they are made
is in the [repository README](../README.md); each organism's README has the
details and the exact pipeline options.

| organism | reference/ | rfeval/ |
|---|---|---|
| human | GRCh38 filtered GTF | `rrna_common` |
| mouse | GRCm39 filtered GTF | `rrna_common` (18S + 5S) |
| yeast | R64-1-1 filtered GTF | `rrna_common` |
| rice | IRGSP-1.0 filtered GTF | none yet |
| ecoli | transcriptome FASTA: 328 K-12 MG1655 records | `rrna_common` |
| dengue | transcriptome FASTA × 4 serotypes | `dengue_elements` |
| influenza | transcriptome FASTA: A/WSN/1933, 8 segments | `iav_elements` |
| rotavirus | transcriptome FASTA: strain RF, 11 segments | `rva_elements` |
| sarscov2 | transcriptome FASTA × 7 strains and variants | `sars2_elements` |
| zika | transcriptome FASTA × 3 strains | `zikv_elements` |

Eukaryote GTFs are gzipped (the pipeline reads `.gtf.gz`). Eukaryote FASTAs are
not stored: the pipeline downloads the Ensembl genome.
Bacteria and viruses run with `--transcriptome --fasta <stem>.transcriptome.fa` and no `--gtf`.

Notes for this release:

- **piRNA cap.** `--pirna-max-loci 5` was applied to every species from
  2026-10-08. Earlier human runs (e.g. rnastruct000038/40) used a GTF without it;
  only piRNA lines differ.
- **Rebuilt from the original scripts.** The filter scripts were reconstructed on
  2026-10-08 and reproduce all four GTFs byte-for-byte.
- **Left out.** Large rf-eval intermediates (Rfam CM databases, cmscan output,
  RDATs, alignments) and the unfiltered GTFs, which come from
  `https://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/27.0/genome_coordinates/gff3/`.
