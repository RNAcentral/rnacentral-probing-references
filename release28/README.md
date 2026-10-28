# RNAcentral release 28

References built from RNAcentral release 28. How they are made
is in the [repository README](../README.md); each organism's README has the
details and the exact pipeline options.

| organism | reference/ | rfeval/ |
|---|---|---|
| homo_sapiens | GRCh38 filtered GTF | `rrna_common` |
| mus_musculus | GRCm39 filtered GTF | `rrna_common` (18S + 5S) |
| saccharomyces_cerevisiae | R64-1-1 filtered GTF | `rrna_common` |
| oryza_sativa | IRGSP-1.0 filtered GTF | `rrna_common` (new) |
| arabidopsis_thaliana | TAIR10 filtered GTF | `rrna_common` |
| danio_rerio | GRCz11 filtered GTF | `rrna_common` |
| plasmodium_falciparum | GCA000002765v3 (3D7) filtered GTF | `rrna_common` |
| triticum_turgidum | Svevo.v1 filtered GTF | `rrna_common` (18S + 26S + 5S, partial) |
| escherichia_coli | transcriptome FASTA: 328 K-12 MG1655 records | `rrna_common` |
| bacillus_subtilis | transcriptome FASTA: 182 strain 168 records | `rrna_common` |
| listeria_monocytogenes | transcriptome FASTA: 125 EGD-e records | `rrna_common` |
| methanosarcina_acetivorans | transcriptome FASTA: 65 C2A records | `rrna_common` (archaeal CMs) |
| pseudomonas_savastanoi | transcriptome FASTA: 79 1448A records | `rrna_common` |
| dengue_virus | transcriptome FASTA × 8 strains (4 serotypes) | `dengue_elements` |
| influenza_a_virus | transcriptome FASTA: A/WSN/1933 and A/Puerto Rico/8/1934, 8 segments each | `iav_elements` |
| rotavirus_a | transcriptome FASTA: strain RF, 11 segments | `rva_elements` |
| sars_cov_2 | transcriptome FASTA × 7 strains and variants | `sars2_elements` |
| zika_virus | transcriptome FASTA × 3 strains | `zikv_elements` |
| yellow_fever_virus | transcriptome FASTA × 3 strains (17D, Asibi, Dakar) | `yellow_fever_elements` |
| west_nile_virus | transcriptome FASTA: NY99 | `west_nile_elements` |
| murine_norovirus | transcriptome FASTA: murine norovirus MNV-1 CW3 | `norovirus_elements` |
| porcine_epidemic_diarrhea_virus | transcriptome FASTA: strain AJ1102 | `pedv_elements` |
| human_bocavirus | transcriptome FASTA: HBoV1 Salvador1 | none (no Rfam family at GA) |
| bluetongue_virus | transcriptome FASTA: BTV-1, 10 segments | none |
| hiv | transcriptome FASTA: NL4-3 provirus | `hiv_elements` |

Eukaryote GTFs are gzipped (the pipeline reads `.gtf.gz`). Eukaryote FASTAs are
not stored: the pipeline downloads the Ensembl genome.
Bacteria and viruses run with `--transcriptome --fasta <stem>.transcriptome.fa` and no `--gtf`.

Notes for this release:

- **Rebuilt from release 27 first.** The release 27 human, mouse, yeast and rice
  GTFs and `filter.log` files were rebuilt from the committed scripts and
  reproduced byte-for-byte before anything was changed.
- **piRNA cap.** `--pirna-max-loci 5` was applied to every species, including the
  new ones.
- **Rfam allowlist.** Regenerated for human (599 IDs) and mouse (386) from release 28;
  598 and 384 of them are shared with release 27. The IDs that left now pass the
  database check on their own, as ncRNA backed by Rfam. New blocks for
  Arabidopsis (217) and zebrafish (208). Rfam release used: 15.1.
- **Human and mouse grew** mainly from RNAcentral content. About 22,000 PirBase-only
  piRNAs now pass the filter (184 in release 27; 93 % of them are 24–32 nt), and
  human has about 13,000 ENA/GeneCards ncRNA records. They are kept on purpose,
  because PirBase is a curated set. Their structure signal is expected to be weak
  (short, PIWI-bound, mostly germline, size-selected libraries). They do not take
  reads from host transcripts (STAR `AllBestScore`, no primary-only, MAPQ 0).
- **Taxids.** The release 27 references already used the corrected metadata
  taxids (yeast 559292, E. coli 511145, the dengue serotype taxids, rice 39947).
  The new organisms use the corrected ones too: P. falciparum 36329 is the 3D7
  assembly GCA000002765v3, durum wheat 4567 is the Svevo assembly.
- **Script changes.** The builders retry on HTTP 5xx and on timeouts. Overlap ties
  are broken by record ID, because EBI Search result order is not stable (E. coli:
  by taxid order 511145 then 83333, then ID). `make_transcriptome_fa.py` writes
  IUPAC codes as N, because rf-count loops on them. `_PRIMARY_CHROM` accepts
  wheat-style `1A` and `Pf3D7_*_v3` chromosomes. Record sequences are read from the
  release's `sequences/rnacentral_species_specific_ids.fasta.gz` instead of
  rnacentral.org, which is slow after a release.
- **rf-eval windows re-pointed.** Human, mouse and yeast renumbered their copies;
  `scripts/check_rfeval_ids.py` moved the windows to the new IDs, all at the same loci.
- **New organisms.** Eukaryotes: arabidopsis_thaliana, danio_rerio,
  plasmodium_falciparum, triticum_turgidum. Bacteria and archaea: bacillus_subtilis,
  listeria_monocytogenes, methanosarcina_acetivorans, pseudomonas_savastanoi. Viruses:
  yellow_fever_virus, west_nile_virus, murine_norovirus, porcine_epidemic_diarrhea_virus,
  human_bocavirus, bluetongue_virus, hiv. New dengue and influenza strains were added to
  the existing folders.
- **No folder: K. lactis (rnastruct00128) and D. hansenii (rnastruct00127).**
  RNAcentral has no GFF3 for them, so the pipeline downloads the Ensembl Fungi
  genome and GTF (ASM251v1 and CBS767 ASM644v2). Ensembl Fungi keeps them in
  `fungi_ascomycota1_collection/`, which the pipeline searches from
  nf-core/rnastructurome `feature/rnacentral-flag` commit `d3cae92`; earlier
  versions fall back to an NCBI search.
- **Left out.** Large rf-eval intermediates (Rfam CM databases, cmscan output,
  RDATs, alignments) and the unfiltered GTFs, which come from
  `https://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/28.0/genome_coordinates/gff3/`.
