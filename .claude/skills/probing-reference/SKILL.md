---
name: probing-reference
description: Build RNAcentral-based reference FASTA/GTF files (and, for new organisms, rf-eval .db + windows.tsv) for nf-core/rnastructurome. Use when asked to regenerate references for a new RNAcentral release (e.g. release 28), add an organism, filter an RNAcentral GFF3, place RNAcentral records on a bacterial/viral genome, or build rf-eval references.
---

# Building probing references

Read the repository `README.md` first (layout and method), then the
`release27/<organism>/README.md` closest to the job.

**Layout rules:**

- Code and fixed inputs live in `scripts/`; each release's data lives in
  `release<N>/<organism>/`.
- `reference/` holds FASTA/GTF only, and `rfeval/` holds `.db` + `windows.tsv` only.
- README, `MANIFEST*.tsv`, `filter.log` and the E. coli `records.fa` cache sit at
  the organism level.
- Never edit an older release's folder.

## A. New RNAcentral release, existing organisms

0. **Rfam allowlist (human, mouse):** regenerate their `RFAM_SEED_ALLOWLIST`
   blocks from the new GFF3 with `scripts/rfam_seed_allowlist.py` (README, step 1)
   before building. The list is per release; stale IDs only miss new records.
1. **Eukaryotes:** `scripts/<org>/build_reference.sh <N>` for human, mouse, yeast
   and rice, run from the repo root. Compare `filter.log` with the previous
   release's; a big jump in any step means the GFF3 changed shape (chromosome
   names, database names).
2. **Bacteria and viruses:** for each one,
   `mkdir -p release<N>/<org> && cd release<N>/<org> && python3 ../../scripts/<org>/build_<x>_refs.py`.
   Then, from the root, `python3 scripts/make_transcriptome_fa.py release<N>/<org>`.
   The builder writes genome FASTA/GTF to `work/` (git-ignored); only
   `reference/<stem>.transcriptome.fa` is kept.
3. **rf-eval:** never rebuild the structures; Rfam structures don't depend on the
   RNAcentral release. Carry each organism's `rfeval/` into `release<N>/` with
   `scripts/check_rfeval_ids.py <old windows.tsv> <old GTF> <new GTF> --outdir release<N>/<org>/rfeval`
   (bacteria/viruses: the builders' `work/<stem>.gtf`). It re-points renumbered
   `URS…_<taxid>.<n>` IDs and exits 1 on anything needing a decision
   (`candidates`, `missing`).
4. **Params:** `sed -i 's#/release<prev>/#/release<N>/#' params/*.yaml` once every
   organism's `release<N>/` has `reference/` and `rfeval/`.
5. **Run the checks** in section D, update each organism README's numbers, and
   write `release<N>/README.md`.

## B. New eukaryote (RNAcentral has a GFF3 for it)

Check `https://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/<N>.0/genome_coordinates/gff3/`.

- **Script:** copy `scripts/yeast/build_reference.sh` to
  `scripts/<org>/build_reference.sh` and change `sp=` and the folder name.
- **Flags:** keep `--pirna-max-loci 5`, the standing choice for every species.
- **`--no-prefilter`:** run once without it. If the pre-filter removes most
  records, the species' sources (e.g. SGD, Ensembl Plants) are not on
  `AUTHORITATIVE_DBS`; add the flag.
- **Chromosome names:** check them against `_PRIMARY_CHROM` (numbered, roman
  numerals, X/Y/MT/Mito/Mt/Pt). Any other naming silently drops whole
  chromosomes. Compare the GFF3's column 1 with what survives.
- **Rfam allowlist:** `RFAM_SEED_ALLOWLIST` only lists human and mouse IDs. If the
  new species uses the pre-filter, run
  `scripts/rfam_seed_allowlist.py Rfam.seed.gz <species>.gff3.gz <species>.fasta`
  and paste its output in. Make the species FASTA by extracting the taxid from
  `rnacentral_species_specific_ids.fasta.gz` (8.7 GB, about 10 minutes).
- **Size:** the build scripts gzip every GTF (`gzip -9 -n`, so reruns are
  byte-identical). The pipeline reads `.gtf.gz`.
- **IDs:** transcript IDs are `URS…_<taxid>.<n>`. Anything linking back to
  RNAcentral splits on the last `.`.

## C. New bacterium or virus (no GFF3)

Copy the closest builder into `scripts/<org>/`: `sarscov2/build_sars_refs.py` for
a virus (the most general), `ecoli/build_ecoli_refs.py` for a bacterium. Builders
write genome FASTA/GTF to `work/` under the current directory;
`make_transcriptome_fa.py` turns them into `reference/<stem>.transcriptome.fa`,
the only file kept.

- **Genome:** the NCBI accession the dataset's probing metadata names. For
  GISAID-only genomes, read the sequence from the authors' GEO ShapeMapper `.map`
  (column 4 = base). Check it against the reference strain for adapter sequence
  and trimmed ends.
- **Records:** RNAcentral records for the taxid, at least 60 nt and at least 99 %
  identical. Where records overlap on a strand, keep the best by identity, then
  length. A record one base off turns true reads into MaP mutations.
- **Genome FASTA** (`work/`): the genome only. Records as extra contigs make every
  read map twice.
- **GTF** (`work/`): holds the records on genome coordinates. For viruses the
  genome is also a transcript, plus `_5UTR`/`_3UTR` cut at the CDS where they
  are long enough. A bacterial chromosome is **never** a GTF transcript (rf-fold
  on 4.6 Mb never finishes).
- **Routes:** these run with `--transcriptome --fasta <stem>.transcriptome.fa` and no
  `--gtf`.
  Don't mix routes: `--transcriptome` with the genome-only `.fa` gave E. coli
  rf-count "Covered 1".
- **rnacentral.org:** it throttles after about 60 quick calls by serving its HTML
  home page. Fetch at 1 per second with backoff, send `Accept-Encoding: gzip`,
  and cache results so a run can resume.

## D. Checks before handing over

- **Only ACGTN.** `tr -d 'ACGTNacgtn\n' < x.fa | wc -c` must be 0. One IUPAC base
  (Y, R, …) makes rf-count loop forever at 100 % CPU, and GISAID genomes often
  carry them.
- **No parentheses in transcript IDs.** RNA Framework's XML parser hangs on them.
- **IDs match.** Every GTF `transcript_id` has its sequence in the FASTA, and
  every `windows.tsv` `ref_seq_id` is a transcript the run will use.
- **Rebuilding a release that already exists:** outputs must match
  byte-for-byte; otherwise explain the difference.

## E. rf-eval for a new organism

The archived builders in `scripts/` show the method: the virus
`build_rfeval_refs.py` and `build_rrna_common.py`. They ran with their inputs
beside them, so adapt the paths, and write outputs to `release<N>/<org>/rfeval/`.
The aim is a small, certain yardstick.

- **GA only.** Use `cmsearch --cut_ga` / `cmscan --cut_ga`. The default `!` flag
  is E ≤ 0.01, not GA. Filter in code; never edit a `.db` by hand.
- **Structure.** Align each hit with `cmalign --mapali <seed> --mapstr` and project
  the family `SS_cons` (`project_ss_cons.py`). Drop pseudoknots unless the family
  *is* the pseudoknot (influenza RF01099; the run then needs
  `rfeval_keep_pseudoknots = true`). Keep only structures with at least
  0.10 base pairs per nucleotide.
- **Common rRNA:** one locus per molecule, inside the SSU's rDNA unit (10 kb).
  Without that rule the search picks pseudogenes and mitochondrial rRNA. If the
  assembly lacks the rDNA arrays (mouse GRCm39), say so in the README rather than
  adding a non-standard contig.
- **Identical copies:** after the first run, pick the copy with the most measured
  positions in its RDATs (as was done for human, from runs rnastruct000038/40).
- **Windows keyed on the run's transcripts.** E. coli hits found on the
  chromosome had to be moved onto URS records; genome-keyed windows gave "no
  windows written".

## F. Finish

- **README per organism:** sources and accessions, filter settings and the
  per-step counts from `filter.log`, what was excluded and why, and the pipeline
  options to run it.
- **Ask before committing.** Commits are a subject line only.
