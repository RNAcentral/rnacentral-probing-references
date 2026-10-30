# rnacentral-probing-reference

Reference files for RNA chemical-probing runs with
[nf-core/rnastructurome](https://github.com/nf-core/rnastructurome), built from
RNAcentral: the FASTA/GTF to map and fold against, and the known structures that
rf-eval scores folds with.

```
scripts/                     reused for every release
  filter_rnacentral_gff3.py  GFF3 filter (eukaryotes)
  gff3_to_gtf_exons.py       GFF3 → GTF
  make_transcriptome_fa.py   <stem>.transcriptome.fa for bacteria and viruses
  build_small_rna_transcriptome.py  small-RNA transcriptome for size-selected libraries
  rfam_seed_allowlist.py     regenerates the filter's Rfam allowlist
  check_rfeval_ids.py        checks rf-eval windows against a new release's GTF
  rfeval_*.py                rf-eval builders for new organisms (rRNA, viral elements, small RNAs)
  rrna_cm/                   the Rfam models the rf-eval builders search with
  <organism>/                that organism's builder(s), plus the archived rf-eval builders
params/                      nf-core/rnastructurome params, one <dataset_id>-params.yaml per dataset
release27/
  README.md                  what this release holds
  <organism>/
    README.md, MANIFEST*.tsv, filter.log
    reference/               FASTA / GTF only (what the pipeline is given)
    work/                    build intermediates, git-ignored
    rfeval/                  rf-eval .db + windows.tsv only
```

File names carry no release number; only the `release<N>/` folder does, so moving
to a new release changes one path segment in the params files.

From release 28, organism folders (`release<N>/<organism>/` and `scripts/<organism>/`)
are named after the metadata `scientific_name`, lowercased with other characters
replaced by `_` (`oryza_sativa`, `sars_cov_2`, `hiv`). That is also how the pipeline
turns a samplesheet organism into a key. Release 27 keeps its common names (`rice`,
`sarscov2`).

## Building a new release

Run everything from the repository root; outputs go to `release<N>/<organism>/`.
The steps below use release 28.

### 1. Refresh the Rfam allowlist (human, mouse)

The filter keeps records that fail its database check if their sequence is in an
Rfam seed alignment. Those IDs are hardcoded in `RFAM_SEED_ALLOWLIST` in
`scripts/filter_rnacentral_gff3.py`, one block per species, and were made from
release 27. New records that match a seed are dropped until the list is
regenerated. Yeast and rice skip the database check (`--no-prefilter`), so only
human and mouse need this. Download the inputs into an empty folder:

```
https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.seed.gz
https://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/28.0/genome_coordinates/gff3/homo_sapiens.GRCh38.gff3.gz
https://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/28.0/sequences/rnacentral_species_specific_ids.fasta.gz
```

then, for each species:

```
python3 scripts/rfam_seed_allowlist.py Rfam.seed.gz homo_sapiens.GRCh38.gff3.gz rnacentral_species_specific_ids.fasta.gz
```

It prints the IDs grouped by RNA type; replace that species' block in
`RFAM_SEED_ALLOWLIST` with them. The species FASTA covers every species
(8.7 GB), but it is read as a stream and is the same file for human and mouse.

### 2. Build the references

**Human, mouse, yeast, rice.** RNAcentral publishes genome coordinates for
these. Each script downloads that release's GFF3 and filters it:

```
scripts/homo_sapiens/build_reference.sh 28
```

The FASTA is the Ensembl genome, which the pipeline downloads itself.

**Bacteria and viruses** (E. coli, dengue, influenza, rotavirus, SARS-CoV-2,
Zika). RNAcentral has no coordinates for these, so the builder fetches the
pinned NCBI genome and places the current RNAcentral records on it:

```
mkdir -p release28/dengue_virus && cd release28/dengue_virus
python3 ../../scripts/dengue_virus/build_dengue_refs.py
cd ../.. && python3 scripts/make_transcriptome_fa.py release28/dengue_virus
```

The builder writes the genome FASTA/GTF to `work/`, and `make_transcriptome_fa.py`
cuts `reference/<stem>.transcriptome.fa` from them. That one file is the reference:
run it with `--transcriptome --fasta <stem>.transcriptome.fa` and no `--gtf`.

**Small-RNA libraries** (size-selected, metadata `rna_type: sRNA`). On the genome
route their tRNA reads are counted on whatever long records span the tRNA genes, and
rf-eval falls to chance; mapped to the mature molecules instead, it does not (see
`release28/homo_sapiens/README.md`). Build the small-RNA transcriptome after the
organism's filtered GTF, from the repository root:

```
python3 scripts/build_small_rna_transcriptome.py homo_sapiens 28 [rnacentral_species_specific_ids.fasta.gz]
python3 scripts/rfeval_small_rna.py homo_sapiens 28
```

The first takes the filtered GTF's tRNA, snRNA, scaRNA, snoRNA, pre-miRNA, Y RNA,
vault, SRP, RNase P/MRP and telomerase records, plus the rRNA records listed for the
organism (the rDNA units rf-eval uses, 5S, mitochondrial rRNA). At each locus it keeps
the best record by the filter's own rule (overlapping, lengths within 2×, best
database rank, then most databases, then length). It then keeps one record per URS,
drops records over 500 nt, and merges records within 2 substitutions of a
better-ranked one, whose reads would otherwise land on the wrong copy and count as
mutations. Sequences come from the release's species FASTA (streamed from the FTP if
no local copy is given; `work/small_rna.seqs.fa` caches them). It writes
`reference/<species>.small_rna.transcriptome.fa` and `MANIFEST-small_rna.tsv`, which
lists each record's class, genomic copies and merged records. The second picks one
representative structure per family (best GA hit) for `rfeval/small_rna.*`. Run with
`transcriptome: true`, `fasta` and `bowtie_all: false`.

### 3. Carry the rf-eval files over

rf-eval structures are built from Rfam, not RNAcentral, so they are not rebuilt.
Their builders are archived in `scripts/` exactly as they ran, and their paths
need adjusting before any rerun. But a windows file names the transcripts it
scores, and a new release can renumber or drop them (human, mouse and yeast use
IDs such as `URS…_9606.1467066`; E. coli, dengue and SARS-CoV-2 also window on
RNAcentral records). For each organism, check the windows against the new GTF and
write the files into the new release:

```
python3 scripts/check_rfeval_ids.py release27/human/rfeval/rrna_common.windows.tsv \
    release27/human/reference/homo_sapiens.GRCh38.filtered.gtf.gz \
    release28/homo_sapiens/reference/homo_sapiens.GRCh38.filtered.gtf.gz \
    --outdir release28/homo_sapiens/rfeval
```

For bacteria and viruses, pass the builder's `work/<stem>.gtf` files (from both
releases) as the GTFs. Each transcript is reported as:

- `ok`: same ID at the same place.
- `moved`: same ID, now another copy of the same sequence.
- `repoint`: the ID is gone, but the same URS is at the same locus under a new
  ID. The written files use the new ID.
- `candidates`: the same URS is only at other loci. Pick one by hand.
- `missing`: the sequence is no longer in the GTF, so the structure has to be rebuilt.

The script exits 1 if anything is left for you to resolve. A URS fixes the
sequence, so window coordinates stay valid under any ID that shares it.

### 4. Point the params files at the new release

```
sed -i 's#/release27/#/release28/#' params/*.yaml
```

This moves both the reference and rf-eval paths, so steps 2 and 3 must have
filled `release28/` for every organism first. Then commit everything.

## Running nf-core/rnastructurome

The run procedure is in section 3 of the [RNAcentral chemical probing
SOP](https://embl.atlassian.net/wiki/spaces/RNAC/pages/229146695/RNAcentral+chemical+probing#3.-Running-nf-core/rnastructurome).
The repo is checked out on codon at
`/hps/nobackup/agb/rnacentral/chemicalprob/rnacentral-probing-references`, and
the params files use absolute paths into that checkout.

Each dataset in
[rnacentral-probing-metadata](https://github.com/RNAcentral/rnacentral-probing-metadata)
gets its own `params/<dataset_id>-params.yaml`, with `rnacentral: true`,
`structextract: true` and the reference options that fit the organism (see
[rnastruct00140-params.yaml](params/rnastruct00140-params.yaml)):

| organism | options |
|---|---|
| eukaryotes with a filtered GTF here | `gtf` (genome route; the pipeline downloads the Ensembl genome) |
| bacteria, archaea, and virus strains with a FASTA here | `transcriptome: true`, `fasta` |
| small-RNA libraries (`rna_type: sRNA`) | `transcriptome: true`, `fasta` (`<species>.small_rna.transcriptome.fa`), `bowtie_all: false`, `rnacentral_min_mapped_pct: 20` (reads from longer RNAs have nowhere to map, so 23–35% map) |
| mRNA libraries (poly(A)-selected or targeted mRNAs) | none (Ensembl route): the RNAcentral GTF holds only ncRNAs, so it leaves almost nothing to count |
| anything else | none; the pipeline downloads the Ensembl/NCBI reference for the samplesheet's `organism` |

Add `rfeval_reference` and `rfeval_windows` if the organism has an `rfeval/`
folder, and put any dataset-specific options (adapters, trimming,
normalisation) in the same file. Datasets whose metadata `comment` starts with
`failed QC` or `skip` are not run and need no params file.

The second line of each params file gives its status. `# status: final` means
the dataset has been run with it and checked. `# status: draft` means it was
generated from the metadata and has not been run yet. Change it to `final` once
a run has been checked.

Run the dataset on codon:

```
nextflow run main.nf \
  -profile codon,singularity \
  --input  "${BASE}/FASTQ/samplesheet/${DATASET_ID}_samplesheet.csv" \
  --outdir "${BASE}/RESULTS/${DATASET_ID}" \
  -params-file "${BASE}/rnacentral-probing-references/params/${DATASET_ID}-params.yaml" \
  -resume
```

If `params/<dataset_id>-config.config` exists (process resources, which a params file
cannot set), add `-c "${BASE}/rnacentral-probing-references/params/${DATASET_ID}-config.config"`.

## How the eukaryote GTF is filtered

`scripts/filter_rnacentral_gff3.py`, in order:

0. **Drop circRNAs.** Mapping and folding assume linear RNA.
1. **piRNA cap** (`--pirna-max-loci 5`, all species). Drop every entry of a piRNA
   sequence with more than 5 entries. These are repeat-derived and their reads
   multimap, so nothing can be measured per locus.
2. **Pre-filter.** Keep main chromosomes only (numbered, roman numerals,
   X/Y/MT/Mito/Mt/Pt). Each record needs a trusted database for its RNA type
   (e.g. miRBase for miRNA, GtRNAdb for tRNA), unless its sequence is in an Rfam
   seed alignment. `--no-prefilter` skips the database check for species whose
   sources are not on the lists (yeast, rice).
3. **Same-locus duplicates.** One sequence listed once per database at the same
   place: keep the best-ranked copy. Copies at other loci all stay, with
   distinct IDs (`URS…_<taxid>.<n>`).
4. **Overlap thinning.** Overlapping records of one RNA type on one strand, with
   lengths within 2× of each other, keep one record per RNAcentral predicted gene
   (or one for the group). Records backed by more than 5 databases always stay.
   The winner is chosen by database rank, then most databases, then length.

## How bacterial and viral references are built

- **Genome:** the NCBI accession the dataset names. GISAID-only SARS-CoV-2
  variants come from the authors' GEO ShapeMapper `.map` files.
- **Records:** RNAcentral records of at least 60 nt that match the genome at
  99 % identity or better. Where records overlap, the best by identity, then
  length, is kept.
- **Build files** (in `work/`): `<stem>.fa` is the genome only; `<stem>.gtf` holds
  the records on genome coordinates. Viruses also get the genome and its
  `_5UTR`/`_3UTR` as transcripts. A bacterial chromosome is never a GTF transcript.
- **Reference:** `<stem>.transcriptome.fa`, every GTF transcript cut from the genome
  (for viruses the genome keeps its accession, so rf-eval windows match).

## Requirements

- Python 3 with Biopython, `curl` and network access (NCBI, Ensembl,
  rnacentral.org, the EBI FTP).
- rnacentral.org throttles after about 60 quick API calls, so the builders fetch
  at 1 per second (E. coli takes about 20 minutes; `records.fa` caches it).
- The archived rf-eval builders also need Docker with `rnacentral/r2dt:latest`
  (for Infernal) and `gffread`.
