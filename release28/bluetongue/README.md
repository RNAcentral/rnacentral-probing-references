# Bluetongue virus serotype 1, South Africa (rnastruct00077) — reference

## Reference

`BTV1.transcriptome.fa`: the ten segment records, one transcript per segment,
FJ969719.1–FJ969728.1 (`build_btv_refs.py`, `MANIFEST.tsv`). The dsRNA segments were
probed as T7 transcripts, so each is a whole record with no UTR split, as for
influenza and rotavirus. RNAcentral has no full-length BTV transcript, so there is no
second record and no `_rnac` comparison.

| segment | protein | accession | nt |
|---|---|---|---|
| 1 | VP1 | FJ969719.1 | 3944 |
| 2 | VP2 | FJ969720.1 | 2940 |
| 3 | VP3 | FJ969721.1 | 2772 |
| 4 | VP4 | FJ969722.1 | 1981 |
| 5 | NS1 | FJ969724.1 | 1772 |
| 6 | VP5 | FJ969723.1 | 1635 |
| 7 | VP7 | FJ969725.1 | 1156 |
| 8 | NS2 | FJ969726.1 | 1125 |
| 9 | VP6 | FJ969727.1 | 1049 |
| 10 | NS3 | FJ969728.1 | 822 |

The GenBank titles of FJ969723 and FJ969724 swap segments 5 and 6. The segments are
mapped by gene: FJ969724.1 is NS1 (segment 5) and FJ969723.1 is VP5 (segment 6).

## rf-eval reference

None.

Rebuild the reference: `cd release<N>/bluetongue && python3 ../../scripts/bluetongue/build_btv_refs.py` (Biopython, network) writes the genome FASTA/GTF to `work/`; then, from the root,
`python3 scripts/make_transcriptome_fa.py release<N>/bluetongue` cuts `reference/<stem>.transcriptome.fa` from them.

## Files

- `reference/BTV1.transcriptome.fa` — the ten segments; `--fasta` with `--transcriptome`.
  `MANIFEST.tsv` lists them.

## Running

```
nextflow run . -profile <…> -resume \
  --transcriptome --fasta release28/bluetongue/reference/BTV1.transcriptome.fa
```
