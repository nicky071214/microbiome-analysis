# Analysis Plan

## Objective

Analyze PacBio full-length 16S rRNA microbiome sequencing data to summarize
community composition, diversity, and group-level patterns.

## Primary Inputs

- Demultiplexed PacBio CCS/HiFi FASTQ files.
- Sample metadata TSV.
- Full-length 16S-compatible taxonomy classifier.

## Primary Outputs

- ASV feature table.
- Representative ASV sequences.
- DADA2 denoising statistics.
- Taxonomy table.
- Taxa bar plots.
- Alpha diversity table and plot.
- Bray-Curtis PCoA table and plot.
- phyloseq R object.
- TDA persistence diagram, barcode, distance matrix, and MDS coordinates.

## Initial QC Questions

- Do all manifest sample IDs exist in metadata?
- Are any samples missing group or batch labels?
- Are read depths adequate after denoising?
- Are read lengths consistent with full-length 16S?
- Are controls included and labeled?
- Are extraction or PCR batches confounded with biological groups?

## Denoising Plan

Use QIIME 2 DADA2 `denoise-ccs`.

Starter parameters:

- `MIN_LEN=1000`
- `MAX_LEN=1800`
- `TRUNC_LEN=0`
- `MAX_EE=2.0`
- `MAX_MISMATCH=2`
- `CHIMERA_METHOD=consensus`
- `POOLING_METHOD=independent`

These should be revisited after reviewing demultiplexed quality summaries and
denoising statistics.

## Taxonomy Plan

Use a classifier trained for full-length 16S reads.

Potential reference options:

- SILVA full-length 16S classifier.
- GTDB-derived 16S classifier, if appropriate for the research question.
- Custom classifier trained with the exact primer region when needed.

Do not use a classifier trained only for V3-V4 or another short amplicon region
unless the data were generated for that region.

## Diversity Plan

Alpha diversity:

- observed ASVs
- Shannon diversity

Beta diversity:

- Bray-Curtis distance
- PCoA visualization

Potential group comparisons:

- Wilcoxon rank-sum test for two groups.
- Kruskal-Wallis test for more than two groups.
- PERMANOVA for beta diversity group effects, after checking design and batch.

## Batch and Confounding Checks

Check whether extraction batch, PCR batch, sequencing run, or collection date is
confounded with biological group. If confounding exists, interpret group
differences cautiously and document the limitation.

## Reproducibility

Each analysis run should record:

- commit SHA
- QIIME 2 version
- R version
- classifier version
- DADA2 parameters
- metadata file version
- raw data checksum file

## Optional TDA Plan

Use `scripts/05_tda_persistent_homology.py` after feature table export.

Starter settings:

- transform: relative abundance
- distance: Bray-Curtis
- maximum homology dimension: 2
- grouping column: `group`

Primary outputs:

- `results/tables/tda_persistence_pairs.tsv`
- `results/tables/tda_distance_matrix.tsv`
- `results/figures/tda_persistence_diagram.png`
- `results/figures/tda_barcode.png`
- `results/figures/tda_mds.png`

TDA results should be interpreted as exploratory and checked against batch,
sequencing run, and sample depth.
