# TDA Analysis

TDA means Topological Data Analysis. In this repository, the starter TDA module
uses persistent homology to summarize the shape of sample-to-sample microbiome
variation.

## What This Adds

The starter TDA workflow:

1. Reads the exported ASV feature table.
2. Matches samples to `config/metadata_template.tsv`.
3. Converts counts to relative abundance by default.
4. Computes a Bray-Curtis sample distance matrix.
5. Runs persistent homology with Ripser.
6. Saves persistence pairs, diagrams, barcode plots, and MDS coordinates.

If a Bray-Curtis distance matrix already exists, the workflow can use that
directly instead of recalculating distances from the ASV table.

## Inputs

```text
results/tables/feature-table.tsv
config/metadata_template.tsv
```

The feature table is created by `scripts/03_taxonomy.sh` after QIIME 2 export.

## Run

```bash
python scripts/05_tda_persistent_homology.py \
  --feature-table results/tables/feature-table.tsv \
  --metadata config/metadata_template.tsv \
  --grouping-column group
```

Optional parameters:

```bash
python scripts/05_tda_persistent_homology.py \
  --transform relative_abundance \
  --metric braycurtis \
  --max-dim 2 \
  --min-feature-prevalence 2
```

If you already have a Bray-Curtis matrix:

```bash
python scripts/05_tda_persistent_homology.py \
  --distance-matrix path/to/bray_curtis.tsv \
  --metadata config/metadata_template.tsv \
  --grouping-column group
```

The distance matrix should be square, symmetric, non-negative, have sample IDs
as row names and column names, and have zeroes on the diagonal.

## Outputs

```text
results/tables/tda_distance_matrix.tsv
results/tables/tda_persistence_pairs.tsv
results/tables/tda_mds_coordinates.tsv
results/tables/tda_summary.tsv
results/figures/tda_persistence_diagram.png
results/figures/tda_barcode.png
results/figures/tda_mds.png
```

## Interpretation Guide

- H0 features describe connected components merging across the filtration.
- H1 features can suggest loop-like structure in the sample distance space.
- H2 features can suggest void-like structure, but they usually require many
  samples and careful interpretation.
- Long persistence is generally more notable than short persistence.

For microbiome data, TDA should be treated as exploratory unless paired with a
clear statistical design, sensitivity checks, and biological validation.

## Recommended Sensitivity Checks

- Compare relative abundance and Hellinger transforms.
- Repeat after filtering rare ASVs.
- Repeat after removing low-depth samples.
- Check whether TDA structure is driven by sequencing run, extraction batch, or
  PCR batch.
- Compare Bray-Curtis with another distance metric when scientifically justified.

## References

- Ripser.py documentation: https://ripser.scikit-tda.org/
- Persim documentation: https://persim.scikit-tda.org/
- scikit-learn pairwise distances: https://scikit-learn.org/stable/modules/generated/sklearn.metrics.pairwise_distances.html
