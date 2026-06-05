# Results Directory

This directory stores reproducible outputs.

By default, large QIIME 2 artifacts and visualizations are ignored by git:

- `results/qiime2/`
- `results/qzv/`
- `results/exported/`
- `results/phyloseq/`

Small final tables and figures can be committed when they do not contain
sensitive information.

Suggested layout:

```text
results/qiime2/     QIIME 2 artifacts
results/qzv/        QIIME 2 visualizations
results/exported/   exported BIOM, FASTA, and taxonomy files
results/tables/     final TSV/CSV summary tables
results/figures/    final plots
results/reports/    final reports
```

