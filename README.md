# Microbiome Analysis

PacBio full-length 16S rRNA microbiome analysis scaffold.

This repository tracks analysis code, metadata templates, configuration,
documentation, and reproducible outputs for PacBio circular consensus sequencing
(CCS/HiFi) full-length 16S amplicon data.

## Data Safety

Do not commit sensitive or large data to GitHub.

- Do not commit raw FASTQ, BAM, or vendor delivery archives.
- Do not commit patient identifiers, names, birthdays, medical record numbers,
  addresses, or contact details.
- Use anonymized sample IDs in metadata.
- Store raw sequencing data in controlled storage, such as institutional
  storage, cloud storage, NCBI SRA, Zenodo, or another approved repository.
- Track file checksums so external data can be verified later.

## Intended Data Type

- Platform: PacBio CCS/HiFi
- Marker: full-length 16S rRNA gene
- Typical target size: about 1,500 bp
- Default denoising approach: QIIME 2 DADA2 `denoise-ccs`
- Downstream analysis: QIIME 2 plus R/phyloseq-style summaries
- Optional exploratory analysis: TDA persistent homology on microbiome distance
  space

## Repository Layout

```text
config/       Project parameters, sample metadata, and manifest template
data/         Data notes only; raw sequencing files are ignored by git
docs/         SOP, analysis plan, and project documentation
references/   Taxonomy classifier notes and external database references
results/      Reproducible outputs, figures, tables, and reports
scripts/      Analysis and utility scripts
workflow/     Optional workflow manager files, such as Snakemake
```

## Quick Start

1. Fill in `config/metadata_template.tsv` with anonymized sample information.
2. Fill in `config/pacbio_manifest.tsv` with absolute paths to demultiplexed
   PacBio CCS FASTQ files.
3. Install the QIIME 2 amplicon distribution by following the official QIIME 2
   installation instructions.
4. Activate your QIIME 2 environment.
5. Replace the primer placeholders in `scripts/02_denoise_dada2_ccs.sh` or pass
   them as environment variables.
6. Run the starter scripts in order:

```bash
bash scripts/01_import_qiime2.sh
FRONT_PRIMER="YOUR_FORWARD_PRIMER" ADAPTER="YOUR_REVERSE_PRIMER" bash scripts/02_denoise_dada2_ccs.sh
CLASSIFIER="references/classifier-full-length-16s.qza" bash scripts/03_taxonomy.sh
Rscript scripts/04_diversity_phyloseq.R
python scripts/05_tda_persistent_homology.py
```

If you already have a Bray-Curtis distance matrix, pass it directly:

```bash
python scripts/05_tda_persistent_homology.py --distance-matrix path/to/bray_curtis.tsv
```

## Standard Workflow

1. Receive PacBio CCS/HiFi reads and the sequencing report.
2. Confirm demultiplexing, primer strategy, and barcode/sample mapping.
3. Prepare metadata and manifest files.
4. Import reads into QIIME 2.
5. Inspect demultiplexed sequence quality.
6. Denoise PacBio CCS reads with DADA2 `denoise-ccs`.
7. Assign taxonomy using a full-length 16S-compatible reference database.
8. Generate feature table, representative sequences, taxonomy table, diversity
   metrics, and visualizations.
9. Optionally run TDA persistent homology on the exported feature table.
10. Export final tables and figures for reporting.

## Important Notes

- Replace the default primer placeholders with the exact primers used by the
  sequencing provider.
- Tune DADA2 length filters after reviewing the QIIME 2 demultiplexed sequence
  summary.
- Use a taxonomy classifier appropriate for full-length 16S data, not a
  short-region V3-V4-only classifier.
- Commit scripts, metadata templates, documentation, and small summary outputs;
  keep raw data external.

## Useful References

- QIIME 2 DADA2 `denoise-ccs`: https://docs.qiime2.org/2024.10/plugins/available/dada2/denoise-ccs/
- QIIME 2 amplicon documentation: https://amplicon-docs.qiime2.org/
- DADA2: https://benjjneb.github.io/dada2/
- Ripser.py: https://ripser.scikit-tda.org/
- Persim: https://persim.scikit-tda.org/
