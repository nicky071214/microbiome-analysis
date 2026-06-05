# Standard Operating Procedure

## Scope

This SOP describes the starter workflow for PacBio CCS/HiFi full-length 16S
rRNA microbiome data.

## 1. Before Analysis

Confirm the following with the sequencing provider or wet-lab record:

- Marker: full-length 16S rRNA.
- Platform: PacBio CCS/HiFi.
- Sample demultiplexing status.
- Exact primer sequences, written 5' to 3'.
- Barcode/sample mapping.
- Expected amplicon size range.
- Sequencing run ID and delivery date.

## 2. Data Handling

Raw FASTQ/BAM files must not be committed to GitHub.

Place raw files in controlled storage and record:

- storage location
- file name
- file size
- checksum
- sequencing run
- sample ID

Use anonymized sample IDs only.

## 3. Metadata Preparation

Start from `config/metadata_template.tsv`.

Check that:

- `sample-id` exactly matches the manifest.
- no identifying human information is present.
- batch columns are filled when relevant.
- group labels are consistent.
- missing values are represented consistently.

Run:

```bash
python scripts/check_metadata.py config/metadata_template.tsv
```

## 4. Manifest Preparation

Start from `config/pacbio_manifest.tsv`.

Use absolute FASTQ paths and `forward` direction for demultiplexed single-end
PacBio CCS reads.

## 5. QIIME 2 Import

Run:

```bash
bash scripts/01_import_qiime2.sh
```

Review:

```text
results/qzv/demux-pacbio-ccs.qzv
results/qzv/metadata.qzv
```

## 6. Denoising

Use DADA2 `denoise-ccs`, which is designed for single-end PacBio CCS reads.

Run after replacing primers:

```bash
FRONT_PRIMER="YOUR_FORWARD_PRIMER" \
ADAPTER="YOUR_REVERSE_PRIMER_OR_ADAPTER" \
bash scripts/02_denoise_dada2_ccs.sh
```

Review:

```text
results/qzv/feature-table.qzv
results/qzv/rep-seqs.qzv
results/qzv/denoising-stats.qzv
```

Tune `MIN_LEN`, `MAX_LEN`, and other filters if many reads are lost.

## 7. Taxonomy

Use a classifier trained for full-length 16S.

Run:

```bash
CLASSIFIER="references/classifier-full-length-16s.qza" bash scripts/03_taxonomy.sh
```

Review:

```text
results/qzv/taxonomy.qzv
results/qzv/taxa-bar-plots.qzv
```

## 8. Diversity and Summary Figures

After exporting feature table and taxonomy:

```bash
Rscript scripts/04_diversity_phyloseq.R
```

Outputs:

```text
results/tables/alpha_diversity.tsv
results/tables/bray_curtis_pcoa.tsv
results/figures/alpha_diversity_shannon.png
results/figures/bray_curtis_pcoa.png
results/phyloseq/pacbio_full_length_16s_phyloseq.rds
```

## 9. Reporting

Record:

- QIIME 2 version.
- DADA2 parameters.
- primer sequences.
- classifier database and version.
- sample exclusion criteria.
- read counts before and after denoising.
- statistical tests used.

