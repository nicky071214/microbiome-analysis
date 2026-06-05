#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(vegan)
})

feature_table_path <- Sys.getenv("FEATURE_TABLE", "results/tables/feature-table.tsv")
taxonomy_path <- Sys.getenv("TAXONOMY", "results/tables/taxonomy.tsv")
metadata_path <- Sys.getenv("METADATA", "config/metadata_template.tsv")
figures_dir <- Sys.getenv("FIGURES_DIR", "results/figures")
tables_dir <- Sys.getenv("TABLES_DIR", "results/tables")
phyloseq_dir <- Sys.getenv("PHYLOSEQ_DIR", "results/phyloseq")
grouping_column <- Sys.getenv("GROUPING_COLUMN", "group")

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(phyloseq_dir, recursive = TRUE, showWarnings = FALSE)

stop_if_missing <- function(path) {
  if (!file.exists(path)) {
    stop("Missing required file: ", path, call. = FALSE)
  }
}

stop_if_missing(feature_table_path)
stop_if_missing(taxonomy_path)
stop_if_missing(metadata_path)

read_feature_table <- function(path) {
  first_line <- readLines(path, n = 1)
  skip_lines <- if (grepl("^# Constructed from biom file", first_line)) 1 else 0
  feature_table <- read.delim(
    path,
    skip = skip_lines,
    check.names = FALSE,
    comment.char = "",
    stringsAsFactors = FALSE
  )
  names(feature_table)[1] <- "FeatureID"
  rownames(feature_table) <- feature_table$FeatureID
  feature_table$FeatureID <- NULL
  as.matrix(feature_table)
}

feature_matrix <- read_feature_table(feature_table_path)
storage.mode(feature_matrix) <- "numeric"

metadata <- read.delim(
  metadata_path,
  check.names = FALSE,
  comment.char = "",
  stringsAsFactors = FALSE
)

if (!"sample-id" %in% names(metadata)) {
  stop("Metadata must contain a sample-id column.", call. = FALSE)
}

rownames(metadata) <- metadata[["sample-id"]]

common_samples <- intersect(colnames(feature_matrix), rownames(metadata))
if (length(common_samples) == 0) {
  stop("No overlapping sample IDs between feature table and metadata.", call. = FALSE)
}

feature_matrix <- feature_matrix[, common_samples, drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]

taxonomy <- read.delim(
  taxonomy_path,
  check.names = FALSE,
  comment.char = "",
  stringsAsFactors = FALSE
)

alpha <- data.frame(
  sample_id = colnames(feature_matrix),
  observed_features = colSums(feature_matrix > 0),
  shannon = vegan::diversity(t(feature_matrix), index = "shannon"),
  total_reads = colSums(feature_matrix),
  stringsAsFactors = FALSE
)

alpha <- merge(alpha, metadata, by.x = "sample_id", by.y = "sample-id", all.x = TRUE)
write.table(
  alpha,
  file = file.path(tables_dir, "alpha_diversity.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

plot_group <- if (grouping_column %in% names(alpha)) grouping_column else "sample_id"

p_alpha <- ggplot(alpha, aes(x = .data[[plot_group]], y = shannon, fill = .data[[plot_group]])) +
  geom_boxplot(outlier.shape = NA, alpha = 0.55) +
  geom_jitter(width = 0.12, size = 2) +
  labs(x = plot_group, y = "Shannon diversity") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 35, hjust = 1))

ggsave(
  filename = file.path(figures_dir, "alpha_diversity_shannon.png"),
  plot = p_alpha,
  width = 7,
  height = 4.5,
  dpi = 300
)

if (ncol(feature_matrix) >= 3) {
  bray <- vegan::vegdist(t(feature_matrix), method = "bray")
  pcoa <- cmdscale(bray, eig = TRUE, k = 2)
  pcoa_df <- data.frame(
    sample_id = rownames(pcoa$points),
    PC1 = pcoa$points[, 1],
    PC2 = pcoa$points[, 2],
    stringsAsFactors = FALSE
  )
  pcoa_df <- merge(pcoa_df, metadata, by.x = "sample_id", by.y = "sample-id", all.x = TRUE)
  write.table(
    pcoa_df,
    file = file.path(tables_dir, "bray_curtis_pcoa.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  color_column <- if (grouping_column %in% names(pcoa_df)) grouping_column else "sample_id"
  p_pcoa <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = .data[[color_column]])) +
    geom_point(size = 3) +
    labs(x = "PCoA 1", y = "PCoA 2", color = color_column) +
    theme_bw(base_size = 12)

  ggsave(
    filename = file.path(figures_dir, "bray_curtis_pcoa.png"),
    plot = p_pcoa,
    width = 6,
    height = 5,
    dpi = 300
  )
}

if (requireNamespace("phyloseq", quietly = TRUE)) {
  tax_matrix <- NULL
  if (all(c("Feature ID", "Taxon") %in% names(taxonomy))) {
    tax_strings <- strsplit(taxonomy[["Taxon"]], ";")
    ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    tax_matrix <- matrix(NA_character_, nrow = nrow(taxonomy), ncol = length(ranks))
    for (i in seq_along(tax_strings)) {
      values <- trimws(tax_strings[[i]])
      values <- sub("^[A-Za-z]__", "", values)
      values <- sub("^D_[0-9]+__", "", values)
      tax_matrix[i, seq_len(min(length(values), length(ranks)))] <- values[seq_len(min(length(values), length(ranks)))]
    }
    colnames(tax_matrix) <- ranks
    rownames(tax_matrix) <- taxonomy[["Feature ID"]]
  }

  ps <- phyloseq::phyloseq(
    phyloseq::otu_table(feature_matrix, taxa_are_rows = TRUE),
    phyloseq::sample_data(metadata)
  )

  if (!is.null(tax_matrix)) {
    ps <- phyloseq::merge_phyloseq(ps, phyloseq::tax_table(tax_matrix))
  }

  saveRDS(ps, file.path(phyloseq_dir, "pacbio_full_length_16s_phyloseq.rds"))
}

message("Wrote alpha diversity table and figures.")
