#!/usr/bin/env python3
"""Run starter TDA on a microbiome feature table or distance matrix.

The script can start from either a BIOM-exported feature table TSV with
features as rows and samples as columns, or a precomputed square sample distance
matrix. For feature tables, it converts counts to relative abundance by default,
builds a Bray-Curtis distance matrix, computes persistent homology with ripser,
and writes publication-friendly starter plots and tables.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from persim import plot_diagrams
from ripser import ripser
from sklearn.manifold import MDS
from sklearn.metrics import pairwise_distances

DEFAULT_METADATA = "config/metadata_template.tsv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run TDA persistent homology on microbiome data.")
    parser.add_argument("--feature-table", default="results/tables/feature-table.tsv")
    parser.add_argument("--distance-matrix", default=None, help="Optional precomputed square sample distance matrix.")
    parser.add_argument("--metadata", default=DEFAULT_METADATA)
    parser.add_argument("--grouping-column", default="group")
    parser.add_argument("--out-tables", default="results/tables")
    parser.add_argument("--out-figures", default="results/figures")
    parser.add_argument("--transform", choices=["relative_abundance", "hellinger", "none"], default="relative_abundance")
    parser.add_argument("--metric", default="braycurtis")
    parser.add_argument("--max-dim", type=int, default=2)
    parser.add_argument("--threshold", type=float, default=None)
    parser.add_argument("--min-feature-total-count", type=float, default=0.0)
    parser.add_argument("--min-feature-prevalence", type=int, default=1)
    parser.add_argument("--random-state", type=int, default=42)
    return parser.parse_args()


def read_feature_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing feature table: {path}")

    with path.open(encoding="utf-8") as handle:
        first_line = handle.readline()

    skiprows = 1 if first_line.startswith("# Constructed from biom file") else 0
    table = pd.read_csv(path, sep="\t", skiprows=skiprows, comment=None)
    if table.empty:
        raise ValueError(f"Feature table has no rows: {path}")

    first_col = table.columns[0]
    table = table.rename(columns={first_col: "FeatureID"}).set_index("FeatureID")
    table = table.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    return table


def read_metadata(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing metadata: {path}")

    metadata = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    if "sample-id" not in metadata.columns:
        raise ValueError("Metadata must contain a sample-id column.")
    if metadata["sample-id"].duplicated().any():
        duplicates = sorted(metadata.loc[metadata["sample-id"].duplicated(), "sample-id"].unique())
        raise ValueError(f"Duplicate sample IDs in metadata: {', '.join(duplicates)}")
    return metadata.set_index("sample-id")


def read_distance_matrix(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing distance matrix: {path}")

    sep = "," if path.suffix.lower() == ".csv" else "\t"
    matrix = pd.read_csv(path, sep=sep, dtype=str)
    if matrix.empty:
        raise ValueError(f"Distance matrix has no rows: {path}")

    first_col = matrix.columns[0]
    matrix = matrix.rename(columns={first_col: "sample_id"}).set_index("sample_id")
    matrix.index = matrix.index.astype(str)
    matrix.columns = matrix.columns.astype(str)
    matrix = matrix.apply(pd.to_numeric, errors="coerce")

    if matrix.isna().any().any():
        raise ValueError("Distance matrix contains non-numeric or missing values.")
    if matrix.shape[0] != matrix.shape[1]:
        raise ValueError("Distance matrix must be square.")
    if set(matrix.index) != set(matrix.columns):
        raise ValueError("Distance matrix row and column sample IDs must match.")

    matrix = matrix.loc[matrix.index, matrix.index]
    if not np.allclose(matrix.values, matrix.values.T, atol=1e-8):
        raise ValueError("Distance matrix must be symmetric.")
    if not np.allclose(np.diag(matrix.values), 0, atol=1e-8):
        raise ValueError("Distance matrix diagonal must be zero.")
    if (matrix.values < 0).any():
        raise ValueError("Distance matrix cannot contain negative values.")

    return matrix


def metadata_for_samples(metadata_path: Path, sample_ids: list[str], grouping_column: str, allow_sample_only: bool) -> pd.DataFrame:
    if metadata_path.exists():
        metadata = read_metadata(metadata_path)
        common_samples = [sample for sample in sample_ids if sample in metadata.index]
        if len(common_samples) >= 3:
            return metadata.loc[common_samples]
        if not allow_sample_only:
            raise ValueError("Fewer than three samples overlap between metadata and the distance/feature table.")

    metadata = pd.DataFrame(index=sample_ids)
    metadata.index.name = "sample-id"
    metadata[grouping_column] = "all_samples"
    return metadata


def filter_features(table: pd.DataFrame, min_total_count: float, min_prevalence: int) -> pd.DataFrame:
    keep_total = table.sum(axis=1) >= min_total_count
    keep_prevalence = (table > 0).sum(axis=1) >= min_prevalence
    filtered = table.loc[keep_total & keep_prevalence]
    if filtered.empty:
        raise ValueError("No features remain after filtering. Relax feature filters.")
    return filtered


def transform_table(sample_by_feature: pd.DataFrame, transform: str) -> pd.DataFrame:
    if transform == "none":
        return sample_by_feature

    row_sums = sample_by_feature.sum(axis=1)
    if (row_sums <= 0).any():
        bad_samples = sample_by_feature.index[row_sums <= 0].tolist()
        raise ValueError(f"Samples with zero total abundance cannot be transformed: {', '.join(bad_samples)}")

    rel = sample_by_feature.div(row_sums, axis=0)
    if transform == "relative_abundance":
        return rel
    if transform == "hellinger":
        return np.sqrt(rel)
    raise ValueError(f"Unsupported transform: {transform}")


def finite_persistence_rows(diagrams: list[np.ndarray]) -> list[dict[str, float | int | str]]:
    rows: list[dict[str, float | int | str]] = []
    for dim, diagram in enumerate(diagrams):
        for birth, death in diagram:
            if math.isinf(death):
                persistence = np.nan
                death_value: float | str = "inf"
            else:
                persistence = float(death - birth)
                death_value = float(death)
            rows.append(
                {
                    "dimension": dim,
                    "birth": float(birth),
                    "death": death_value,
                    "persistence": persistence,
                }
            )
    return rows


def plot_barcode(diagrams: list[np.ndarray], output_path: Path) -> None:
    fig, axes = plt.subplots(len(diagrams), 1, figsize=(8, max(3, 2.4 * len(diagrams))), squeeze=False)
    for dim, diagram in enumerate(diagrams):
        ax = axes[dim, 0]
        finite = diagram[np.isfinite(diagram[:, 1])]
        finite = finite[np.argsort(finite[:, 0])]
        for i, (birth, death) in enumerate(finite):
            ax.plot([birth, death], [i, i], linewidth=1.8)
        ax.set_title(f"H{dim} barcode")
        ax.set_xlabel("Filtration value")
        ax.set_ylabel("Feature")
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def plot_mds(
    distance_matrix: np.ndarray,
    metadata: pd.DataFrame,
    grouping_column: str,
    output_path: Path,
    random_state: int,
) -> pd.DataFrame:
    sample_ids = metadata.index.to_list()
    if len(sample_ids) < 3:
        raise ValueError("At least three samples are needed for an MDS plot.")

    mds = MDS(
        n_components=2,
        dissimilarity="precomputed",
        normalized_stress="auto",
        random_state=random_state,
    )
    coords = mds.fit_transform(distance_matrix)
    coord_table = pd.DataFrame(coords, columns=["MDS1", "MDS2"], index=sample_ids)
    coord_table.index.name = "sample_id"
    coord_table = coord_table.join(metadata, how="left")

    color_values = coord_table[grouping_column] if grouping_column in coord_table.columns else coord_table.index
    groups = pd.Series(color_values, index=coord_table.index, dtype="category")
    palette = plt.get_cmap("tab20")
    category_codes = groups.cat.codes

    fig, ax = plt.subplots(figsize=(6, 5))
    scatter = ax.scatter(
        coord_table["MDS1"],
        coord_table["MDS2"],
        c=category_codes,
        cmap=palette,
        s=55,
        edgecolor="black",
        linewidth=0.4,
    )
    ax.set_xlabel("MDS1")
    ax.set_ylabel("MDS2")
    ax.set_title("TDA input space by microbiome distance")

    categories = list(groups.cat.categories)
    handles = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            label=str(category),
            markerfacecolor=palette(i % palette.N),
            markeredgecolor="black",
            markersize=7,
        )
        for i, category in enumerate(categories)
    ]
    if handles:
        ax.legend(handles=handles, title=grouping_column, bbox_to_anchor=(1.04, 1), loc="upper left")
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)
    return coord_table


def main() -> None:
    args = parse_args()
    feature_table_path = Path(args.feature_table)
    distance_matrix_path = Path(args.distance_matrix) if args.distance_matrix else None
    metadata_path = Path(args.metadata)
    out_tables = Path(args.out_tables)
    out_figures = Path(args.out_figures)
    out_tables.mkdir(parents=True, exist_ok=True)
    out_figures.mkdir(parents=True, exist_ok=True)

    if distance_matrix_path is not None:
        distance_df = read_distance_matrix(distance_matrix_path)
        sample_ids = distance_df.index.to_list()
        metadata = metadata_for_samples(
            metadata_path=metadata_path,
            sample_ids=sample_ids,
            grouping_column=args.grouping_column,
            allow_sample_only=args.metadata == DEFAULT_METADATA,
        )
        common_samples = [sample for sample in sample_ids if sample in metadata.index]
        if len(common_samples) < 3:
            raise ValueError("TDA requires at least three samples.")
        distance_df = distance_df.loc[common_samples, common_samples]
        feature_count: int | str = "not_applicable_precomputed_distance"
        transform_label = "precomputed_distance"
    else:
        feature_table = read_feature_table(feature_table_path)
        metadata = read_metadata(metadata_path)

        common_samples = [sample for sample in feature_table.columns if sample in metadata.index]
        if len(common_samples) < 3:
            raise ValueError("TDA requires at least three samples shared by feature table and metadata.")

        feature_table = feature_table[common_samples]
        metadata = metadata.loc[common_samples]
        feature_table = filter_features(feature_table, args.min_feature_total_count, args.min_feature_prevalence)

        sample_by_feature = feature_table.T
        transformed = transform_table(sample_by_feature, args.transform)

        distance_matrix = pairwise_distances(transformed, metric=args.metric)
        distance_df = pd.DataFrame(distance_matrix, index=transformed.index, columns=transformed.index)
        feature_count = feature_table.shape[0]
        transform_label = args.transform

    distance_df.index.name = "sample_id"
    distance_df.to_csv(out_tables / "tda_distance_matrix.tsv", sep="\t")
    distance_matrix = distance_df.values

    threshold = np.inf if args.threshold is None else args.threshold
    result = ripser(
        distance_matrix,
        distance_matrix=True,
        maxdim=args.max_dim,
        thresh=threshold,
    )
    diagrams = result["dgms"]

    persistence_table = pd.DataFrame(finite_persistence_rows(diagrams))
    persistence_table.to_csv(out_tables / "tda_persistence_pairs.tsv", sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(6, 5))
    plot_diagrams(diagrams, ax=ax, show=False)
    ax.set_title("Persistence diagram")
    fig.tight_layout()
    fig.savefig(out_figures / "tda_persistence_diagram.png", dpi=300)
    plt.close(fig)

    plot_barcode(diagrams, out_figures / "tda_barcode.png")

    coord_table = plot_mds(
        distance_matrix=distance_matrix,
        metadata=metadata,
        grouping_column=args.grouping_column,
        output_path=out_figures / "tda_mds.png",
        random_state=args.random_state,
    )
    coord_table.to_csv(out_tables / "tda_mds_coordinates.tsv", sep="\t")

    summary = pd.DataFrame(
        {
            "metric": [args.metric],
            "transform": [transform_label],
            "samples": [len(common_samples)],
            "features_after_filtering": [feature_count],
            "max_dimension": [args.max_dim],
            "threshold": [args.threshold if args.threshold is not None else "none"],
        }
    )
    summary.to_csv(out_tables / "tda_summary.tsv", sep="\t", index=False)

    print(f"TDA complete for {len(common_samples)} samples.")
    print(f"Wrote {out_tables / 'tda_persistence_pairs.tsv'}")
    print(f"Wrote {out_figures / 'tda_persistence_diagram.png'}")


if __name__ == "__main__":
    main()
