#!/usr/bin/env python3
"""
test_time_split_FC_v2.py

Refactored MEANtools metafingerprint script.

Major changes compared with the original script:
  1. Keeps the original command-line flags (-ft and -qm), but internally uses only -ft.
  2. Treats the input as one integrated feature table, not separate transcript/metabolite tables.
  3. Renames the first identifier column to 'feature' if needed.
  4. Uses a fast feature lookup instead of repeatedly scanning data frames.
  5. Robustly parses cluster members separated by spaces, commas, semicolons, or brackets.
  6. Skips empty/problematic clusters instead of crashing during SVD.
  7. Actually uses the normalized feature matrix for fingerprint calculation.
  8. Filters cluster tables by the requested decay rate when DR labels are available.

Example:
python test_time_split_FC_v2.py \
  -ft trial5_allsamples_batch_quantification_pos_label.csv \
  -qm trial5_allsamples_batch_quantification_pos_label_2.csv \
  -dr 25 \
  -dn nwogroot_coabun.sqlite \
  -g timepoint \
  -m sample_info_metabolomics.csv \
  -o correlation
"""

import os
import re
import sys
import argparse
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy.stats import zscore

import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import to_hex

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import (
    KMeans,
    DBSCAN,
    AgglomerativeClustering,
    SpectralClustering,
    Birch,
    AffinityPropagation,
    MeanShift,
    estimate_bandwidth,
)
from sklearn.mixture import GaussianMixture

import gizmos


# -----------------------------------------------------------------------------
# Arguments
# -----------------------------------------------------------------------------


def get_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Calculate SVD-based metafingerprints from an integrated feature table."
    )
    parser.add_argument(
        "-ft", "--feature_table",
        required=True,
        help="Integrated feature table. First column should be feature/gene/metabolite ID."
    )
    parser.add_argument(
        "-qm", "--quantitation_matrix",
        default="",
        required=False,
        help="Deprecated/backward-compatible flag. This script ignores -qm and uses only -ft."
    )
    parser.add_argument(
        "-dr", "--decay_rate",
        default=25,
        type=int,
        required=False,
        help="Decay rate to filter DR tables/clusters. Default: 25."
    )
    parser.add_argument(
        "-g", "--groups",
        default="treatment",
        choices=["treatment", "timepoint"],
        help="Group samples by treatment or timepoint."
    )
    parser.add_argument(
        "-m", "--metadata",
        required=True,
        help="Metadata file with at least sample and time/timepoint or condition/treatment columns."
    )
    parser.add_argument(
        "-o", "--outfile",
        default="results",
        help="Output prefix. Default: results."
    )
    parser.add_argument(
        "-dn", "--sqlite_db_name",
        required=True,
        help="SQLite database containing cluster tables."
    )
    parser.add_argument(
        "--no_zscore",
        action="store_true",
        help="Do not z-score features before SVD. By default, row-wise feature z-scoring is applied."
    )
    parser.add_argument(
        "--cluster_method",
        default="kmeans",
        choices=["kmeans", "dbscan", "agglomerative", "gmm", "spectral", "meanshift", "birch", "affinity"],
        help="Clustering method for PCA coordinates. Default: kmeans."
    )
    parser.add_argument(
        "--n_clusters",
        default=3,
        type=int,
        help="Number of PCA clusters for applicable methods. Default: 3."
    )
    parser.add_argument(
        "--plot_format",
        default="png",
        choices=["png", "pdf", "svg"],
        help="Format for PCA plots. Default: png."
    )
    return parser.parse_args()


# -----------------------------------------------------------------------------
# Utility functions
# -----------------------------------------------------------------------------


def natural_key(value: object) -> List[object]:
    """Natural sort key: 2 before 10, T2 before T10."""
    text = str(value)
    return [int(x) if x.isdigit() else x.lower() for x in re.split(r"(\d+)", text)]


def find_column(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    """Return first matching column name from a candidate list, case-insensitive."""
    lower_map = {str(c).lower(): c for c in df.columns}
    for c in candidates:
        if c.lower() in lower_map:
            return lower_map[c.lower()]
    return None


def safe_output_prefix(prefix: object) -> str:
    prefix = str(prefix) if prefix else "results"
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", prefix)


def orient_vector(vec: np.ndarray) -> np.ndarray:
    """
    SVD vectors have arbitrary sign. Orient the vector so that the largest
    absolute entry is positive. This makes outputs more consistent.
    """
    if vec.size == 0:
        return vec
    idx = int(np.nanargmax(np.abs(vec)))
    if vec[idx] < 0:
        vec = -vec
    return vec


# -----------------------------------------------------------------------------
# Input data
# -----------------------------------------------------------------------------


def read_feature_table(feature_table: str, quantitation_matrix: str = "") -> pd.DataFrame:
    """
    Read only the -ft file. The -qm flag is kept for compatibility but ignored.
    The first identifier column is normalized to 'feature'.
    """
    if quantitation_matrix:
        try:
            ft_abs = os.path.abspath(feature_table)
            qm_abs = os.path.abspath(quantitation_matrix)
            if ft_abs != qm_abs:
                print(
                    "WARNING: -qm/--quantitation_matrix was supplied but is ignored in v2.\n"
                    f"         Using only -ft/--feature_table: {feature_table}\n"
                    f"         Ignored -qm file: {quantitation_matrix}"
                )
        except Exception:
            print("WARNING: -qm/--quantitation_matrix is ignored in v2.")

    df = pd.read_csv(feature_table)
    df.columns = [str(c).strip() for c in df.columns]

    if df.shape[1] < 2:
        raise ValueError("Feature table must contain one feature column and at least one sample column.")

    # Normalize ID column to 'feature'
    if "feature" not in df.columns:
        if "gene" in df.columns:
            df = df.rename(columns={"gene": "feature"})
        elif "metabolite" in df.columns:
            df = df.rename(columns={"metabolite": "feature"})
        else:
            first_col = df.columns[0]
            print(f"WARNING: No 'feature', 'gene', or 'metabolite' column found. Using first column '{first_col}' as feature ID.")
            df = df.rename(columns={first_col: "feature"})

    df["feature"] = df["feature"].astype(str).str.strip()
    df = df[df["feature"].notna() & (df["feature"] != "")].copy()

    # Coerce sample columns to numeric where possible.
    for col in df.columns:
        if col != "feature":
            df[col] = pd.to_numeric(df[col], errors="coerce")

    # Remove completely non-numeric columns apart from feature.
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    keep_cols = ["feature"] + numeric_cols
    dropped = [c for c in df.columns if c not in keep_cols]
    if dropped:
        print(f"WARNING: Dropping non-numeric columns from feature table: {dropped}")
    df = df[keep_cols]

    # Collapse duplicate feature IDs by averaging numeric columns.
    if df["feature"].duplicated().any():
        n_dup_rows = int(df["feature"].duplicated().sum())
        n_dup_ids = int(df.loc[df["feature"].duplicated(), "feature"].nunique())
        print(f"WARNING: Found {n_dup_rows} duplicated rows across {n_dup_ids} duplicated feature IDs. Averaging duplicates.")
        df = df.groupby("feature", as_index=False).mean(numeric_only=True)

    print(f"Loaded integrated feature table: {df.shape[0]} features x {df.shape[1] - 1} numeric sample columns")
    return df


def normalize_feature_table(df: pd.DataFrame, do_zscore: bool = True) -> pd.DataFrame:
    """
    Row-wise z-score each feature across all numeric samples.
    This makes SVD fingerprints describe shape rather than absolute magnitude.
    """
    out = df.copy()
    numeric_cols = out.select_dtypes(include=[np.number]).columns.tolist()

    if not numeric_cols:
        raise ValueError("No numeric sample columns found in feature table.")

    if do_zscore:
        values = out[numeric_cols].to_numpy(dtype=float)
        scaled = zscore(values, axis=1, nan_policy="omit")
        scaled = np.nan_to_num(scaled, nan=0.0, posinf=0.0, neginf=0.0)
        out[numeric_cols] = scaled
        print("Applied row-wise z-score normalization across sample columns.")
    else:
        out[numeric_cols] = out[numeric_cols].replace([np.inf, -np.inf], np.nan).fillna(0.0)
        print("Skipped z-score normalization; using numeric values as provided.")

    return out


def read_metadata(metadata_file: str) -> pd.DataFrame:
    metadata = pd.read_csv(metadata_file)
    metadata.columns = [str(c).strip() for c in metadata.columns]

    sample_col = find_column(metadata, ["sample", "Sample", "sample_id", "SampleID", "sample_name"])
    if sample_col is None:
        raise ValueError("Metadata must contain a sample column, e.g. 'sample'.")

    if sample_col != "sample":
        metadata = metadata.rename(columns={sample_col: "sample"})

    metadata["sample"] = metadata["sample"].astype(str).str.strip()
    print(f"Loaded metadata: {metadata.shape[0]} rows")
    return metadata


def get_group_column(metadata: pd.DataFrame, group_flag: str) -> str:
    if group_flag == "timepoint":
        col = find_column(metadata, ["time", "timepoint", "dpi", "Time", "Timepoint"])
        if col is None:
            raise ValueError("For -g timepoint, metadata must contain a time/timepoint column.")
        return col

    col = find_column(metadata, ["condition", "treatment", "Treatment", "group", "Group"])
    if col is None:
        raise ValueError("For -g treatment, metadata must contain a condition/treatment column.")
    return col


def group_samples(metadata: pd.DataFrame, feature_df: pd.DataFrame, group_flag: str) -> Dict[str, pd.DataFrame]:
    """Create one feature matrix per timepoint/treatment group."""
    group_col = get_group_column(metadata, group_flag)
    print(f"Grouping samples by metadata column: {group_col}")

    # Work with strings to avoid mismatches such as 12 vs '12'.
    metadata = metadata.copy()
    metadata[group_col] = metadata[group_col].astype(str).str.strip()

    feature_columns = set(feature_df.columns)
    groups = sorted(metadata[group_col].dropna().unique().tolist(), key=natural_key)
    grouped = {}

    for group in groups:
        subset = metadata[metadata[group_col] == group]
        requested_samples = subset["sample"].astype(str).str.strip().tolist()
        selected_samples = [s for s in requested_samples if s in feature_columns]
        missing_samples = [s for s in requested_samples if s not in feature_columns]

        if missing_samples:
            print(f"WARNING: group {group}: {len(missing_samples)} metadata samples not present in feature table.")

        if not selected_samples:
            print(f"WARNING: group {group}: no matching sample columns; skipping group.")
            continue

        grouped[group] = feature_df[["feature"] + selected_samples].copy()
        print(f"Group {group}: using {len(selected_samples)} samples")

    if not grouped:
        raise ValueError("No groups contained matching samples between metadata and feature table.")

    return grouped


# -----------------------------------------------------------------------------
# Cluster import and parsing
# -----------------------------------------------------------------------------


def load_clusters(sqlite_db_name: str, decay_rate: int) -> pd.DataFrame:
    """Import cluster tables from SQLite and keep clusters matching DR_<decay_rate> when available."""
    print(f"Importing clusters from SQLite database: {sqlite_db_name}")
    all_tables = gizmos.import_from_sql(
        sqlite_db_name,
        sqlite_tablename="",
        df_columns=[],
        conditions={},
        structures=False,
        clone=True,
    )

    if not all_tables:
        raise ValueError(f"No tables were imported from {sqlite_db_name}")

    frames = []
    for table_name, df in all_tables.items():
        df = df.copy()
        match = re.search(r"(DR_\d+)", str(table_name))
        df["Source"] = match.group(1) if match else "unknown"
        df["Table"] = table_name
        frames.append(df)

    frame = pd.concat(frames, axis=0, ignore_index=True)

    if "Cluster" not in frame.columns:
        raise ValueError("Cluster table must contain a 'Cluster' column.")
    if "Members" not in frame.columns:
        raise ValueError("Cluster table must contain a 'Members' column.")

    requested_source = f"DR_{decay_rate}"
    if requested_source in set(frame["Source"]):
        frame = frame[frame["Source"] == requested_source].copy()
        print(f"Filtered clusters to {requested_source}: {frame.shape[0]} clusters")
    else:
        print(f"WARNING: {requested_source} not found in table names. Using all imported clusters: {frame.shape[0]} clusters")

    frame = frame.reset_index(drop=True)
    frame["cluster_id"] = [f"cluster{i:03d}" for i in range(1, frame.shape[0] + 1)]
    frame["ID"] = frame["cluster_id"] + "_" + frame["Source"].astype(str) + "_" + frame["Cluster"].astype(str)

    clusters_frame = frame[["ID", "Members"]].copy()
    print(f"Prepared cluster frame: {clusters_frame.shape[0]} clusters")
    return clusters_frame


def parse_members(member_string: object) -> List[str]:
    """Parse cluster members separated by whitespace, comma, semicolon, or brackets."""
    if pd.isna(member_string):
        return []
    text = str(member_string).strip()
    text = text.replace("[", " ").replace("]", " ").replace("{", " ").replace("}", " ")
    text = text.replace("'", "").replace('"', "")
    members = [x.strip() for x in re.split(r"[,;\s]+", text) if x.strip()]
    return members


# -----------------------------------------------------------------------------
# Fingerprints
# -----------------------------------------------------------------------------


def build_feature_lookup(feature_df: pd.DataFrame) -> pd.DataFrame:
    lookup = feature_df.set_index("feature", drop=True)
    # After duplicate collapsing, index should be unique. If not, keep first as a fallback.
    if not lookup.index.is_unique:
        print("WARNING: feature index is still not unique after preprocessing. Keeping first occurrence.")
        lookup = lookup[~lookup.index.duplicated(keep="first")]
    return lookup


def extract_cluster_matrix(row: pd.Series, feature_lookup: pd.DataFrame) -> Tuple[str, pd.DataFrame, List[str]]:
    cluster_id = row["ID"]
    members = parse_members(row["Members"])

    found = []
    missing = []

    for member in members:
        if member in feature_lookup.index:
            found.append(member)
        else:
            missing.append(member)

    if not found:
        return cluster_id, pd.DataFrame(), missing

    extracted = feature_lookup.loc[found].copy()
    extracted.insert(0, "feature", found)
    return cluster_id, extracted, missing


def calculate_sigma(cluster_df: pd.DataFrame) -> Optional[np.ndarray]:
    """Calculate first right singular vector of cluster feature x sample matrix."""
    if cluster_df.empty:
        return None

    numeric_cols = cluster_df.select_dtypes(include=[np.number]).columns.tolist()
    if not numeric_cols:
        return None

    X = cluster_df[numeric_cols].to_numpy(dtype=float)
    X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)

    # Remove all-zero rows because they do not contribute to the fingerprint.
    keep_rows = np.linalg.norm(X, axis=1) > 0
    X = X[keep_rows]

    if X.shape[0] == 0 or X.shape[1] == 0:
        return None

    if X.shape[0] == 1:
        return orient_vector(X.flatten())

    try:
        _, _, VT = np.linalg.svd(X, full_matrices=False)
    except np.linalg.LinAlgError:
        return None

    if VT.shape[0] == 0:
        return None

    return orient_vector(VT[0])


def calc_metafingerprint(row: pd.Series, feature_lookup: pd.DataFrame) -> Tuple[Optional[dict], List[dict]]:
    cluster_id, cluster_df, missing = extract_cluster_matrix(row, feature_lookup)

    missing_records = [
        {"ID": cluster_id, "missing_feature": m}
        for m in missing
    ]

    sigma = calculate_sigma(cluster_df)
    if sigma is None:
        return None, missing_records

    numeric_cols = cluster_df.select_dtypes(include=[np.number]).columns.tolist()
    if len(numeric_cols) != len(sigma):
        # This should not happen, but keep a guard for safety.
        return None, missing_records

    record = {"ID": cluster_id, "n_features_found": int(cluster_df.shape[0]), "n_features_missing": int(len(missing))}
    for col, value in zip(numeric_cols, sigma):
        record[col] = float(value)

    return record, missing_records


def calculate_fingerprints_by_group(grouped_data: Dict[str, pd.DataFrame], clusters_frame: pd.DataFrame) -> Tuple[Dict[str, pd.DataFrame], pd.DataFrame]:
    fingerprints_by_group = {}
    all_missing_records = []

    for group, group_df in grouped_data.items():
        print("\n" + "=" * 70)
        print(f"Processing group: {group}")
        print("=" * 70)

        feature_lookup = build_feature_lookup(group_df)
        results = []

        for _, row in tqdm(clusters_frame.iterrows(), total=clusters_frame.shape[0], desc=f"Processing Group {group}"):
            fp, missing_records = calc_metafingerprint(row, feature_lookup)
            all_missing_records.extend([{**r, "Group": group} for r in missing_records])
            if fp is not None:
                results.append(fp)

        if not results:
            print(f"WARNING: No valid fingerprints were generated for group {group}.")
            continue

        fingerprints_df = pd.DataFrame(results)
        fingerprints_by_group[group] = fingerprints_df
        print(f"Generated fingerprints for group {group}: {fingerprints_df.shape[0]} clusters")

    missing_df = pd.DataFrame(all_missing_records)
    return fingerprints_by_group, missing_df


# -----------------------------------------------------------------------------
# PCA, clustering, and plotting
# -----------------------------------------------------------------------------


def get_fingerprint_numeric_columns(df: pd.DataFrame) -> List[str]:
    meta_cols = {"ID", "Group", "n_features_found", "n_features_missing", "Cluster", "Cluster_Color"}
    return [c for c in df.select_dtypes(include=[np.number]).columns if c not in meta_cols]


def pca_on_fingerprints(df: pd.DataFrame, output_prefix: str, label: str = "") -> Tuple[pd.DataFrame, np.ndarray]:
    numeric_cols = get_fingerprint_numeric_columns(df)
    if len(numeric_cols) < 2:
        raise ValueError("Need at least two numeric fingerprint columns for PCA.")

    X = df[numeric_cols].to_numpy(dtype=float)
    X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)

    n_components = min(5, X.shape[0], X.shape[1])
    if n_components < 1:
        raise ValueError("Not enough data for PCA.")

    X_scaled = StandardScaler().fit_transform(X)
    pca = PCA(n_components=n_components)
    pcs = pca.fit_transform(X_scaled)
    explained = pca.explained_variance_ratio_ * 100

    pca_df = pd.DataFrame(pcs, columns=[f"PC{i+1}" for i in range(n_components)])
    pca_df.insert(0, "ID", df["ID"].values)

    for col in ["Group", "n_features_found", "n_features_missing"]:
        if col in df.columns:
            pca_df[col] = df[col].values

    pca_file = f"{output_prefix}{label}_pca.csv"
    pca_df.to_csv(pca_file, index=False)
    print(f"Saved PCA table: {pca_file}")

    return pca_df, explained


def perform_clustering_and_extract_ids(
    df: pd.DataFrame,
    x: int,
    y: int,
    n_clusters: int,
    method: str = "kmeans",
    **kwargs,
) -> pd.DataFrame:
    pcx = f"PC{x}"
    pcy = f"PC{y}"
    if pcx not in df.columns or pcy not in df.columns:
        print(f"WARNING: {pcx}/{pcy} not available; skipping clustering.")
        return df

    data = df[[pcx, pcy]].to_numpy(dtype=float)
    n_rows = data.shape[0]
    if n_rows < 2:
        print("WARNING: fewer than two rows; skipping clustering.")
        return df

    n_clusters = min(max(1, int(n_clusters)), n_rows)

    if method == "kmeans":
        model = KMeans(n_clusters=n_clusters, random_state=42, n_init="auto")
        labels = model.fit_predict(data)
    elif method == "dbscan":
        eps = kwargs.get("eps", 0.5)
        min_samples = kwargs.get("min_samples", 5)
        model = DBSCAN(eps=eps, min_samples=min_samples)
        labels = model.fit_predict(data)
    elif method == "agglomerative":
        model = AgglomerativeClustering(n_clusters=n_clusters)
        labels = model.fit_predict(data)
    elif method == "gmm":
        model = GaussianMixture(n_components=n_clusters, random_state=42)
        labels = model.fit_predict(data)
    elif method == "spectral":
        model = SpectralClustering(n_clusters=n_clusters, affinity="nearest_neighbors", random_state=42)
        labels = model.fit_predict(data)
    elif method == "meanshift":
        bandwidth = estimate_bandwidth(data)
        if bandwidth <= 0 or np.isnan(bandwidth):
            bandwidth = None
        model = MeanShift(bandwidth=bandwidth)
        labels = model.fit_predict(data)
    elif method == "birch":
        model = Birch(n_clusters=n_clusters)
        labels = model.fit_predict(np.ascontiguousarray(data))
    elif method == "affinity":
        damping = kwargs.get("damping", 0.9)
        preference = kwargs.get("preference", -50)
        model = AffinityPropagation(damping=damping, preference=preference, random_state=42)
        labels = model.fit_predict(data)
    else:
        raise ValueError(f"Unknown clustering method: {method}")

    out = df.copy()
    out["Cluster"] = labels

    unique_clusters = sorted(pd.unique(labels), key=natural_key)
    palette = sns.color_palette("hsv", len(unique_clusters))
    cluster_colors = dict(zip(unique_clusters, palette))
    out["Cluster_Color"] = out["Cluster"].map(lambda c: to_hex(cluster_colors[c]))
    return out


def plot_pca(
    pca_df: pd.DataFrame,
    explained: np.ndarray,
    output_file: str,
    title: str,
    color_by: Optional[str] = None,
) -> None:
    if "PC1" not in pca_df.columns:
        print("WARNING: PC1 not found; skipping PCA plot.")
        return

    plt.figure(figsize=(8, 6))

    if "PC2" in pca_df.columns:
        x = pca_df["PC1"]
        y = pca_df["PC2"]
        xlabel = f"PC1 ({explained[0]:.2f}% variance)" if len(explained) > 0 else "PC1"
        ylabel = f"PC2 ({explained[1]:.2f}% variance)" if len(explained) > 1 else "PC2"
    else:
        x = pca_df["PC1"]
        y = np.zeros(pca_df.shape[0])
        xlabel = f"PC1 ({explained[0]:.2f}% variance)" if len(explained) > 0 else "PC1"
        ylabel = "0"

    if color_by and color_by in pca_df.columns:
        categories = sorted(pca_df[color_by].astype(str).unique(), key=natural_key)
        palette = sns.color_palette("tab10", len(categories))
        color_map = dict(zip(categories, palette))
        for cat in categories:
            mask = pca_df[color_by].astype(str) == cat
            plt.scatter(x[mask], y[mask], label=cat, alpha=0.75, s=35)
        plt.legend(title=color_by, bbox_to_anchor=(1.05, 1), loc="upper left")
    else:
        plt.scatter(x, y, alpha=0.75, s=35)

    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()
    print(f"Saved PCA plot: {output_file}")


def combine_groups(fingerprints_by_group: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    frames = []
    for group, df in fingerprints_by_group.items():
        tmp = df.copy()
        tmp["Group"] = group
        frames.append(tmp)
    if not frames:
        raise ValueError("No fingerprint groups to combine.")
    return pd.concat(frames, axis=0, ignore_index=True)


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------


def main(options: argparse.Namespace) -> None:
    output_prefix = safe_output_prefix(options.outfile)

    feature_df = read_feature_table(options.feature_table, options.quantitation_matrix)
    feature_df = normalize_feature_table(feature_df, do_zscore=not options.no_zscore)

    metadata = read_metadata(options.metadata)
    grouped_data = group_samples(metadata, feature_df, options.groups)

    clusters_frame = load_clusters(options.sqlite_db_name, options.decay_rate)
    clusters_file = f"{output_prefix}_clusters_members_DR_{options.decay_rate}.csv"
    clusters_frame.to_csv(clusters_file, index=False)
    print(f"Saved cluster members table: {clusters_file}")

    fingerprints_by_group, missing_df = calculate_fingerprints_by_group(grouped_data, clusters_frame)

    if not missing_df.empty:
        missing_file = f"{output_prefix}_missing_features.csv"
        missing_df.to_csv(missing_file, index=False)
        print(f"Saved missing feature report: {missing_file}")

    # Per-group PCA and clustering.
    for group, fp_df in fingerprints_by_group.items():
        group_safe = safe_output_prefix(group)
        fp_file = f"{output_prefix}_fingerprints_group_{group_safe}.csv"
        fp_df.to_csv(fp_file, index=False)
        print(f"Saved fingerprints for group {group}: {fp_file}")

        try:
            pca_df, explained = pca_on_fingerprints(fp_df, output_prefix, label=f"_group_{group_safe}")
            pca_df = perform_clustering_and_extract_ids(
                pca_df,
                x=1,
                y=2,
                n_clusters=options.n_clusters,
                method=options.cluster_method,
            )
            labelled_file = f"{output_prefix}_pca_group_{group_safe}_labelled.csv"
            pca_df.to_csv(labelled_file, index=False)
            print(f"Saved labelled PCA table for group {group}: {labelled_file}")

            plot_file = f"{output_prefix}_pca_group_{group_safe}.{options.plot_format}"
            plot_pca(
                pca_df,
                explained,
                output_file=plot_file,
                title=f"PCA of metafingerprints: group {group}",
                color_by="Cluster" if "Cluster" in pca_df.columns else None,
            )
        except Exception as exc:
            print(f"WARNING: PCA/clustering failed for group {group}: {exc}")

    # Combined PCA across all groups. This is useful for comparing groups in one PCA space.
    try:
        combined = combine_groups(fingerprints_by_group)
        combined_file = f"{output_prefix}_fingerprints_combined.csv"
        combined.to_csv(combined_file, index=False)
        print(f"Saved combined fingerprints: {combined_file}")

        combined_pca, combined_explained = pca_on_fingerprints(combined, output_prefix, label="_combined")
        combined_pca = perform_clustering_and_extract_ids(
            combined_pca,
            x=1,
            y=2,
            n_clusters=options.n_clusters,
            method=options.cluster_method,
        )
        combined_labelled_file = f"{output_prefix}_combined_pca_labelled.csv"
        combined_pca.to_csv(combined_labelled_file, index=False)
        print(f"Saved combined labelled PCA table: {combined_labelled_file}")

        plot_pca(
            combined_pca,
            combined_explained,
            output_file=f"{output_prefix}_combined_pca_by_group.{options.plot_format}",
            title="Combined PCA of metafingerprints by group",
            color_by="Group" if "Group" in combined_pca.columns else None,
        )
        plot_pca(
            combined_pca,
            combined_explained,
            output_file=f"{output_prefix}_combined_pca_by_cluster.{options.plot_format}",
            title="Combined PCA of metafingerprints by PCA cluster",
            color_by="Cluster" if "Cluster" in combined_pca.columns else None,
        )
    except Exception as exc:
        print(f"WARNING: Combined PCA failed: {exc}")

    print("\nDone.")


if __name__ == "__main__":
    Options = get_args()
    main(Options)
