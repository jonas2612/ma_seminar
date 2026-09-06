import gzip
import logging
import os
import re
import shutil
import subprocess
import tarfile
import tempfile
import zipfile

from pathlib import Path
from typing import List, Optional
from matplotlib import pyplot as plt
from scipy.stats import median_abs_deviation
from urllib.request import urlopen

import anndata as ad
import decoupler as dc
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from tqdm import tqdm

ad.settings.allow_write_nullable_strings = True

def geo_family(gse_id):
    return gse_id[:-3] + "nnn"


def geo_suppl_url(gse_id):
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{geo_family(gse_id)}/{gse_id}/suppl"


def download_file(url, output_path, logger: logging.Logger):
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if output_path.exists() and output_path.stat().st_size > 0:
        logger.warning(f"Already exists: {output_path}")
        return True

    logger.info(f"Downloading: {url}")

    try:
        subprocess.run(
            ["curl", "-L", "-s", "-C", "-", url, "-o", str(output_path)],
            check=True
        )
        return output_path.exists() and output_path.stat().st_size > 0

    except subprocess.CalledProcessError:
        if output_path.exists():
            output_path.unlink()
        return False


def is_valid_tar(path):
    try:
        return tarfile.is_tarfile(path)
    except Exception:
        return False


def exctract_archive(arch_path, out_dir = None):
    archive_path = Path(arch_path)
    if (archive_path.name.endswith(".tar.gz")):
        name = archive_path.name[:-7]
    elif (archive_path.name.endswith(".tgz")):
        name = archive_path.name[:-4]
    elif (archive_path.name.endswith(".tar")):
        name = archive_path.name[:-4]
    elif (archive_path.name.endswith(".zip")):
        name = archive_path.name[:-4]
    else:
        return None

    outdir = archive_path.parent / name if out_dir is None else out_dir
    outdir.mkdir(parents=True, exist_ok=True)

    if archive_path.name.endswith((".tar", ".tar.gz", ".tgz")):
        with tarfile.open(archive_path, "r:*") as tar:
            tar.extractall(outdir)
    elif archive_path.name.endswith(".zip"):
        with zipfile.ZipFile(archive_path, "r") as zip_file:
            zip_file.extractall(outdir)
    return outdir


def convert_genes_to_features(sample_out, logger: logging.Logger):
    genes_path = sample_out / "genes.tsv.gz"
    features_path = sample_out / "features.tsv.gz"

    if not genes_path.exists():
        return

    needs_conversion = not features_path.exists()

    if features_path.exists():
        try:
            with gzip.open(features_path, "rt") as f:
                line_count = sum(1 for _ in f)
            if line_count < 1000:
                needs_conversion = True
        except Exception:
            needs_conversion = True

    if not needs_conversion:
        return

    with gzip.open(genes_path, "rt") as infile, gzip.open(features_path, "wt") as outfile:
        for line in infile:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                outfile.write(f"{parts[0]}\t{parts[1]}\tGene Expression\n")

    logger.info(f"Converted genes.tsv.gz to features.tsv.gz in {sample_out}")

def determine_level(dataset_dir: Path, logger: logging.Logger):
    rel_suffixes = ("barcodes.tsv.gz", "features.tsv.gz", "genes.tsv.gz", "genes.tsv.gz", "matrix.mtx.gz")
    arch_suffixes = (".tar", ".tar.gz", ".tgz", ".zip")

    for path in list(dataset_dir.rglob("*")):
        if not path.is_file():
            continue
        if path.name.endswith(arch_suffixes):
            logger.info(f"Extracting {path}")
            exctract_archive(path)

    found_files = []

    for path in dataset_dir.rglob("*"):
        if not path.is_file():
            continue
        if not path.name.endswith(rel_suffixes):
            continue
        if path.parent == dataset_dir:
            continue
        rel_path = path.relative_to(dataset_dir)
        out_name = "_".join(rel_path.parts)
        output_path = dataset_dir / out_name

        if output_path.exists():
            logger.warning(f"Already exists: {output_path.name}")
            found_files.append(output_path)
            continue

        shutil.copy2(path, output_path)
        found_files.append(output_path)
    return found_files


def organize_dataset(dataset_dir, dataset_id, out_dir, logger: logging.Logger):
    """
    dataset_dir: directory containing the contents of downloaded files.
    dataset_id: Dataset identifier
    out_dir: output directory. Here a folder for the dataset_id gets created
    """
    dataset_out = out_dir / dataset_id
    dataset_out.mkdir(parents=True, exist_ok=True)

    files = list(dataset_dir.glob("*.gz"))
    sample_prefixes = set()

    suffixes = [
            "barcodes.tsv.gz",
            "features.tsv.gz",
            "genes.tsv.gz",
            "matrix.mtx.gz",
        ]

    for f in files:
        for suffix in suffixes:
            if f.name.endswith(suffix):
                prefix = f.name.removesuffix(suffix)
                sample_prefixes.add(prefix)
                break

    for prefix in sorted(sample_prefixes):
        sample_out = dataset_out / prefix
        sample_out.mkdir(parents=True, exist_ok=True)

        barcodes = dataset_dir / f"{prefix}_barcodes.tsv.gz"
        features = dataset_dir / f"{prefix}_features.tsv.gz"
        genes = dataset_dir / f"{prefix}_genes.tsv.gz"
        matrix = dataset_dir / f"{prefix}_matrix.mtx.gz"

        if not matrix.exists():
            logger.warning(f"Skipped {prefix}: matrix missing")
            continue
        
        if not barcodes.exists():
            logger.warning(f"Skipped {prefix}: barcodes missing")
            continue
        
        if not features.exists() and not genes.exists():
            logger.warning(f"Skipped {prefix}: features/genes missing")
            continue

        shutil.copy2(matrix, sample_out / matrix.name)
        shutil.copy2(barcodes, sample_out / barcodes.name)
        if features.exists():
            shutil.copy2(features, sample_out / features.name) 
        if genes.exists():
            shutil.copy2(genes, sample_out / genes.name) 
        logger.info(f"Organized {dataset_id}/{prefix}")



def list_supplementary_files(gse_id, logger: logging.Logger):
    url = geo_suppl_url(gse_id) + "/"

    try:
        html = urlopen(url).read().decode("utf-8")
    except Exception as e:
        logger.error(f"Could not read supplementary page for {gse_id}: {e}")
        return []

    files = re.findall(r'href="([^"]+)"', html)

    return [
        f for f in files
        if f.endswith((".tsv.gz", ".mtx.gz"))
    ]


def download_individual_10x_files(gse_id, out_dir, logger: logging.Logger):
    files = list_supplementary_files(gse_id, logger)

    if not files:
        logger.warning(f"No individual 10x-like files found for {gse_id}")
        return

    dataset_out = out_dir / gse_id
    dataset_out.mkdir(parents=True, exist_ok=True)

    prefixes = set()

    for filename in files:
        for suffix in [
            "barcodes.tsv.gz",
            "features.tsv.gz",
            "genes.tsv.gz",
            "matrix.mtx.gz",
        ]:
            if filename.endswith(suffix):
                prefixes.add(filename.replace(suffix, ""))

    for prefix in sorted(prefixes):
        sample_out = dataset_out / prefix
        sample_out.mkdir(parents=True, exist_ok=True)

        mapping = {
            f"{prefix}barcodes.tsv.gz": "_barcodes.tsv.gz",
            f"{prefix}features.tsv.gz": "_features.tsv.gz",
            f"{prefix}genes.tsv.gz": "_genes.tsv.gz",
            f"{prefix}matrix.mtx.gz": "_matrix.mtx.gz",
        }

        downloaded = 0

        for source_name, target_name in mapping.items():
            if source_name in files:
                url = f"{geo_suppl_url(gse_id)}/{source_name}"
                output_path = sample_out / target_name

                if download_file(url, output_path, logger):
                    downloaded += 1

        convert_genes_to_features(sample_out, logger)
        logger.info(f"Downloaded {downloaded} files for {gse_id}/{prefix}")


def process_gse(gse_id, data_dir, out_dir, logger: logging.Logger):
    raw_tar_url = f"{geo_suppl_url(gse_id)}/{gse_id}_RAW.tar"
    raw_tar_path = data_dir / f"{gse_id}_RAW.tar"

    logger.info(f"\nProcessing {gse_id}")

    downloaded_tar = download_file(raw_tar_url, raw_tar_path, logger)

    if downloaded_tar and is_valid_tar(raw_tar_path):
        logger.info(f"extracted to {exctract_archive(raw_tar_path, out_dir / gse_id)}")
        organize_dataset(data_dir, gse_id, out_dir, logger)
        return

    if raw_tar_path.exists():
        raw_tar_path.unlink()

    logger.warning(f"No valid RAW.tar for {gse_id}. Trying individual supplementary files.")
    download_individual_10x_files(gse_id, out_dir, logger)


def load_enabled_gse_ids(data_dir):
    ids_path = data_dir / "gse_datasets.csv"

    df = pd.read_csv(ids_path)

    df["enabled"] = (
        df["enabled"]
        .astype(str)
        .str.lower()
        .isin(["true", "1", "yes"])
    )

    df = df[df["enabled"]]

    return [
        str(dataset_id).strip()
        for dataset_id in df["dataset_id"].dropna()
    ]

def read_raw_data(path: Path, metadata: pd.DataFrame, logger: logging.Logger, save_path: Optional[Path | str] = None) -> ad.AnnData:
    path = Path(path)
    if not path.name.endswith("matrix.mtx.gz") and not path.name.endswith(".h5"):
        raise ValueError("The path doesn't point to a 10x formatted file.")
    if "GSM" not in str(path):
        raise NotImplementedError("This is currently restricted to GEO files.")

    if path.name.endswith("matrix.mtx.gz"):
        pre = path.name[:-13]
        logger.info(f"Loading from {path}. Assume prefix {pre}")
        adata = sc.read_10x_mtx(path.parent, prefix=pre)
    else:
        adata = sc.read_10x_h5(path)
    logger.info("Adding metadata")
    meta = metadata[metadata['sample'].isin(path.name.split("_"))]
    assert len(meta) == 1, "Identified metadata contains more than one rows."
    for col in meta.columns:
        adata.obs[col] = meta.loc[:, col].values[0]
        if col in ["accession", "chip", "species", "sample", "cell_type", "location", "sex", "diet", "medication", "KO", "combined_condition", "Original Name"]:
            adata.obs[col] = adata.obs[col].astype(str)
        elif col in ["age", "sample_nr"]:
            adata.obs[col] = adata.obs[col].astype("Int64")
        elif col in ["symptomatic_atherosclerosis"]:
            adata.obs[col] = adata.obs[col].astype(bool)
        else:
            adata.obs[col] = adata.obs[col].astype(object)

    if save_path is not None:
        adata.write(Path(save_path))
    return adata 

def qc_statistical(adata: ad.AnnData, logger: logging.Logger, save_path: Optional[Path | str] = None) -> ad.AnnData:
    def is_outlier(adata, metric: str, nmads: int):
        M = adata.obs[metric]
        outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
            np.median(M) + nmads * median_abs_deviation(M) < M
        )
        return outlier
    # agnostic for Human and Mouse
    adata.var["mt"] = adata.var_names.str.startswith(("mt-", "MT-"))
    adata.var["ribo"] = adata.var_names.str.startswith(("rps", "rpl", "RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.contains(r"^hb[abdegmqz](?:[-][a-z0-9]+|[0-9]*)?$", case=False, regex=True)
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True
    )
    adata = adata[adata.obs['total_counts'] >= 50]
    adata.obs["outlier"] = (
        is_outlier(adata, "log1p_total_counts", 5)
        | is_outlier(adata, "log1p_n_genes_by_counts", 5)
        | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
    )
    adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (
        adata.obs["pct_counts_mt"] > 8
    )
    logger.info(f"Total number of cells: {adata.n_obs}")
    adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()
    logger.info(f"Number of cells after filtering of low quality cells: {adata.n_obs}")
    if save_path is not None:
        logger.info(f"saving adata in {save_path}")
        adata.write(save_path)
    return adata

def doublet_detection(adata: ad.AnnData, logger: logging.Logger, save_path: Optional[Path | str] = None) -> ad.AnnData:
    logger.info("Starting doublet detection.")
    sc.external.pp.scrublet(adata) #scanpy 1.9.6 convention
    logger.info(f"Number of cells after filtering doublets: {adata[~adata.obs['predicted_doublet']].n_obs}")
    adata = adata[~adata.obs['predicted_doublet']]
    if save_path is not None:
        logger.info(f"saving adata in {save_path}")
        adata.write(save_path)
    return adata

def combine_samples(adatas: List[ad.AnnData], logger: logging.Logger, save_path: Optional[Path | str] = None) -> ad.AnnData:
    logger.info("combining samples")
    adata = ad.concat(adatas)
    adata.obs_names_make_unique() #important for cell type annotation
    adata.layers['counts'] = adata.X.copy()
    if save_path is not None:
        logger.info(f"saving adata in {save_path}")
        adata.write(save_path)
    return adata

def normalize_data(adata: ad.AnnData, logger: logging.Logger, save_path: Optional[Path | str] = None) -> ad.AnnData:
    logger.info("running normalization")
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    if save_path is not None:
        logger.info(f"saving adata in {save_path}")
        adata.write(save_path)
    return adata

def pca(adata: ad.AnnData, logger: logging.Logger, n_comps: int = 50, save_path: Optional[Path | str] = None) -> ad.AnnData:
    logger.info(f"running pca computation with {n_comps} components")
    sc.pp.pca(adata, n_comps=n_comps)
    if save_path is not None:
        logger.info(f"saving adata in {save_path}")
        adata.write(save_path)
    return adata

def neighbors(adata: ad.AnnData, logger: logging.Logger, n_neighbors: int = 15, save_path: Optional[Path | str] = None) -> ad.AnnData:
    logger.info(f"running neighbor computation with {n_neighbors} neighbors")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors)
    if save_path is not None:
        logger.info(f"saving adata in {save_path}")
        adata.write(save_path)
    return adata

def umap(adata: ad.AnnData, logger: logging.Logger, save_path: Optional[Path | str] = None, save_path_plot: Optional[Path | str] = None) -> ad.AnnData:
    logger.info("running UMAP")
    sc.tl.umap(adata)
    if save_path is not None:
        logger.info(f"saving adata in {save_path}")
        adata.write(save_path)
    if save_path_plot is not None:
        logger.info(f"saving UMAP plot in {save_path_plot}")
        ax = sc.pl.umap(adata, show=False)
        plt.savefig(save_path_plot)
        plt.close()
    return adata

def annotate_adata(
        adata: ad.AnnData, 
        logger: logging.Logger,
        path_to_annotation_script,
        max_split_size: int = 50_000, 
        output_dir_images: Optional[Path | str] = None, 
        save_path: Optional[Path | str] = None,
        save_path_plot: Optional[Path | str] = None
    ) -> ad.AnnData:
    assert adata.obs.index.is_unique, "adata index is not unique. This is necessary for the annotation."
    n_splits = np.ceil(adata.n_obs/max_split_size)
    split_size = int(np.floor(adata.n_obs/n_splits))

    logger.info(f"creating temporary directory in {os.getcwd()}")
    tmp_path = Path(f"{os.getcwd()}/.tmp")
    tmp_path.mkdir()

    new_obs = []
    for enum, i in tqdm(enumerate(range(0, adata.n_obs, split_size))):
        curr_save_path = tmp_path / f"adata_{enum}.h5ad"
        curr_output_dir = tmp_path / f"outdir_{enum}"
        adata[i:(i+split_size)].write(curr_save_path)
        logger.info("Running annotation ...")
        if output_dir_images is None:
            subprocess.run([
                "python", 
                str(path_to_annotation_script),
                str(curr_save_path),
                str(curr_output_dir),
                "--lognorm_bool"
            ], check=True)
        else:
            subprocess.run([
                        "python", 
                        str(path_to_annotation_script),
                        str(curr_save_path),
                        str(curr_output_dir),
                        "output_dir_images", 
                        str(output_dir_images), 
                        "--lognorm_bool"
                    ], check=True)
        new_obs.append(sc.read_h5ad(curr_output_dir / "full_level1.h5ad").obs[['cell_type_level1', 'cell_type_uncert']])
    logger.info("Removing temporary directory")
    subprocess.run(["rm", "-r", tmp_path], check=True)
    logger.info("Adding cell type annotation information.")
    new_obs = pd.concat(new_obs)
    logger.debug(f"new obs: {new_obs.columns}")
    adata_obs = adata.obs.copy()
    adata_obs_new = adata_obs.merge(new_obs, left_index=True, right_index=True)
    logger.debug(f"new adata obs: {adata_obs_new.columns}")
    adata.obs = adata_obs_new
    if save_path is not None:
        logger.info(f"saving adata in {save_path}")
        adata.write(save_path)
    if save_path_plot is not None:
        logger.info(f"saving UMAP plot in {save_path_plot}")
        sc.pl.umap(adata, color="cell_type_level1")
        plt.savefig(save_path_plot)
    return adata

def run_limma_dea(
        adata: ad.AnnData,
        cell_type,
        condition_levels,
        output_path,
        logger: logging.Logger, 
        sample_col = "sample", 
        cell_type_col = "cell_type", 
        condition_col="...", 
        layer="counts", 
        mode="sum",
        covariate_cols=None,
        min_cells=20,
        min_samples_per_group=2
        ):
    logger.info("start DEA")
    sample_to_condition = (
        adata.obs[[sample_col, condition_col]]
        .dropna()
        .drop_duplicates()
    )
    condition_lookup = (
        sample_to_condition
        .drop_duplicates(subset=sample_col)
        .set_index(sample_col)[condition_col]
    )
    adata_pb = dc.pp.pseudobulk(
        adata, sample_col=sample_col,
        groups_col=cell_type_col,
        layer=layer,
        mode=mode
    )
    adata_pb.obs[condition_col] = adata_pb.obs[sample_col].map(condition_lookup)
    if "psbulk_n_cells" in adata_pb.obs.columns:
        adata_pb.obs["n_cells"] = adata_pb.obs["psbulk_n_cells"]
    elif "psbulk_n_cells" not in adata_pb.obs.columns:
        n_cells_lookup = (
            adata.obs
            .groupby([sample_col, cell_type_col], observed=True)
            .size()
            .rename("n_cells")
        )

        adata_pb.obs["n_cells"] = [
            n_cells_lookup.get((sample, cell_type), np.nan)
            for sample, cell_type in zip(
                adata_pb.obs[sample_col],
                adata_pb.obs[cell_type_col],
            )
        ]

    covariate_cols = list(covariate_cols or [])
    keep = (
        adata_pb.obs[cell_type_col].astype(str).eq(str(cell_type)) 
        & adata_pb.obs["n_cells"].ge(min_cells)
        & adata_pb.obs[condition_col].isin(condition_levels)
    )
    logger.info(adata_pb)
    logger.info(keep)
    adata_pb_ct = adata_pb[keep].copy()

    meta = adata_pb_ct.obs[[sample_col, cell_type_col, condition_col, "n_cells", *covariate_cols]].copy()
    meta.index = meta.index.astype(str)
    meta.index.name = "pseudobulk_id"

    if sp.issparse(adata_pb_ct.X):
        counts = adata_pb_ct.X.toarray().T
    else:
        counts = np.asarray(adata_pb_ct.X).T

    counts = pd.DataFrame(counts, index=adata_pb_ct.var_names.astype(str), columns=meta.index)
    logger.info(counts.head())

    r_script_path = Path(os.path.dirname(os.path.realpath(__file__))) / "scRNAseq_dea.R"

    with tempfile.TemporaryDirectory(prefix="limma_voom_", dir=Path(__file__).resolve().parent) as tmp_dir_name:
        tmp_dir = Path(tmp_dir_name)
        counts_file = tmp_dir / "pseudobulk_counts.tsv"
        metadata_file = tmp_dir / "pseudobulk_metadata.tsv"

        counts.to_csv(
            counts_file,
            sep="\t",
            index=True,
            index_label="gene",
        )

        meta.to_csv(
            metadata_file,
            sep="\t",
            index=True,
            index_label="pseudobulk_id",
        )

        cmd = [
            "Rscript",
            str(r_script_path),
            "--counts", str(counts_file),
            "--metadata", str(metadata_file),
            "--output", str(output_path),
            "--group-col", condition_col,
            "--group-levels", ",".join(condition_levels),
            "--celltype-col", cell_type_col,
            "--celltype", str(cell_type),
            "--n-cells-col", "n_cells",
            "--min-cells", str(min_cells),
            "--min-samples-per-group", str(min_samples_per_group),
        ]

        if covariate_cols:
            cmd.extend([
                "--covariates",
                ",".join(covariate_cols),
            ])

        try:
            completed = subprocess.run(
                cmd,
                check=True,
                text=True,
                capture_output=True,
            )
        except subprocess.CalledProcessError as exc:
            logger.error("Rscript command failed: %s", " ".join(cmd))
            logger.error("Rscript stdout:\n%s", exc.stdout)
            logger.error("Rscript stderr:\n%s", exc.stderr)
            raise RuntimeError(
                "limma-voom pseudobulk DEA failed. "
                f"R stderr:\n{exc.stderr}"
            ) from exc

    
    results = pd.read_csv(output_path, sep="\t")

    return results.sort_values(
        ["adj.P.Val", "P.Value"],
        ascending=True,
        kind="stable",
    ).reset_index(drop=True)