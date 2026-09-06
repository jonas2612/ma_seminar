from os import listdir

from code.scRNAseq_pipeline import *
import pandas as pd
import logging
from tqdm import tqdm

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s [%(levelname)s]: %(message)s"
)
logger = logging.getLogger("shape_mapping")

base_dir = Path("/usr/local/storage/data_scRNAseq")
umap_dir = base_dir / "umaps"
dea_dir = base_dir / "dea"
annot_dir = base_dir / "annotated_data"

datasets_and_contrasts = {
    "GSE309462": {"cond_col": "symptomatic_atherosclerosis", "contrast": [("FALSE", "TRUE")]},
    "GSE260657": {"cond_col": "symptomatic_atherosclerosis", "contrast": [("FALSE", "TRUE")]},
    "GSE253903": {"cond_col": "symptomatic_atherosclerosis", "contrast": [("FALSE", "TRUE")]},
    "GSE159677": {"cond_col": "cell_type", "contrast": [("plaque adjacent", "plaque")]},
}

for key, value in datasets_and_contrasts.items():
    process_gse(key, base_dir / "raw_data_download", base_dir / "raw_data", logger)
    data_dir = base_dir / "raw_data" / key
    data_dirs = determine_level(data_dir, logger)
    adatas = []
    metadata = pd.read_csv(r"/home/f/flor/metadata_all_samples.txt", sep = "\t")
    for path in tqdm([x for x in listdir(data_dir) if x.endswith("matrix.mtx.gz") and "GSM" in x]):
        logger.info(f"Working on {path}")
        adata = read_raw_data(data_dir / path, metadata, logger)
        adata = qc_statistical(adata, logger)
        adata = doublet_detection(adata, logger)
        adatas.append(adata)
    adata = combine_samples(adatas, logger)
    adata = normalize_data(adata, logger)
    adata = pca(adata, logger)
    adata = neighbors(adata, logger)
    adata = umap(adata, logger, save_path_plot = umap_dir / f"{key}_umap.png")
    adata = annotate_adata(adata, logger, "/home/f/flor/gitrepos/plaque-atlas-mapping/automatic_mapping_level1.py",
                        save_path_plot = umap_dir / f"{key}_umap_annot.png", save_path=annot_dir / f"{key}_annot.h5ad")
    for ct in adata.obs['cell_type_level1'].unique():
        logger.info(f"Working on {ct}")
        for contr in value['contrast']:
            run_limma_dea(adata, ct, contr, dea_dir / f"{key}_{ct}_{value['cond_col']}_{contr[1]}-{contr[0]}.tsv",
                      condition_col=value['cond_col'], logger=logger, cell_type_col="cell_type_level1")
    