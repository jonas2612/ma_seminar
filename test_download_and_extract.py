from os import listdir

from code.scRNAseq_pipeline import *
import pandas as pd
import logging
from tqdm import tqdm

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s [%(levelname)s]: %(message)s"
)
logger = logging.getLogger("shape_mapping")


gse_number = "GSE253903"
base_dir = Path(r"/usr/local/storage/data_scRNAseq")
#process_gse(gse_number, base_dir / "raw_data_download", base_dir / "raw_data")


gse_number = "GSE253903"
process_gse(gse_number, base_dir / "raw_data_download", base_dir / "raw_data", logger)
data_dir = base_dir / "raw_data" / gse_number
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
adata = umap(adata, logger, save_path_plot = r"/home/f/flor/test_umap.png")
adata = annotate_adata(adata, logger, "/home/f/flor/gitrepos/plaque-atlas-mapping/automatic_mapping_level1.py",
                       save_path_plot = r"/home/f/flor/test_umap_annot.png", save_path=f"/home/f/flor/{gse_number}_annot.h5ad")
print(adata)
print(run_limma_dea(adata, "ECs", ("TRUE", "FALSE"), "ECs_atherosclerosis.tsv", 
                    condition_col="symptomatic_atherosclerosis", logger=logger,
                     cell_type_col="cell_type_level1"))