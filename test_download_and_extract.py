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
base_dir = Path(r"C:\Users\jonas\OneDrive\Desktop\Semester_4_Bio\Masterpraktikum\test_data")
#process_gse(gse_number, base_dir / "raw_data_download", base_dir / "raw_data")


gse_number = "GSE159677"
#process_gse(gse_number, base_dir / "raw_data_download", base_dir / "raw_data")
#data_dirs = determine_level(base_dir / "raw_data" / gse_number)
adatas = []
metadata = pd.read_csv(r"C:\Users\jonas\OneDrive\Desktop\Semester_4_Bio\Masterpraktikum\metadata_all_samples.txt", sep = "\t")
parent = Path(r"C:\Users\jonas\OneDrive\Desktop\Semester_4_Bio\Masterpraktikum\test_data\raw_data\GSE159677")
for path in tqdm([x for x in listdir(parent) if x.endswith("matrix.mtx.gz") and "GSM" in x]):
    logger.info(f"Working on {path}")
    adata = read_raw_data(parent / path, metadata, logger)
    adata = qc_statistical(adata, logger)
    #adata = doublet_detection(adata, logger)
    adatas.append(adata)
adata = combine_samples(adatas, logger)
adata = normalize_data(adata, logger)
adata = pca(adata, logger)
adata = neighbors(adata, logger)
adata = umap(adata, logger, save_path_plot = r"C:\Users\jonas\OneDrive\Desktop\Semester_4_Bio\Masterpraktikum\test_umap.png")
adata = annotate_adata(adata, logger, r"C:\Users\jonas\OneDrive\Desktop\gitrepos\plaque-atlas-mapping\automatic_mapping_level1.py",
                       save_path_plot = r"C:\Users\jonas\OneDrive\Desktop\Semester_4_Bio\Masterpraktikum\test_umap_annot.png")
print(adata)