from os import listdir

from code.scRNAseq_pipeline import *
import pandas as pd
import logging
from tqdm import tqdm

logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s [%(levelname)s]: %(message)s"
)
logger = logging.getLogger("shape_mapping")

base_dir = Path("/usr/local/storage/data_scRNAseq")
umap_dir = base_dir / "umaps"
dea_dir = base_dir / "dea"
annot_dir = base_dir / "annotated_data"

datasets_and_contrasts = {
   # "GSE309462": {"cond_col": "symptomatic_atherosclerosis", "contrast": [("False", "True")]},
   # "GSE260657": {"cond_col": "symptomatic_atherosclerosis", "contrast": [("False", "True")]}, # smartseq2 data
   # "GSE253903": {"cond_col": "symptomatic_atherosclerosis", "contrast": [("False", "True")]},
   # "GSE159677": {"cond_col": "cell_type", "contrast": [("plaque adjacent", "plaque")]},
   # "GSE260656": {"cond_col": "combined_condition", "contrast":[("|10w|||WT||", "|30w|||LDLR -/-, ApoB100/100||"), # smartseq2
   #                                                             ("|10w|||WT||", "|45w|||LDLR -/-, ApoB100/100||"),
   #                                                             ("|10w|||WT||", "|60w|||LDLR -/-, ApoB100/100||"),
   #                                                             ("|30w|||LDLR -/-, ApoB100/100||", "|20w|||LDLR -/-, ApoB100/100||"),
   #                                                             ("|30w|||LDLR -/-, ApoB100/100||", "|60w|||LDLR -/-, ApoB100/100||"),
   #                                                             ("|60w|||LDLR -/-, ApoB100/100||", "|20w|||LDLR -/-, ApoB100/100||"),
   #                                                             ("|60w|||LDLR -/-, ApoB100/100||", "|20w|||LDLR -/-, ApoB100/100||"),
   #                                                             ("|20w|||LDLR -/-, ApoB100/100||", "|20w|||LDLR -/-, ApoB100/100||")]}
   #"GSE131776": {"cond_col": "combined_condition", "contrast": [("||high-fat (16w)||Tcf21 flx/flx, Apoe -/-||", "||high-fat (8w)||Tcf21 flx/flx, Apoe -/-||"),
   #                                                             ("||high-fat (16w)||Apoe -/-||", "||high-fat (8w)||Apoe -/-||"),
   #                                                             ("||high-fat (16w)||Apoe -/-||", "||chow||Apoe -/-||"),
   #                                                             ("||high-fat (8w)||Apoe -/-||", "||chow||Apoe -/-||")]}
   #"GSE246083": {"cond_col": "combinded_condition", "contrast": [("||||WT||", "||||Atgl ECKO||")]}
   "GSE150644": {"cond_col": "combinded_condition", "contrast": [("||western diet (18w)||WT||", "||western diet (18w)||Kfl4 KO||")]}
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
    logger.debug("species of datasets:")
    logger.debug(", ".join([x.obs['species'].unique()[0] for x in adatas]))
    if key in ["GSE260656", "GSE150644", "GSE246083", "GSE131776"]:
        adata = combine_samples(adatas, logger, species="Mouse")
    else:
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
            try:
                run_limma_dea(adata, ct, contr, dea_dir / f"{key}_{ct}_{value['cond_col']}_{contr[1]}-{contr[0]}.tsv",
                      condition_col=value['cond_col'], logger=logger, cell_type_col="cell_type_level1")
            except RuntimeError as exc:
                error_message = str(exc)
                if "Insufficient independent pseudobulk profiles per group" in error_message:
                    message = (
                        f"Skipping pseudobulk DEA: dataset={key}; "
                        f"cell type={ct!r}; "
                        f"contrast={contr[1]!r} vs {contr[0]!r}; "
                        f"condition column={value['cond_col']!r}. "
                        f"Reason: insufficient independent pseudobulk profiles.\n"
                        f"{error_message}"
                    )

                    if logger is not None:
                        logger.warning(message)
                    else:
                        print(message)
                elif "Count matrix is empty." in error_message:
                    message = (
                        f"Skipping pseudobulk DEA: dataset={key}; "
                        f"cell type={ct!r}; "
                        f"contrast={contr[1]!r} vs {contr[0]!r}; "
                        f"condition column={value['cond_col']!r}. "
                        f"Reason: Count matrix is empty.\n"
                        f"{error_message}"
                    )
                    if logger is not None:
                        logger.warning(message)
                    else:
                        print(message)

                else:
                    raise
    