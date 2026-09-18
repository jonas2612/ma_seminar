import pandas as pd
from pathlib import Path
from tqdm import tqdm
import json
import logging

logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s [%(levelname)s]: %(message)s"
)
logger = logging.getLogger("shape_mapping")

def read_fuzzy_values(values_path):
    values_sep = [ ]
    for file in tqdm(values_path.glob("*.tsv")):
        tmp = pd.read_csv(file, sep="\t", index_col=0)
        tmp.index = ["_".join(".".join(file.name.split(".")[:-1]).split("_")[1:])]
        values_sep.append(tmp)
    return pd.concat(values_sep)

def determine_gene_list(value, level):
    if not isinstance(level, list):
        level = [level]
    maxes = value.max(axis=1)
    values_bool = value.copy()
    for col in values_bool.columns:
        values_bool[col] = (values_bool[col] == maxes)
    logger.info(values_bool.sum())

    mask = values_bool[level[0]]
    for l in level[1:]:
        mask |= values_bool[l]
    gene_list = values_bool[mask].index.tolist()
    return gene_list

base_dir = Path("/usr/local/storage/data_microarray/fuzzy_res")

logger.info("Determine D views with significant results")
rel_views = {}
for dir in base_dir.glob("*_D_padj"):
    miRNA_values = read_fuzzy_values(dir / "fuzzy_values")
    miRNA_padj_candidated = determine_gene_list(miRNA_values, ["****", "***", "**", "*"])
    if miRNA_padj_candidated:
        logger.info("{} contains significant results".format(dir.name[:-7]))
        rel_views[dir.name[:-7]] = {"padj": miRNA_padj_candidated}

for dir in base_dir.glob("*_D_logfc"):
    if dir.name[:-8] in rel_views.keys():
        miRNA_values = read_fuzzy_values(dir / "fuzzy_values")
        miRNA_logfc_candidated = determine_gene_list(miRNA_values, ["--", "-", "+", "++"])
        if miRNA_logfc_candidated:
            rel_views[dir.name[:-8]]["logfc"] = miRNA_logfc_candidated
            rel_views[dir.name[:-8]]["Overlap"] = set(miRNA_logfc_candidated).intersection(set(rel_views[dir.name[:-8]]['padj']))
            logger.info("{} contains results with higher logfc. Overlap with significant results: {}".format(dir.name[:-8], len(rel_views[dir.name[:-8]]["Overlap"])))

logger.info("Saving results in /home/f/flor/relevant_views.json")
with open('/home/f/flor/relevant_views.json', 'w') as fp:
    json.dump(rel_views, fp)