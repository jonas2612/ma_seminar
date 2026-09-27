# For miRNA datasets
source("code/microarray_functions.R")

dataset_accessions <- c(
  #'GSE10000', dataloading
  #'GSE101126',
  #'GSE111782',
  #'GSE113969', Identification of Probe names
  #'GSE12261',
  #'GSE13139',
  #'GSE135626',
  #'GSE137578',
  #'GSE137581',
  #'GSE143162', Gene name could not be determined
  #'GSE15062',
  #'GSE152625', too little replicated
  #'GSE152884',
  #'GSE15914',
  #'GSE17556',
  #'GSE20060',
  #'GSE205120',
  #'GSE20739',
  #'GSE24487',
  #'GSE26295',
  #'GSE28117', DEA not full rank
  #'GSE29903', Data loading
  #'GSE30004',
  #'GSE34262',
  #'GSE35676', Data loading
  #'GSE39264',
  #'GSE42419',
  #'GSE46097',
  'GSE48006',
  'GSE49519',
  'GSE50250',
  'GSE50251',
  'GSE66624',
  'GSE67182',
  'GSE70126',
  'GSE72180',
  'GSE72633'
)

for (ds in dataset_accessions) download_microarray_data(ds, "/usr/local/storage/data_microarray/raw_data")

contrasts <- list(
  GSE10000 = list(
    c("WT", "Apoe-/-"),
    c("WT", "Apoe-/-", "age", "location"),
    c("tissue|aorta||6|chow||WT|||||||", "tissue|aorta||6|chow||Apoe-/-|||||||"),
    c("tissue|aorta||32|chow||WT|||||||", "tissue|aorta||32|chow||Apoe-/-|||||||"),
    c("tissue|aorta||78|chow||WT|||||||", "tissue|aorta||78|chow||Apoe-/-|||||||"),
    c("tissue|aorta||78|chow||Apoe-/-|||||||", "tissue|advanced lesion||78|chow||Apoe-/-|||||||"),
    c("tissue|aorta||78|chow||Apoe-/-|||||||", "tissue|adventita||78|chow||Apoe-/-|||||||"),
    c("tissue|aorta||78|chow||Apoe-/-|||||||", "tissue|lesion||78|chow||Apoe-/-|||||||"),
    c("tissue|aorta||78|chow||Apoe-/-|||||||", "tissue|adventita||78|chow||WT|||||||"),
    c("tissue|aorta||78|chow||Apoe-/-|||||||", "tissue|media||78|chow||WT|||||||"),
    c("tissue|advanced lesion||78|chow||Apoe-/-|||||||", "tissue|adventita||78|chow||Apoe-/-|||||||"),
    c("tissue|advanced lesion||78|chow||Apoe-/-|||||||", "tissue|lesion||78|chow||Apoe-/-|||||||"),
    c("tissue|advanced lesion||78|chow||Apoe-/-|||||||", "tissue|adventita||78|chow||WT|||||||"),
    c("tissue|advanced lesion||78|chow||Apoe-/-|||||||", "tissue|media||78|chow||WT|||||||"),
    c("tissue|adventita||78|chow||Apoe-/-|||||||", "tissue|adventita||78|chow||WT|||||||"),
    c("tissue|adventita||78|chow||Apoe-/-|||||||", "tissue|media||78|chow||WT|||||||"),
    c("tissue|adventita||78|chow||WT|||||||", "tissue|media||78|chow||WT|||||||")
  ),
  GSE101126 = list(
    c("shCtrl", "shADK")
  ),
  GSE111782 = list(
    c("FALSE", "TRUE")
  ),
  GSE113969 = list(
    c("ECs|iPS|||||||||TRUE||FALSE|", "ECs|iPS|||||||||TRUE||TRUE|"),
    c("SMCS|iPS|||||||||TRUE||FALSE|", "SMCS|iPS|||||||||TRUE||TRUE|")
  ),
  GSE12261 = list(
    c("Ctrl", "2-Methoxyestradiol", "duration"),
    c("SMCs|aorta||||Ctrl|||4|||||", "SMCs|aorta||||2-Methoxyestradiol|||4|||||"),
    c("SMCs|aorta||||Ctrl|||30|||||", "SMCs|aorta||||2-Methoxyestradiol|||30|||||"),
    c("SMCs|aorta||||Ctrl|||4|||||", "SMCs|aorta||||Ctrl|||30|||||"),
    c("SMCs|aorta||||Ctrl|||4|||||", "SMCs|aorta||||2-Methoxyestradiol|||30|||||"),
    c("SMCs|aorta||||Ctrl|||30|||||", "SMCs|aorta||||2-Methoxyestradiol|||4|||||"),
    c("SMCs|aorta||||2-Methoxyestradiol|||4|||||", "SMCs|aorta||||2-Methoxyestradiol|||30|||||")
  ),
  GSE13139 = list(
    c("GFP_OE", "LOX1_OE", "duration", "medication"),
    c("HAECT|aorta||||oxLDL|GFP_OE||0|||||", "HAECT|aorta||||Ctrl|GFP_OE||6|||||"),
    c("HAECT|aorta||||oxLDL|GFP_OE||0|||||", "HAECT|aorta||||Ctrl|GFP_OE||12|||||"),
    c("HAECT|aorta||||oxLDL|GFP_OE||0|||||", "HAECT|aorta||||Ctrl|GFP_OE||24|||||")
  ),
  GSE135626 = list(
    c("SMCs|aorta|||chow||Apoe-/-|||||||", "SMCs|aorta|||HFD||Apoe-/-|||||||"),
    c("SMCs|aorta|||chow||Apoe-/-|||||||", "SMCs|aorta|||HFD||Apoe-/-,AAS|||||||"),
    c("SMCs|aorta|||HFD||Apoe-/-|||||||", "SMCs|aorta|||HFD||Apoe-/-,AAS|||||||")
  ),
  GSE137578 = list(
    c("Ctrl", "oxLDL")
  ),
  GSE137581 = list(
    c("WT", "LDLR(W483STOP)+/+")
  ),
  GSE143162 = list(
    c("WT", "LAMP2A"),
    c("Ctrl", "LDL")
  ),
  GSE15062 = list(
    c("Ctrl", "LTBR", "duration"),
    c("Ctrl", "TNF", "duration"),
    c("Ctrl", "TNF,LTBR", "duration"),
    c("LTBR", "TNF", "duration"),
    c("LTBR", "TNF,LTBR", "duration"),
    c("TNF", "TNF,LTBR", "duration")
  ),
  GSE152625 = list(
    c("direct coculture monocytes (CD16+)", "direct coculture monocytes (CD14+)")
  ),
  GSE152884 = list(
    c("WT", "Mef2a, Mef2c, Mef2d")
  ),
  GSE15914 = list(
    c("HDAd-0", "HDAd-gE3"),
    c("HDAd-0", "PBS"),
    c("HDAd-gE3", "PBS")
  ),
  GSE17556 = list(
    c("highGlucose", "highGlucose_TSP1"),
    c("highGlucose", "highGlucose_TSP1"),
    c("highGlucose", "lowGlucose"),
    c("highGlucose", "Manose"),
    c("highGlucose", "Manose_TSP1"),
    c("highGlucose", "TSP1")
  ),
  GSE20060 = list(
    c("Ctrl", "oxPAPC", "sample_nr")
  ),
  GSE205120 = list(
    c("hsa_circ_0122319", "hsa_circ_0002457")
  ),
  GSE20739 = list(
    c("WT", "miR-663-LNA"),
    c("oscillatory wall shear stress", "stable laminar shear stress"),
    c("EC|umbilical vein||||oscillatory wall shear stress|WT|||||||", "EC|umbilical vein||||oscillatory wall shear stress|miR-663-LNA|||||||"),
    c("EC|umbilical vein||||oscillatory wall shear stress|WT|||||||", "EC|umbilical vein||||stable laminar shear stress|WT|||||||"),
    c("EC|umbilical vein||||oscillatory wall shear stress|WT|||||||", "EC|umbilical vein||||stable laminar shear stress|miR-663-LNA|||||||"),
    c("EC|umbilical vein||||oscillatory wall shear stress|miR-663-LNA|||||||", "EC|umbilical vein||||stable laminar shear stress|WT|||||||"),
    c("EC|umbilical vein||||oscillatory wall shear stress|miR-663-LNA|||||||", "EC|umbilical vein||||stable laminar shear stress|miR-663-LNA|||||||"),
    c("EC|umbilical vein||||stable laminar shear stress|WT|||||||", "EC|umbilical vein||||stable laminar shear stress|miR-663-LNA|||||||")
  ),
  GSE26295 = list(
    c("medium", "3OC12Hsl")
  ),
  GSE28117 = list(
    c("WT", "siCtrl", "duration", "medication"),
    c("WT", "siSTAT6 (oligo1)", "duration", "medication")
  ),
  GSE29903 = list(
    c("PAPC", "oxPAPC", "dose", "sample_nr")
  ),
  GSE30004 = list(
    c("Ctrl", "TGFb")
  ),
  GSE34262 = list(
    c("Apoe-/-,miR-126+/+", "Apoe-/-,miR-126+-/-")
  ),
  GSE35676 = list(
    c("ECs||m|||||||||||", "Macrophages||m|||||||||||")
  ),
  GSE39264 = list(
    c("ECs|aorta|||||WT|||||||", "ECs|aorta|||||Apoe-/-|||||||"),
    c("ECs|aorta||||media|WT||4|||||", "ECs|aorta||||LPS|WT||4|||||")
  ),
  GSE42419 = list(
    c("EC|heart|||chow|Ctrl|WT|||||||", "EC|heart|||chow|Ctrl|PPARg|||||||"),
    c("EC|heart|||chow|Ctrl|WT|||||||", "EC|heart|||atherogenic|Ctrl|WT|||||||"),
    c("EC|heart|||chow|Ctrl|WT|||||||", "EC|heart|||atherogenic|Ctrl|PPARg|||||||"),
    c("EC|heart|||chow|Ctrl|PPARg|||||||", "EC|heart|||atherogenic|Ctrl|WT|||||||"),
    c("EC|heart|||chow|Ctrl|PPARg|||||||", "EC|heart|||atherogenic|Ctrl|PPARg|||||||"),
    c("EC|heart|||atherogenic|Ctrl|WT|||||||", "EC|heart|||atherogenic|Ctrl|PPARg|||||||")
  ),
  GSE46097 = list(
    c("FALSE", "TRUE", "duration", "sex", "age"),
    c("FALSE", "TRUE", "duration", "sex", "age")
  ),
  GSE48006 = list(
    c("miRCtrl", "miR-21-3p"),
    c("miRCtrl", "miR-27a-5p"),
    c("miR-21-3p", "miR-27a-5p"),
  ),
  GSE49519 = list(
    c("WT", "STAT1-/-", "medication"),
    c("media", "INFg", "KO")
  ),
  GSE50250 = list(
    c("ascending aorta", "coronary artery"),
    c("ascending aorta", "descending aorta"),
    c("ascending aorta", "mesenteric artery"),
    c("coronary artery", "descending aorta"),
    c("coronary artery", "mesenteric artery"),
    c("descending aorta", "mesenteric artery")
  ),
  GSE50251 = list(
    c("ascending aorta", "descending aorta")
  ),
  GSE66624 = list(
    c("WT", "V3+")
  ),
  GSE67182 = list(
    c("WT", "Smad4-/-")
  ),
  GSE72180 = list(
    c("KRR+", "WT", "medication"),
    c("Ctrl", "E2", "KO"),
    c("EC|||||E2|KRR+|||||||", "EC|||||Ctrl|KRR+|||||||"),
    c("EC|||||E2|KRR+|||||||", "EC|||||E2|WT|||||||"),
    c("EC|||||E2|KRR+|||||||", "EC|||||Ctrl|WT|||||||"),
    c("EC|||||Ctrl|KRR+|||||||", "EC|||||E2|WT|||||||"),
    c("EC|||||Ctrl|KRR+|||||||", "EC|||||Ctrl|WT|||||||"),
    c("EC|||||E2|WT|||||||", "EC|||||Ctrl|WT|||||||")
  ),
  GSE72633 = list(
    c("siCtrl", "siNotch1"),
    c("siCtrl", "oxPAPC"),
    c("siNotch1", "oxPAPC")
  )
)

dea_col <- list(
  GSE10000 = c("KO", "KO", "combinded_condition", "combinded_condition", "combinded_condition", "combinded_condition", "combinded_condition", "combinded_condition", "combinded_condition", "combinded_condition", "combinded_condition", "combinded_condition", "combinded_condition", "combinded_condition", "combinded_condition", "combinded_condition", "combinded_condition"),
  GSE101126 = c("KO"),
  GSE111782 = c("symptomatic_atherosclerosis"),
  GSE113969 = c("combinded_condition", "combined_condition"),
  GSE12261 = c("medication", "combined_condition", "combined_condition", "combined_condition", "combined_condition", "combined_condition", "combined_condition"),
  GSE13139 = c("KO", "combined_condition", "combined_condition", "combined_condition"),
  GSE135626 = c("combined_condition", "combined_condition", "combined_condition"),
  GSE137578 = c("medication"),
  GSE137581 = c("KO"),
  GSE143162 = c("KO", "medication"),
  GSE15062 = c("medication", "medication", "medication", "medication", "medication", "medication"),
  GSE152625 = c("medication"),
  GSE152884 = c("KO"),
  GSE15914 = c("medication", "medication", "medication"),
  GSE17556 = c("medication", "medication", "medication", "medication", "medication", "medication"),
  GSE20060 = c("medication"),
  GSE205120 = c("KO"),
  GSE20739 = c("KO", "medication", "combined_condition", "combined_condition", "combined_condition", "combined_condition", "combined_condition", "combined_condition"),
  GSE26295 = c("medication"),
  GSE28117 = c("KO", "KO"),
  GSE29903 = c("medication"),
  GSE30004 = c("medication"),
  GSE34262 = c("KO"),
  GSE35676 = c("combined_condition"),
  GSE39264 = c("combined_condition", "combined_condition"),
  GSE42419 = c("combined_condition", "combined_condition", "combined_condition", "combined_condition", "combined_condition", "combined_condition"),
  GSE46097 = c("diabetis", "CAD"),
  GSE48006 = c("KO", "KO", "KO", "KO", "KO", "KO", "KO", "KO", "KO", "KO"),
  GSE49519 = c("KO", "medication"),
  GSE50250 = c("location", "location", "location", "location", "location", "location"),
  GSE50251 = c("location"),
  GSE66624 = c("KO"),
  GSE67182 = c("KO"),
  GSE72180 = c("KO", "medication", "combined_condition", "combined_condition", "combined_condition", "combined_condition", "combined_condition", "combined_condition"),
  GSE72633 = c("KO", "KO", "KO")
)


for (ds in dataset_accessions) {
  message("\n==============================")
  message("Working on ", ds)
  message("==============================")
  raw_data <- load_data(paste0("/usr/local/storage/data_microarray/raw_data/", ds), "/home/f/flor/metadata_all_samples_new.xlsx", sep = if (ds == "GSE72180") "_" else ".")
  data1 <- background_correction(raw_data)
  data1 <- annotate_data(data1)
  data1 <- clean_genes(data1)
  data1 <- collapse_duplicate_genes(
    data1,
    symbol_col = if (ds %in% c("GSE137580", "GSE205119", "GSE137582", "GSE130486", "GSE74755", "GSE45433")) "GeneName"  else "Symbol",
    method = "mean"
  )
  save_expression_matrix(data1, save.dir = "/usr/local/storage/data_microarray/260927_background_corrected_mRNA", save.name = paste0(ds, "_background_corr_data.tsv"))
  data2 <- normalization(raw_data)
  data2 <- annotate_data(data2)
  data2 <- clean_genes(data2)
  data2 <- collapse_duplicate_genes(
    data2,
    symbol_col = if (ds %in% c("GSE137580", "GSE205119", "GSE137582", "GSE130486", "GSE74755", "GSE45433")) "GeneName"  else "Symbol",
    method = "mean"
  )
  save_expression_matrix(data2, save.dir = paste0("/usr/local/storage/data_microarray/260927_normalized_mRNA"), save.name = paste0(ds, "_norm_data.tsv"))
  
  dataset_contrasts <- contrasts[[ds]]
  dea_groups <- dea_col[[ds]]
  metadata <- if (inherits(data2, c("ExpressionSet", "ExpressionFeatureSet"))) {
    Biobase::pData(data2)
  } else if (inherits(data2, c("EList", "EListRaw", "uRNAList"))) {
    data2$targets
  } else {
    stop(
      "Could not extract metadata from normalized object for ",
      ds, ". Unsupported class: ",
      paste(class(data2), collapse = ", ")
    )
  }
  missing_group_columns <- setdiff(unique(dea_groups), colnames(metadata))
  if (length(missing_group_columns)) {
    stop(
      "Cannot run DEA for ", ds,
      ". Grouping-column(s) absent from metadata:\n  - ",
      paste(missing_group_columns, collapse = "\n  - ")
    )
  }
  for (i in seq_along(dataset_contrasts)) {
    comb <- dataset_contrasts[[i]]
    group_col <- dea_groups[[i]]
    
    if (length(comb) < 2L) {
      stop(
        "Invalid contrast ", i, " for ", ds,
        ": require c(reference_condition, test_condition, ...optional_covariates)."
      )
    }
    
    condition_pair <- trimws(as.character(comb[1:2]))
    
    if (anyNA(condition_pair) || any(!nzchar(condition_pair))) {
      stop(
        "Contrast ", i, " for ", ds,
        " has missing or empty reference/test condition labels."
      )
    }
    
    available_conditions <- trimws(as.character(metadata[[group_col]]))
    available_conditions <- unique(
      available_conditions[
        !is.na(available_conditions) & nzchar(available_conditions)
      ]
    )
    
    absent_conditions <- setdiff(condition_pair, available_conditions)
    
    if (length(absent_conditions)) {
      warning(
        "Skipping DEA contrast ", i, " for ", ds, ".\n",
        "Grouping column: `", group_col, "`.\n",
        "Missing requested condition(s):\n  - ",
        paste(absent_conditions, collapse = "\n  - "),
        "\nAvailable values:\n  - ",
        paste(sort(available_conditions), collapse = "\n  - ")
      )
      next
    }
    
    covariate_cols <- if (length(comb) > 2L) {
      trimws(as.character(comb[3:length(comb)]))
    } else {
      NULL
    }
    
    if (!is.null(covariate_cols)) {
      if (anyNA(covariate_cols) || any(!nzchar(covariate_cols))) {
        warning(
          "Skipping DEA contrast ", i, " for ", ds,
          ": covariate list contains missing or empty values."
        )
        next
      }
      
      missing_covariates <- setdiff(covariate_cols, colnames(metadata))
      
      if (length(missing_covariates)) {
        warning(
          "Skipping DEA contrast ", i, " for ", ds, ".\n",
          "Metadata does not contain requested covariate(s):\n  - ",
          paste(missing_covariates, collapse = "\n  - ")
        )
        next
      }
    }
    
    message(
      "DEA ", ds, ", contrast ", i, ": ",
      condition_pair[2L], " vs ", condition_pair[1L],
      " [group_col = ", group_col, "]",
      if (is.null(covariate_cols)) {
        " [unadjusted]"
      } else {
        paste0(" [covariates: ", paste(covariate_cols, collapse = ", "), "]")
      }
    )
    
    dea <- run_dea(
      data = data2,
      group_col = group_col,
      contrast_str = condition_pair,
      covariate_cols = covariate_cols,
      save.view = TRUE,
      save.dir = "/usr/local/storage/data_microarray/260927_dea_mRNA"
    )
  }
}