# For miRNA datasets
source("code/microarray_functions.R")

dataset_accession <- c(
  "GSE111794",
  "GSE137582",
  "GSE137580",
  "GSE205119",
  "GSE130486",
  "GSE74755",
  "GSE45433",
  "GSE52243"
)

for (ds in dataset_accessions) download_microarray_data(ds, "/usr/local/storage/data_microarray/raw_data")

contrasts <- list(
  GSE111794 = list( 
    c("FALSE", "TRUE")
  ),
  GSE137582 = list( 
    c("WT", "LDLR(W483STOP)+/+")
  ),
  GSE137580 = list( 
    c("Ctrl", "oxLDL")
  ),
  GSE205119 = list(
    c("hsa_circ_0122319", "hsa_circ_0002457")
  ),
  GSE130486 = list( 
    c("Ctrl", "Pi", "covariates"),
    c("VSMCs|thoratic aorta||||Ctrl|WT||0|||||", "VSMCs|thoratic aorta||||Pi|WT||3|||||"),
    c("VSMCs|thoratic aorta||||Ctrl|WT||0|||||", "VSMCs|thoratic aorta||||Pi|WT||6|||||"),
    c("VSMCs|thoratic aorta||||Pi|WT||3|||||", "VSMCs|thoratic aorta||||Pi|WT||6|||||")
  ),
  GSE74755 = list( 
    c("VSMCs|thoratic aorta||||Pi|WT||3|||||", "VSMCs|thoratic aorta||||Pi|WT||6|||||")
  ),
  GSE45433 = list( 
    c("Ctrl", "injury")
  ),
  GSE52243 = list( 
    c("left carotid artery", "right carotid artery", "duration"),
    c("tissue|left carotid artery||||postPartialLigation|||12|||||", "tissue|left carotid artery||||postPartialLigation|||48|||||"),
    c("tissue|left carotid artery||||postPartialLigation|||12|||||", "tissue|right carotid artery||||postPartialLigation|||12|||||"),
    c("tissue|left carotid artery||||postPartialLigation|||12|||||", "tissue|right carotid artery||||postPartialLigation|||48|||||"),
    c("tissue|right carotid artery||||postPartialLigation|||12|||||", "tissue|left carotid artery||||postPartialLigation|||48|||||"),
    c("tissue|right carotid artery||||postPartialLigation|||12|||||", "tissue|right carotid artery||||postPartialLigation|||48|||||"),
    c("tissue|left carotid artery||||postPartialLigation|||48|||||", "tissue|right carotid artery||||postPartialLigation|||48|||||")
  )
)

dea_col <- list(
  GSE111794 = c("symptomatic_atherosclerosis"),
  GSE137582 = c("KO"),
  GSE137580 = c("medication"),
  GSE205119 = c("KO"),
  GSE130486 = c("medication", "combinded_condition", "combinded_condition", "combinded_condition"),
  GSE74755 = c("combined_condition"),
  GSE45433 = c("medication"),
  GSE52243 = c("location", "combinded_condition", "combinded_condition", "combinded_condition", "combinded_condition", "combinded_condition", "combinded_condition")
)


for (ds in dataset_accessions) {
  message("\n==============================")
  message("Working on ", ds)
  message("==============================")
  raw_data <- load_data(paste0("/usr/local/storage/data_microarray/raw_data/", ds), "/home/f/flor/metadata_all_samples.txt", sep = if (ds == "GSE72180") "_" else ".")
  data1 <- background_correction(raw_data)
  data1 <- annotate_data(data1)
  data1 <- clean_genes(data1)
  data1 <- collapse_duplicate_genes(
    data1,
    symbol_col = if (ds %in% c("GSE137580", "GSE205119", "GSE137582", "GSE130486", "GSE74755", "GSE45433")) "GeneName"  else "Symbol",
    method = "mean"
  )
  save_expression_matrix(data1, save.dir = "/usr/local/storage/data_microarray/260927_background_corrected", save.name = paste0(ds, "_background_corr_data.tsv"))
  data2 <- normalization(raw_data)
  data2 <- annotate_data(data2)
  data2 <- clean_genes(data2)
  data2 <- collapse_duplicate_genes(
    data2,
    symbol_col = if (ds %in% c("GSE137580", "GSE205119", "GSE137582", "GSE130486", "GSE74755", "GSE45433")) "GeneName"  else "Symbol",
    method = "mean"
  )
  save_expression_matrix(data2, save.dir = paste0("/usr/local/storage/data_microarray/260927_normalized"), save.name = paste0(ds, "_norm_data.tsv"))
  
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
      save.dir = "/usr/local/storage/data_microarray/260927_dea"
    )
  }
}