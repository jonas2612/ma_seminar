# For miRNA datasets
source("code/microarray_functions.R")

dataset_accessions <- c(
  "GSE105449",
  "GSE168149",
  "GSE89858"
)

for (ds in dataset_accessions) download_microarray_data(ds, "/usr/local/storage/data_microarray/raw_data")

  contrasts <- list(
    GSE105449=list(
      c("Ctrl", "after_medication", "CAD", "CVD", "age", "sex")
    ),
    GSE168149=list(
      c("TRUE", "FALSE", "symptomatic_atherosclerosis", "cell_type"),
      c("TRUE", "FALSE", "CAD", "cell_type")
    ),
    GSE89858=list(
      c("chow", "HFD", "duration")
    )
  )
  
  dea_col <- list(
    GSE105449=c("medication"),
    GSE168149=c("CAD", "symptomatic_atherosclerosis"),
    GSE89858=c("diet")
  )


for (ds in dataset_accessions) {
  message("\n==============================")
  message("Working on ", ds)
  message("==============================")
  raw_data <- load_data(paste0("/usr/local/storage/data_microarray/raw_data/", ds), "/home/f/flor/metadata_mirna_arrays_mewis.xlsx", sep = if (ds == "GSE72180") "_" else ".")
  data1 <- background_correction(raw_data)
  data1 <- annotate_data(data1)
  data1 <- clean_genes(data1)
  data1 <- collapse_duplicate_genes(
    data1,
    symbol_col = if (ds %in% c("GSE105449", "GSE89858")) "GeneName"  else "Symbol",
    method = "mean"
  )
  save_expression_matrix(data1, save.dir = "/usr/local/storage/data_microarray/261004_background_corrected_mewis", save.name = paste0(ds, "_background_corr_data.tsv"))
  data2 <- normalization(raw_data)
  data2 <- annotate_data(data2)
  data2 <- clean_genes(data2)
  data2 <- collapse_duplicate_genes(
    data2,
    symbol_col = if (ds %in% c("GSE105449", "GSE89858")) "GeneName"  else "Symbol",
    method = "mean"
  )
  save_expression_matrix(data2, save.dir = paste0("/usr/local/storage/data_microarray/261004_normalized_mewis"), save.name = paste0(ds, "_norm_data.tsv"))
  
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
      save.dir = "/usr/local/storage/data_microarray/261004_dea_mewis"
    )
  }
}