source("code/microarray_functions.R")

dataset_accessions <- c(
  "GSE152884",#
  "GSE72633",#
  "GSE42419",#
  "GSE152625",#
  "GSE28117",#
  "GSE111782",#
  "GSE111794",#
  "GSE137578",#
  "GSE137580",#
  "GSE137581",#
  "GSE137582",#
  "GSE205119",#
  "GSE205120",#
  "GSE48006",# contains three unknown samples!! GSMXXXn
  "GSE20739",#
  "GSE72180",#
  "GSE66624",#
  "GSE20060"#
)

for (ds in dataset_accessions) download_microarray_data(ds, "/usr/local/storage/data_microarray/raw_data")

contrasts = list(
    GSE53999 = list(
        c("||atherogenic|Ctrl|LDL -/||", "||atherogenic|JQ1|LDL -/||"),
        c("||atherogenic|Ctrl|LDL -/||", "||atherogenic|TNF|LDL -/||"),
        c("||atherogenic|Ctrl|LDL -/||", "||atherogenic|JQ1, TNF|LDL -/||"),
        c("||atherogenic|JQ1|LDL -/||", "||atherogenic|TNF|LDL -/||"),
        c("||atherogenic|JQ1|LDL -/||", "||atherogenic|JQ1, TNF|LDL -/||"),
        c("||atherogenic|TNF|LDL -/||", "||atherogenic|JQ1, TNF|LDL -/||")
    ),
    GSE152884 = list(
        c("||chow|Ctrl|WT||", "||chow|Ctrl|Mef2a, Mef2c, Mef2d||")
    ),
    GSE72633 = list(
        c("|||Ctrl|siCtrl||", "|||Ctrl|siNotch1||"),
        c("|||Ctrl|siCtrl||", "|||Ctrl|ox-PAPC||"),
        c("|||Ctrl|siNotch1||", "|||Ctrl|ox-PAPC||")
    ),
    GSE42419 = list(
        c("||chow|Ctrl|WT||", "||chow|Ctrl|PPARg||"),
        c("||chow|Ctrl|WT||", "||atherogenic|Ctrl|WT||"),
        c("||chow|Ctrl|WT||", "||atherogenic|Ctrl|PPARg||"),
        c("||chow|Ctrl|PPARg||", "||atherogenic|Ctrl|WT||"),
        c("||chow|Ctrl|PPARg||", "||atherogenic|Ctrl|PPARg||"),
        c("||atherogenic|Ctrl|PPARg||", "||atherogenic|Ctrl|PPARg||")
    ),
    GSE152625 = list(
        c("|||direct coculture monocytes (CD16+)|WT||", "|||transwell coculture monocytes (CD16+)|WT||")
    ),
    GSE28117 = list(
        c("|||Ctrl|WT||", "|||Il-4 (1h)|WT||"),
        c("|||Ctrl|WT||", "|||Il-4 (2h)|WT||"),
        c("|||Ctrl|WT||", "|||Il-4 (4h)|WT||"),
        c("|||Ctrl|WT||", "|||Il-4 (8h)|WT||"),
        c("|||Ctrl|WT||", "|||Il-4 (16h)|WT||"),
        c("|||Ctrl|WT||", "|||Ctrl|siCtrl||"),
        c("|||Ctrl|WT||", "|||Il-4 (4h)|siCtrl||"),
        c("|||Ctrl|WT||", "|||Ctrl|siSTAT6 (oligo1)||"),
        c("|||Ctrl|WT||", "|||Il-4 (4h)|siSTAT6 (oligo1)||"),
        c("|||Ctrl|WT||", "|||Ctrl|adenovirusCtrl||"),
        c("|||Ctrl|WT||", "|||Ctrl|adenovirusSTAT6||"),
        c("|||Il-4 (1h)|WT||", "|||Il-4 (2h)|WT||"),
        c("|||Il-4 (1h)|WT||", "|||Il-4 (4h)|WT||"),
        c("|||Il-4 (1h)|WT||", "|||Il-4 (8h)|WT||"),
        c("|||Il-4 (1h)|WT||", "|||Il-4 (16h)|WT||"),
        c("|||Il-4 (1h)|WT||", "|||Ctrl|siCtrl||"),
        c("|||Il-4 (1h)|WT||", "|||Il-4 (4h)|siCtrl||"),
        c("|||Il-4 (1h)|WT||", "|||Ctrl|siSTAT6 (oligo1)||"),
        c("|||Il-4 (1h)|WT||", "|||Il-4 (4h)|siSTAT6 (oligo1)||"),
        c("|||Il-4 (1h)|WT||", "|||Ctrl|adenovirusCtrl||"),
        c("|||Il-4 (1h)|WT||", "|||Ctrl|adenovirusSTAT6||"),
        c("|||Il-4 (2h)|WT||", "|||Il-4 (4h)|WT||"),
        c("|||Il-4 (2h)|WT||", "|||Il-4 (8h)|WT||"),
        c("|||Il-4 (2h)|WT||", "|||Il-4 (16h)|WT||"),
        c("|||Il-4 (2h)|WT||", "|||Ctrl|siCtrl||"),
        c("|||Il-4 (2h)|WT||", "|||Il-4 (4h)|siCtrl||"),
        c("|||Il-4 (2h)|WT||", "|||Ctrl|siSTAT6 (oligo1)||"),
        c("|||Il-4 (2h)|WT||", "|||Il-4 (4h)|siSTAT6 (oligo1)||"),
        c("|||Il-4 (2h)|WT||", "|||Ctrl|adenovirusCtrl||"),
        c("|||Il-4 (2h)|WT||", "|||Ctrl|adenovirusSTAT6||"),
        c("|||Il-4 (4h)|WT||", "|||Il-4 (8h)|WT||"),
        c("|||Il-4 (4h)|WT||", "|||Il-4 (16h)|WT||"),
        c("|||Il-4 (4h)|WT||", "|||Ctrl|siCtrl||"),
        c("|||Il-4 (4h)|WT||", "|||Il-4 (4h)|siCtrl||"),
        c("|||Il-4 (4h)|WT||", "|||Ctrl|siSTAT6 (oligo1)||"),
        c("|||Il-4 (4h)|WT||", "|||Il-4 (4h)|siSTAT6 (oligo1)||"),
        c("|||Il-4 (4h)|WT||", "|||Ctrl|adenovirusCtrl||"),
        c("|||Il-4 (4h)|WT||", "|||Ctrl|adenovirusSTAT6||"),
        c("|||Il-4 (8h)|WT||", "|||Il-4 (16h)|WT||"),
        c("|||Il-4 (8h)|WT||", "|||Ctrl|siCtrl||"),
        c("|||Il-4 (8h)|WT||", "|||Il-4 (4h)|siCtrl||"),
        c("|||Il-4 (8h)|WT||", "|||Ctrl|siSTAT6 (oligo1)||"),
        c("|||Il-4 (8h)|WT||", "|||Il-4 (4h)|siSTAT6 (oligo1)||"),
        c("|||Il-4 (8h)|WT||", "|||Ctrl|adenovirusCtrl||"),
        c("|||Il-4 (8h)|WT||", "|||Ctrl|adenovirusSTAT6||"),
        c("|||Il-4 (16h)|WT||", "|||Ctrl|siCtrl||"),
        c("|||Il-4 (16h)|WT||", "|||Il-4 (4h)|siCtrl||"),
        c("|||Il-4 (16h)|WT||", "|||Ctrl|siSTAT6 (oligo1)||"),
        c("|||Il-4 (16h)|WT||", "|||Il-4 (4h)|siSTAT6 (oligo1)||"),
        c("|||Il-4 (16h)|WT||", "|||Ctrl|adenovirusCtrl||"),
        c("|||Il-4 (16h)|WT||", "|||Ctrl|adenovirusSTAT6||"),
        c("|||Ctrl|siCtrl||", "|||Il-4 (4h)|siCtrl||"),
        c("|||Ctrl|siCtrl||", "|||Ctrl|siSTAT6 (oligo1)||"),
        c("|||Ctrl|siCtrl||", "|||Il-4 (4h)|siSTAT6 (oligo1)||"),
        c("|||Ctrl|siCtrl||", "|||Ctrl|adenovirusCtrl||"),
        c("|||Ctrl|siCtrl||", "|||Ctrl|adenovirusSTAT6||"),
        c("|||Il-4 (4h)|siCtrl||", "|||Ctrl|siSTAT6 (oligo1)||"),
        c("|||Il-4 (4h)|siCtrl||", "|||Il-4 (4h)|siSTAT6 (oligo1)||"),
        c("|||Il-4 (4h)|siCtrl||", "|||Ctrl|adenovirusCtrl||"),
        c("|||Il-4 (4h)|siCtrl||", "|||Ctrl|adenovirusSTAT6||"),
        c("|||Ctrl|siSTAT6 (oligo1)||", "|||Il-4 (4h)|siSTAT6 (oligo1)||"),
        c("|||Ctrl|siSTAT6 (oligo1)||", "|||Ctrl|adenovirusCtrl||"),
        c("|||Ctrl|siSTAT6 (oligo1)||", "|||Ctrl|adenovirusSTAT6||"),
        c("|||Il-4 (4h)|siSTAT6 (oligo1)||", "|||Ctrl|adenovirusCtrl||"),
        c("|||Il-4 (4h)|siSTAT6 (oligo1)||", "|||Ctrl|adenovirusSTAT6||"),
        c("|||Ctrl|adenovirusCtrl||", "|||Ctrl|adenovirusSTAT6||")
    ),
    GSE111782 = list(
        c("|||Ctrl|WT|TRUE|", "|||Ctrl|WT|FALSE|")
    ),
    GSE111794 = list(
        c("|||Ctrl|WT|TRUE|", "|||Ctrl|WT|FALSE|")
    ),
    GSE137582 = list(
        c("||HFD|Ctrl|WT||", "||HFD|Ctrl|LDLR(W483STOP) +/+||")
    ),
    GSE137581 = list(
        c("||HFD|Ctrl|WT||", "||HFD|Ctrl|LDLR(W483STOP) +/+||")
    ),
    GSE137580 = list(
        c("|||Ctrl|WT||", "|||oxLDL|WT||")
    ),
    GSE137578 = list(
        c("|||Ctrl|WT||", "|||oxLDL|WT||")
    ),
    GSE205119 = list(
        c("|||Ctrl|WT||", "|||Ctrl|hsa_circ_0122319||"),
        c("|||Ctrl|WT||", "|||Ctrl|hsa_circ_0002457||"),
        c("|||Ctrl|hsa_circ_0122319||", "|||Ctrl|hsa_circ_0002457||")
    ),
    GSE205120 = list(
        c("|||Ctrl|WT||", "|||Ctrl|hsa_circ_0122319||"),
        c("|||Ctrl|WT||", "|||Ctrl|hsa_circ_0002457||"),
        c("|||Ctrl|hsa_circ_0122319||", "|||Ctrl|hsa_circ_0002457||")
    ),
    GSE48006 = list(
        c("|||Ctrl|miR-21-3p||", "|||Ctrl|miR-27a-5p||"),
        c("|||Ctrl|miR-21-3p||", "|||Ctrl|miRCtrl||"),
        c("|||Ctrl|miR-27a-5p||", "|||Ctrl|miRCtrl||")
    ),
    GSE20739 = list(
        c("|||oscillatory wall shear stress|WT||", "|||oscillatory wall shear stress|miR-663-LNA||"),
        c("|||oscillatory wall shear stress|WT||", "|||stable laminar shear stress|WT||"),
        c("|||oscillatory wall shear stress|WT||", "|||stable laminar shear stress|miR-663-LNA||"),
        c("|||oscillatory wall shear stress|miR-663-LNA||", "|||stable laminar shear stress|WT||"),
        c("|||oscillatory wall shear stress|miR-663-LNA||", "|||stable laminar shear stress|miR-663-LNA||"),
        c("|||stable laminar shear stress|WT||", "|||stable laminar shear stress|miR-663-LNA||")
    ),
    GSE72180 = list(
        c("|||E2|KRR+||", "|||Ctrl|KRR+||"),
        c("|||E2|KRR+||", "|||E2|WT||"),
        c("|||E2|KRR+||", "|||Ctrl|WT||"),
        c("|||Ctrl|KRR+||", "|||E2|WT||"),
        c("|||Ctrl|KRR+||", "|||Ctrl|WT||"),
        c("|||E2|WT||", "|||Ctrl|WT||")
    ),
    GSE70126 = c(
        c("||HFD|Ctrl|Apoe -/-||", "||chow|Ctrl|WT||")
    ),
    GSE20060 = list(
        c("||chow|Ctrl|WT||", "||chow|oxPAPC|WT||")
    )
)

for (ds in dataset_accessions) {
  message("\n==============================")
  message("Working on ", ds)
  message("==============================")
  raw_data <- load_data(paste0("/usr/local/storage/data_microarray/raw_data/", ds), "/home/f/flor/metadata_all_samples.txt", sep = if (ds == "GSE72180") "_" else ".")
  data1 <- background_correction(raw_data, save.view = T, 
  save.dir = paste0("/usr/local/storage/data_microarray/background_corrected")
  )
  data2 <- normalization(raw_data, save.view = T, 
  save.dir = paste0("/usr/local/storage/data_microarray/normalized")
  )

  dataset_contrasts <- contrasts[[ds]]
  available_conditions <- unique(
    as.character(data2$pheno_data$combined_condition)
  )

  requested_conditions <- unique(unlist(dataset_contrasts, use.names = FALSE))
  absent_conditions <- setdiff(requested_conditions, available_conditions)
  if (length(absent_conditions) > 0L) {
    warning(
      "Skipping DEA for ", ds,
      ". The following requested combined_condition values are absent:\n  - ",
      paste(absent_conditions, collapse = "\n  - ")
    )
    next
  }
  
  for (comb in dataset_contrasts){
    contrast_str <- paste0(comb[2L], "-", comb[1L])
    dea <- run_dea(data2, "combined_condition", contrast_str = comb, save.view = T, 
    save.dir = paste0("/usr/local/storage/data_microarray/dea/", ds, "_", contrast_str, ".csv")
    )
  }
}