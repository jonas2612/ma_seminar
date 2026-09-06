library(tools)
library(stringr)
library(affy)
library(affyPLM)
library(oligo)
library(limma)
library(AgiMicroRna)
library(GEOquery)
library(ArrayExpress)
library(tibble)

read_metadata <- function(metadata_file_path) {
  if (!file.exists(metadata_file_path)) stop("Metadata file not found")
  meta <- read.delim(metadata_file_path, sep="\t", stringsAsFactors = F,
                     check.names = F, row.names = 4)
  meta
}

match_metadata_to_files <- function(files, meta) {
  samples <- rownames(meta)
  matched <- vapply(samples, function(x){
    hit <- files[stringr::str_detect(basename(files), fixed(x))]
    if (length(hit)==0) stop("No file found for sample: ", x)
    if (length(hit)>1) warning("Multiple files match ", x, ". Use the first file: ", hit[1])
    hit[1]
  }, character(1))
  matched
}

unpack_archives <- function(dataset_dir) {
  tar_files <- list.files(dataset_dir, pattern = "\\.tar$", full.names = TRUE, ignore.case = TRUE)
  for (f in tar_files) {
    utils::untar(f, exdir = dataset_dir)
  }
  
  gz_files <- list.files(dataset_dir, pattern = "\\.(cel|idat|txt|chp)\\.gz$", 
                         full.names = TRUE, ignore.case = TRUE)
  for (f in gz_files) {
    R.utils::gunzip(f, overwrite = FALSE, remove = FALSE)
  }
  
  zip_files <- list.files(dataset_dir, pattern = "\\.zip$", full.names = TRUE, ignore.case = TRUE)
  for (f in zip_files) {
    utils::unzip(f, exdir = dataset_dir)
  }
}

download_microarray_data <- function(accession, out_dir, force = F) {
  dataset_dir <- file.path(out_dir, accession)
  dir.create(dataset_dir, recursive = T, showWarnings = F)
  
  raw_files_exist <- function(path) {
    patterns <- c(
      "\\.cel(\\.gz)?$",
      "\\.idat(\\.gz)?$",
      "\\.txt(\\.gz)?$",
      "\\.RAW(\\.gz)?$",
      "\\.zip$",
      "\\.chp(\\.gz)?$"
    )
    any(vapply(patterns, function(p) {
      length(list.files(path, pattern = p, full.names = T, ignore.case = T)) > 0
    }, logical(1)))
  }
  
  if (!force && raw_files_exist(dataset_dir)) {
    warning("Raw files exist. Skipping download and unpacking.")
    return(invisible(dataset_dir))
  } 
  
  message("Download Data")
  if (grepl("^GSE", accession, ignore.case = T)) {
    GEOquery::getGEOSuppFiles(
      GEO = accession,
      makeDirectory = F,
      baseDir = dataset_dir
    )
  } else if (grepl("^E-MTAB-", accession, ignore.case = T)) {
    ArrayExpress::getAE(
      accession,
      path = dataset_dir,
      type = "raw"
    )
  } else {
    stop("Unsupported accession type.")
  }
  
  message("Unpack data")
  unpack_archives(dataset_dir)
  message("Downloaded ", accession, ".")
  
  return(invisible(dataset_dir))
}

unpack_geo_seq_archives <- function(
    dataset_dir,
    remove_archives = FALSE,
    force = FALSE,
    max_depth = 10L,
    quiet = FALSE
) {
  dataset_dir <- normalizePath(dataset_dir, mustWork = TRUE)

  is_raw_tar <- function(path) {
    grepl(
      "_RAW\\.tar(\\.gz|\\.bz2)?$",
      basename(path),
      ignore.case = TRUE,
      perl = TRUE
    )
  }

  is_archive <- function(path) {
    grepl(
      "\\.(tar|tar\\.gz|tgz|tar\\.bz2|tbz|tbz2|zip)$",
      basename(path),
      ignore.case = TRUE,
      perl = TRUE
    )
  }

  archive_stem <- function(path) {
    sub(
      "\\.(tar\\.gz|tar\\.bz2|tgz|tbz2?|tar|zip)$",
      "",
      basename(path),
      ignore.case = TRUE,
      perl = TRUE
    )
  }

  extract_archive <- function(archive, exdir) {
    dir.create(exdir, recursive = TRUE, showWarnings = FALSE)

    if (!quiet) {
      message(
        "Extracting ", basename(archive),
        " -> ", normalizePath(exdir)
      )
    }

    if (grepl("\\.zip$", archive, ignore.case = TRUE)) {
      utils::unzip(
        zipfile = archive,
        exdir = exdir,
        overwrite = force
      )
    } else {
      utils::untar(
        tarfile = archive,
        exdir = exdir,
      )
    }

    if (remove_archives) {
      unlink(archive, force = TRUE)
    }

    invisible(exdir)
  }

    is_gz_matrix_file <- function(path) {
    grepl(
      "\\.(mtx|tsv)\\.gz$",
      basename(path),
      ignore.case = TRUE,
      perl = TRUE
    )
  }

  decompress_gz_file <- function(gz_file) {
    out_file <- sub("\\.gz$", "", gz_file, ignore.case = TRUE)

    if (file.exists(out_file) && !force) {
      if (!quiet) {
        message("Skipping existing decompressed file: ", out_file)
      }
      return(invisible(out_file))
    }

    if (!quiet) {
      message(
        "Decompressing ", basename(gz_file),
        " -> ", basename(out_file)
      )
    }

    in_con <- gzfile(gz_file, open = "rb")
    out_con <- file(out_file, open = "wb")

    on.exit({
      close(in_con)
      close(out_con)
    }, add = TRUE)

    repeat {
      chunk <- readBin(in_con, what = "raw", n = 1024L * 1024L)

      if (!length(chunk)) {
        break
      }

      writeBin(chunk, out_con)
    }

    if (remove_archives) {
      unlink(gz_file, force = TRUE)
    }

    invisible(out_file)
  }

  decompress_gz_matrix_files <- function(path) {
    gz_files <- list.files(
      path,
      full.names = TRUE,
      recursive = TRUE,
      ignore.case = TRUE
    )

    gz_files <- gz_files[
      vapply(gz_files, is_gz_matrix_file, logical(1))
    ]

    if (!length(gz_files)) {
      return(character())
    }

    if (!quiet) {
      message("Decompressing ", length(gz_files), " .mtx.gz/.tsv.gz file(s).")
    }

    vapply(gz_files, decompress_gz_file, character(1))
  }

  extracted <- character()

  raw_archives <- list.files(
    dataset_dir,
    pattern = "_RAW\\.tar(\\.gz|\\.bz2)?$",
    full.names = TRUE,
    recursive = FALSE,
    ignore.case = TRUE
  )

  if (length(raw_archives) > 1L) {
    warning(
      "Found multiple top-level *_RAW archives. Extracting all of them.",
      call. = FALSE
    )
  }

  for (archive in raw_archives) {
    extract_archive(archive, exdir = dataset_dir)
    extracted <- c(extracted, archive)
  }

  for (depth in seq_len(max_depth)) {
    archives <- list.files(
      dataset_dir,
      full.names = TRUE,
      recursive = TRUE,
      ignore.case = TRUE
    )

    archives <- archives[
      vapply(archives, is_archive, logical(1)) &
        !vapply(archives, is_raw_tar, logical(1))
    ]

    archives <- setdiff(archives, extracted)

    if (length(archives) == 0L) {
      break
    }

    if (!quiet) {
      message(
        "Nested archive layer ", depth, ": ",
        length(archives), " archive(s)."
      )
    }

    for (archive in archives) {
      # GSM..._Sample.tar -> dataset_dir/GSM..._Sample/
      sample_dir <- file.path(dirname(archive), archive_stem(archive))
      extract_archive(archive, exdir = sample_dir)
      extracted <- c(extracted, archive)
    }
  }

  decompressed_files <- decompress_gz_matrix_files(dataset_dir)

  # 3. Report usable Matrix Market directories.
  is_10x_dir <- function(path) {
    files <- list.files(path, full.names = FALSE, ignore.case = TRUE)

    any(grepl("^matrix\\.mtx$", files, ignore.case = TRUE)) &&
      any(grepl("^barcodes\\.tsv$", files, ignore.case = TRUE)) &&
      any(grepl(
        "^(features|genes)\\.tsv$",
        files,
        ignore.case = TRUE
      ))
  }

  candidate_dirs <- unique(c(
    dataset_dir,
    list.dirs(dataset_dir, recursive = TRUE, full.names = TRUE)
  ))

  tenx_dirs <- candidate_dirs[
    vapply(candidate_dirs, is_10x_dir, logical(1))
  ]

  remaining_archives <- list.files(
    dataset_dir,
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
  remaining_archives <- remaining_archives[
    vapply(remaining_archives, is_archive, logical(1))
  ]

  if (!quiet) {
    message(
      "Extracted ", length(extracted), " archive(s); found ",
      length(tenx_dirs), " complete 10x directory/directories."
    )
  }

  invisible(list(
    dataset_dir = dataset_dir,
    extracted_archives = extracted,
    decompressed_files = decompressed_files,
    tenx_dirs = tenx_dirs,
    remaining_archives = remaining_archives
  ))
}

download_geo_seq_data <- function(accession, out_dir, force = F, unpack = F, pattern = NULL, quiet = F) {
  accession <- toupper(trimws(accession))

  dataset_dir <- file.path(out_dir, accession)
  dir.create(dataset_dir, recursive = T, showWarnings = F)

  has_seq_data <- function(path) {
    patterns <- c(
      "\\.(mtx|h5|h5ad|loom|rds|rda|rdata)(\\.gz)?$",
      "\\.(csv|tsv|txt)(\\.gz)?$",
      "(matrix|counts?|expression|metadata|barcodes?|features?|genes).*",
      "\\.(fastq|fq|vam|cram)(\\.gz)?$",
      "\\.(zip|tar|tgz|tar\\.gz|tar\\.bz|7z)$"
    )
    any(vapply(patterns, function(pattern) {
      length(list.files(
        path,
        pattern = pattern,
        full.names = T,
        recursive = T,
        ignore.case = T
      )) > 0L
    }, logical(1)))
  }

  if (!force && has_seq_data(dataset_dir)) {
    warning(sprintf("Sequencing-asso. files already in '%s'. Skip download.", dataset_dir), call. = F)
    return(invisible(dataset_dir))
  }

  if (!requireNamespace("GEOquery", quietly = T)) {
    stop("Install 'GEOquery' with BiocManager::install('GEOquery')", call. = F)
  }
  if (!quiet) message("Downloading GEO supp. files for ", accession, "...")
  
  GEOquery::getGEOSuppFiles(
    GEO = accession,
    makeDirectory = F,
    baseDir = dataset_dir,
    filter_regex = pattern
  )

  archive_pattern <- "\\.(zip|tar|tgz|tar\\.gz|tar\\.bz2|7z)$"

  if (unpack) {
    if (!quiet) message("Unpacking downloaded archives ...")
    unpack_geo_seq_archives(dataset_dir, remove_archives = T)
  }
  files <- list.files(
    dataset_dir, full.names = T, recursive = unpack
  )

  if (!quiet) message("Finished ", accession, ".")
  invisible(dataset_dir)
}

create_agilent_targets <- function(data_dir) {
  targets_path <- file.path(data_dir, "targets.txt")
  txt_files <- list.files(data_dir, pattern = "\\.txt$", full.names = F, ignore.case = T)
  if (length(txt_files) == 0) {
    stop("No Agilent .txt files were found. No targets.txt file could be generated")
  }
  
  targets <- data.frame(
    FileName = txt_files,
    stringsAsFactors = F
  )
  utils::write.table(
    targets, file = targets_path, sep = "\t", quote = F, row.names = F
  )
  invisible(targets_path)
}

make_annotation_df <- function(meta_df, file_paths) {
  tmp <- meta_df
  rownames(tmp) <- basename(file_paths)
  Biobase::AnnotatedDataFrame(data=tmp)
}

make_targets <- function(meta_df, file_paths, data_dir, col = "FileName") {
  treatment <- apply(meta_df[, c("age", "diet", "medication", "KO", "symptomatic_atherosclerosis"), drop=F], 1, function(row) {paste(row, collapse="|")})
  gerep <- as.integer(factor(treatment, levels=unique(treatment)))
  paths <- file.path(data_dir, basename(file_paths))
  df <- data.frame(setNames(list(paths), col),
                  Treatment = treatment,
                  GErep = gerep,
                   meta_df,
                   stringsAsFactors = F,
                   row.names = NULL)
  rownames(df) <- tools::file_path_sans_ext(df[[col]])
  df
}

read_chp_files <- function(files, meta) {
  if (!requireNamespace("affxparser", quietly = T)) {
    stop("package 'affxparser' required for .chp files.")
  }
  extract_chp_table <- function(chp, file) {
  candidates <- c(
    "QuantificationEntries",
    "ProbeSetResults",
    "ProbeSetName",
    "ProbeSets",
    "probeSets",
    "probesets",
    "results",
    "Results",
    "expression",
    "Expression",
    "quantification",
    "Quantification"
  )

  component_names <- names(chp)

  preferred <- match(
    tolower(candidates),
    tolower(component_names),
    nomatch = 0L
  )

  preferred <- preferred[preferred > 0L]

  if (!length(preferred)) {
    stop(
      "Could not find an expression-result component in CHP file: ",
      basename(file), "\n",
      "Available top-level components: ",
      paste(component_names, collapse = ", ")
    )
  }

  table_name <- component_names[preferred[1L]]
  result <- chp[[table_name]]

  if (is.data.frame(result) || is.matrix(result)) {
    tbl <- as.data.frame(
      result,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

  } else if (is.list(result)) {
    field_lengths <- vapply(result, length, integer(1))

    nonzero <- field_lengths[field_lengths > 0L]

    if (!length(nonzero)) {
      stop(
        "CHP component '", table_name,
        "' contains no non-empty fields in ", basename(file)
      )
    }

    expected_length <- max(nonzero)

    keep <- names(field_lengths)[field_lengths == expected_length]

    if (!length(keep)) {
      stop(
        "No equally sized vector fields found in CHP component '",
        table_name, "' for ", basename(file)
      )
    }

    tbl <- as.data.frame(
      result[keep],
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

  } else {
    stop(
      "Unsupported class for CHP component '", table_name,
      "' in ", basename(file), ": ",
      paste(class(result), collapse = "/")
    )
  }

  if (!nrow(tbl)) {
    stop(
      "Zero features extracted from CHP component '", table_name,
      "' in ", basename(file), ".\n",
      "Detected fields: ",
      paste(colnames(tbl), collapse = ", ")
    )
  }

  list(
    table = tbl,
    table_name = table_name
  )
}
  choose_column <- function(x, candidates, file, role) {
    idx <- match(
      tolower(candidates),
      tolower(colnames(x)),
      nomatch = 0L
    )

    idx <- idx[idx > 0L]

    if (!length(idx)) {
      stop(
        "Could not find a ", role,
        " field in CHP file: ", basename(file), "\n",
        "Available fields: ", paste(colnames(x), collapse = ", ")
      )
    }

    colnames(x)[idx[1L]]
  }

  one_file <- function(file) {
    chp <- affxparser::readChp(file)
    extracted <- extract_chp_table(chp, file)
    tab <- extracted$table

    id_col <- choose_column(
      tab,
      candidates = c(
        "ProbeSetName",
        "ProbeSetID",
        "ProbeSet",
        "Probe_ID",
        "probe_set_id",
        "ID",
        "Name"
      ),
      file = file,
      role = "feature-ID"
    )

    signal_col <- choose_column(
      tab,
      candidates = c(
        "QuantificationValue",
        "Quantification",
        "Signal",
        "Quant",
        "Expression",
        "Expression_Value",
        "Value",
        "Intensity"
      ),
      file = file,
      role = "expression-signal"
    )

    ids <- as.character(tab[[id_col]])
    signal <- suppressWarnings(as.numeric(tab[[signal_col]]))

    list(ids = ids, signal = signal, file = file, table_name = extracted$table_name, id_col = id_col, signal_col = signal_col)
  }

  parsed <- lapply(files, one_file)
  reference_ids <- parsed[[1L]]$ids

  for (i in seq_along(parsed)) {
    if (!setequal(reference_ids, parsed[[i]]$ids)) {
      stop(
        "Feature sets differ between CHP files:\n",
        "Reference: ", basename(files[1L]), "\n",
        "Different: ", basename(files[i]), "\n",
        "Mismatch: ", paste(setdiff(reference_ids, parsed[[i]]$ids)), "\n",
        "CHP-derived samples cannot be safely combined without ",
        "explicit feature reconciliation."
      )
    }
  }

  mat <- vapply(parsed, function(x) {
    x$signal[match(reference_ids, x$ids)]
  }, numeric(length(reference_ids)))
  if (is.null(dim(mat))) {
    mat <- matrix(mat, ncol=1L, dimnames=list(reference_ids, NULL))
  }

  rownames(mat) <- reference_ids
  colnames(mat) <- rownames(meta)
  
  if (anyNA(mat)) {
    warning(
      sum(is.na(mat)),
      " CHP signal values are missing after feature alignment."
    )
  }

  eset <- Biobase::ExpressionSet(
    assayData = mat,
    phenoData = Biobase::AnnotatedDataFrame(meta),
    featureData = Biobase::AnnotatedDataFrame(
      data.frame(
        ProbeName = reference_ids,
        row.names = reference_ids,
        stringsAsFactors = FALSE
      )
    )
  )

  attr(eset, "preprocessing_done") <- TRUE
  attr(eset, "preprocessing_status") <- "vendor_processed_CHP"
  attr(eset, "available_views") <- "normalized_only"
  attr(eset, "chp_fields") <- list(
    result_table = parsed[[1L]]$table_name,
    feature_id = parsed[[1L]]$id_col,
    signal = parsed[[1L]]$signal_col
  )

  eset
}

get_supplier <- function(data) {
  raw_supplier <- if (inherits(data, c("AffyBatch", "ExpressionSet", "ExpressionFeatureSet"))) {
    Biobase::pData(data)[["chip"]][1]
  } else {
    data$targets[["chip"]][1]   # EList, EListRaw, RGList, uRNAList
  }
  # Strip replicate suffixes added in metadata.tsv (e.g. "Array 1.3" → "Array 1.0")
  sub("\\.[0-9]+$", "", raw_supplier)
}

load_data <- function(data_dir, metadata_file_path, sep = ".", other.columns = "Detection") {
  if (!dir.exists(data_dir)) stop(paste0(data_dir, " does not exists"))
  accession <- basename(data_dir)
  metadata <- read_metadata(metadata_file_path)
  meta <- metadata[metadata[['accession']]==accession,]
  supplier <- meta[['chip']][1]
  if (supplier %in% c(
    "GeneChip PrimeView Human Gene Expression Array",
    "Affymetrix Human genome U133 Plus 2.0 Array",
    "Affymetrix Human genome U133A 2.0 Array",
    "Affymetrix HT Human Genome U133A",
    "Affymetrix GeneChip Human Genome U133 Plus 2.0 Array",
    "Affymetrix Human Genome U133 Plot 2.0 Array",
    "Affymetrix Rat Genome 230A"
  )){
    cel_files <- list.celfiles(data_dir, full.names = T)
    matched_files <- match_metadata_to_files(cel_files, meta)
    annot_df <- make_annotation_df(meta, matched_files)
    raw_data <- affy::ReadAffy(filenames = matched_files, phenoData = annot_df)
  }
  else if (supplier %in% c(
    "Affymetrix Mouse Transcriptome Array 1.0",
    "Affymetrix Human Gene 2.0 ST",
    "Affymetrix GeneChip Mouse Gene 1.0 ST Array",
    "Affymetrix Mouse Exon-Junction Array",
    "Affymetrix Rat Gene 1.0 ST",
    "Affymetrix GeneChip miRNA 3.0"
  )) {
    chp_files <- list.files(data_dir, pattern = "\\.chp$", full.names = T, ignore.case = T)
    if (length(chp_files) > 0) {
      message("reading .chp files")
      matched_files <- match_metadata_to_files(chp_files, meta)
      raw_data <- read_chp_files(files = matched_files, meta = meta)
      return(raw_data)
    }
    cel_files <- list.celfiles(data_dir, full.names=T)
    matched_files <- match_metadata_to_files(cel_files, meta)
    annot_df <- make_annotation_df(meta, matched_files)
    raw_data <- oligo::read.celfiles(filenames = matched_files, phenoData = annot_df) 
  }
  else if (supplier %in% c(
    "Illumina HumanHT-12 v4.0",
    "Illumina HumanHT-12 v3.0",
    "Illumina MouseRef-8 v2.0"
  )) {
    # add distinction for existance of .bgx files (manifest)
    idat_files <- list.files(data_dir, pattern = "\\.idat(\\.gz)?$", full.names = T)
    txt_files <- list.files(data_dir, pattern = "\\.txt(\\.gz)?$", full.names = T)
    bgx_files <- list.files(data_dir, pattern = "\\.bgx(\\.gz)?$", full.names = T)
    matrix_files <- list.files(data_dir, pattern = "\\.txt(\\.gz)?$", full.names = T, ignore.case = T)
    matrix_files <- matrix_files[stringr::str_detect(basename(matrix_files),
                                                     stringr::regex("non.?normalized|non.?norm", ignore_case = T))]
    gpl_files <- list.files(data_dir, pattern = "\\.txt(\\.gz)?$",
                          full.names = T, ignore.case = T)
    gpl_files <- gpl_files[stringr::str_detect(basename(gpl_files),
                                           stringr::regex("^GPL", ignore_case = T))]
    
    if (length(idat_files) != 0) {
      message("Illumina data is presented as .idat data")
      if (length(bgx_files)!=1) stop("No or more than 1 manifest files found")
      matched_files <- match_metadata_to_files(idat_files, meta)
      raw_data <- limma::read.idat(idatfiles = matched_files, bgxfile = bgx_files[1])
      id_candidates = c(
        "Probe_Id",
        "ProbeID",
        "ProbeName"
      )
      id_col <- intersect(id_candidates, colnames(raw_data$genes))
      ids <- trimws(as.character(raw_data$genes[[id_col[1L]]]))
      rownames(raw_data$E) <- ids
      raw_data$targets <- make_targets(meta, matched_files, data_dir)
    }
    else if (length(matrix_files) != 0) {
      message("Illumina data is presented as non-normalized matrix file")
      raw_data <- read_ilmn_matrix_safe(matrix_files[1], sep = sep, other.columns = other.columns)
      rownames(meta) <- meta[["Original Name"]]
      gsm_ids <- rownames(meta)
      message(gsm_ids)
      message(colnames(raw_data$E))
      col_match <- vapply(gsm_ids, function(gsm){
        hit <- grep(gsm, colnames(raw_data$E), value = T)
        if (length(hit) == 0) stop("sample column not found in matrix for ", gsm)
        if (length(hit) > 1) warning("multiple columns matched ", gsm, ". Using first hit")
        hit[1]
      }, character(1))
      raw_data$E <- raw_data$E[, col_match, drop = F]
      colnames(raw_data$E) <- gsm_ids
      id_candidates = c(
        "ProbeName",
        "PROBE_ID",
        "ProbeID",
        "ID_REF"
      )
      id_col <- intersect(id_candidates, colnames(raw_data$genes))
      ids <- trimws(as.character(raw_data$genes[[id_col[1L]]]))
      rownames(raw_data$E) <- ids
      if (!is.null(raw_data$other)) { #necessary to ensure matching signal/detection values
        raw_data$other <- lapply(
          raw_data$other, function(m) {
            m <- m[, col_match, drop=F]
            colnames(m) <- gsm_ids
            rownames(m) <- rownames(raw_data$E)
            m
          }
        )
      }
      raw_data$targets <- meta
    }
    else {
      message("Illumina data is presented as .txt data")
      matched_files <- match_metadata_to_files(txt_files, meta)
      raw_data <- limma::read.ilmn(files = txt_files,
                                   probeid = "ProbeID", 
                                   expr = "AVG_Signal", 
                                   other.columns = "Detection"
      )
      id_candidates = c(
        "ProbeName",
        "PROBE_ID",
        "ProbeID",
        "ID_REF"
      )
      id_col <- intersect(id_candidates, colnames(raw_data$genes))
      ids <- trimws(as.character(raw_data$genes[[id_col[1L]]]))
      rownames(raw_data$E) <- ids
      raw_data$targets <- make_targets(meta, matched_files, data_dir)
    }
  }
  else if (supplier %in% c(
    "Agilent-084555 026652QM_RCUG_HomoSapiens",
    "Agilent-074809 SurePrint G3 Mouse GE v2",
    "Agilent-072363 SurePrint G3 Human GE v3"
  )) {
    txt_files <- list.files(data_dir, pattern = "\\.txt$", full.names = T, ignore.case = T)
    matched_files <- match_metadata_to_files(txt_files, meta)
    targets <- make_targets(meta, matched_files, data_dir)
    raw_data <- read.maimages(targets, source="agilent", green.only=T)
    rownames(raw_data$E) <- raw_data$genes$ProbeName
  }
  else if (supplier %in% c(
    "Agilent miRNA microarray",
    "Agilent-070155 Mouse miRNA Microarray",
    "Agilent-070156 Human_miRNA_V21.0_Microarray",
    "Agilent-070156 Human miRNA",
    "Agilent-050340 Custom Rat miRNA Microarray",
    "Agilent-019159 Rat miRNA"
  )) {
    txt_files <- list.files(data_dir, pattern = "\\.txt$", full.names = T, ignore.case = T)
    matched_files <- match_metadata_to_files(txt_files, meta)
    targets <- make_targets(meta, matched_files, data_dir)
    targets_path <- file.path(data_dir, "targets.txt")
    utils::write.table(targets, file = targets_path, sep = "\t", quote = F, row.names = F)
    targets <- AgiMicroRna::readTargets(targets_path)
    raw_data <- AgiMicroRna::readMicroRnaAFE(targets = targets, verbose = F)
    sample_ids <- colnames(raw_data$TGS)
    raw_data$targets <- targets[sample_ids, , drop = FALSE]
    rownames(raw_data$TGS) <- raw_data$genes$ProbeName
    if (!is.null(raw_data$E)) {
    rownames(raw_data$E) <- raw_data$genes$ProbeName
  }
    #raw_data <- read_agilent_mirna_safe(targets = targets, verbose = F)
  }
  else {
    message("Please implement the microarray version within this function. Also, check the spelling of 'Microarray supplier'. Microarray supplier: ", supplier)
    raw_data <- NULL
  }
  return(raw_data)
}

read_ilmn_matrix_safe <- function(
  file_path,
  probe_id = NULL,
  expr = "AVG_Signal",
  sep = "_",
  other.columns = "Detection"
) {
  stopifnot(file.exists(file_path))

  lines <- readLines(file_path, warn = FALSE)

  probe_candidates <- c("ID_REF", "PROBE_ID", "ProbeID", "Probe")
  comment_line <- stringr::str_detect(lines, "^\\s*[\"']?[!#]")

  header_line_idx <- NA_integer_

  for (candidate in probe_candidates) {
    hit <- which(
      !comment_line &
        stringr::str_detect(
          lines,
          stringr::fixed(candidate)
        )
    )[1]

    if (!is.na(hit)) {
      header_line_idx <- hit
      if (is.null(probe_id)) probe_id <- candidate
      break
    }
  }

  if (is.na(header_line_idx)) {
    stop(
      "Could not find a header row containing one of: ",
      paste(probe_candidates, collapse = ", ")
    )
  }

  cleaned <- lines[header_line_idx:length(lines)]
  cleaned <- stringr::str_trim(cleaned, side = "right")
  cleaned <- cleaned[nzchar(cleaned)]

  footer_idx <- which(
    stringr::str_detect(
      cleaned[-1],
      "^\\s*[\"']?[!#]"
    )
  )

  if (length(footer_idx) > 0L) {
    cleaned <- cleaned[seq_len(footer_idx[1L])]
  }

  header_names <- strsplit(cleaned[1L], "\t", fixed = TRUE)[[1L]]
  header_names <- trimws(header_names)

  if (!probe_id %in% header_names) {
    stop(
      "Detected probe-ID field '", probe_id,
      "' is absent from parsed header: ",
      paste(header_names, collapse = ", ")
    )
  }

  sep_esc <- gsub(
    "([][{}()+*^$.|?\\\\])",
    "\\\\\\1",
    sep
  )

  expr_esc <- gsub(
    "([][{}()+*^$.|?\\\\])",
    "\\\\\\1",
    expr
  )

  other_esc <- gsub(
    "([][{}()+*^$.|?\\\\])",
    "\\\\\\1",
    other.columns
  )

  probe_col <- header_names == probe_id

  suffix_sep_pattern <- paste0(
    "(?:",
    sep_esc,
    "|\\.)"
  )

  signal_pattern <- paste0(
    suffix_sep_pattern,
    expr_esc,
    "$"
  )

  detection_pattern <- paste0(
    suffix_sep_pattern,
    "(?:",
    other_esc,
    ")",
    "(?:",
    "[ _.]?",
    "P(?:VAL|VALUE)",
    ")?$"
  )

  signal_col <- grepl(
    signal_pattern,
    header_names,
    ignore.case = TRUE,
    perl = TRUE
  )

  detection_col <- grepl(
    detection_pattern,
    header_names,
    ignore.case = TRUE,
    perl = TRUE
  )

  suffix_present <- any(signal_col)

  if (!suffix_present) {
    message(
      "No '", expr,
      "' suffix found: interpreting as a GEO-style sample matrix."
    )

    signal_col <- !probe_col & !detection_col
  }

  annotation_col <- !probe_col & !signal_col & !detection_col

  if (!any(signal_col)) {
    stop(
      "No signal columns detected. Header fields are:\n",
      paste(header_names, collapse = "\n")
    )
  }

  message(
    "Probe ID column: ", probe_id, "\n",
    "Signal columns: ", paste(header_names[signal_col], collapse = ", "), "\n",
    "Detection columns: ",
    if (any(detection_col)) paste(header_names[detection_col], collapse = ", ") else "<none>", "\n",
    "Annotation columns: ",
    if (any(annotation_col)) paste(header_names[annotation_col], collapse = ", ") else "<none>"
  )

  tmp <- tempfile(fileext = ".txt")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(cleaned, tmp, useBytes = TRUE)

  full <- utils::read.delim(
    tmp,
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = "",
    fill = TRUE
  )

  if (!identical(colnames(full), header_names)) {
    stop(
      "Parsed column names differ from the detected header. ",
      "This usually indicates quoted tabs or malformed rows."
    )
  }

  expr_mat <- as.matrix(full[, signal_col, drop = FALSE])
  storage.mode(expr_mat) <- "double"

  signal_names <- colnames(expr_mat)

  if (suffix_present) {
    sample_ids <- sub(
      paste0(suffix_sep_pattern, expr_esc, "$"),
      "",
      signal_names,
      ignore.case = TRUE,
      perl = TRUE
    )
  } else {
    sample_ids <- signal_names
  }

  if (anyDuplicated(sample_ids)) {
    stop(
      "Signal-column parsing generated duplicate sample IDs: ",
      paste(unique(sample_ids[duplicated(sample_ids)]), collapse = ", ")
    )
  }

  colnames(expr_mat) <- sample_ids

  other_list <- list()

  if (any(detection_col)) {
    detection_mat <- as.matrix(full[, detection_col, drop = FALSE])
    storage.mode(detection_mat) <- "double"

    detection_names <- colnames(detection_mat)

    detection_sample_ids <- sub(
      paste0(
        suffix_sep_pattern,
        "(?:",
        other_esc,
        ")",
        "(?:",
        "[ _.]?",
        "P(?:VAL|VALUE)",
        ")?$"
      ),
      "",
      detection_names,
      ignore.case = TRUE,
      perl = TRUE
    )

    if (setequal(detection_sample_ids, sample_ids)) {
      detection_mat <- detection_mat[
        ,
        match(sample_ids, detection_sample_ids),
        drop = FALSE
      ]
      colnames(detection_mat) <- sample_ids
    } else {
      warning(
        "Detection-p-value columns cannot be aligned one-to-one with ",
        "signal columns. Detection data will be retained but not renamed."
      )
    }

    other_list[["Detection PValue"]] <- detection_mat
  }

  genes <- full[, probe_col | annotation_col, drop = FALSE]
  colnames(genes)[colnames(genes) == probe_id] <- "ProbeName"

  if (anyDuplicated(genes$ProbeName)) {
    warning(
      "Duplicate probe identifiers found: ",
      sum(duplicated(genes$ProbeName)),
      ". Keep probe-level rows for now; resolve duplicates during annotation."
    )
  }

  annotation_numeric_fraction <- function(x) {
    x <- trimws(as.character(x))

    nonempty <- !is.na(x) & nzchar(x)

    if (!any(nonempty)) {
      return(0)
    }

    mean(
      !is.na(
        suppressWarnings(as.numeric(x[nonempty]))
      )
    )
  }

  annotation_fields <- intersect(
    c(
      "SYMBOL",
      "SEARCH_KEY",
      "TRANSCRIPT",
      "ILMN_GENE",
      "SOURCE_REFERENCE_ID",
      "REFSEQ_ID"
    ),
    colnames(full)
  )

  numeric_fraction <- vapply(
    full[, annotation_fields, drop = FALSE],
    annotation_numeric_fraction,
    numeric(1)
  )

  if (any(numeric_fraction > 0.80)) {
    stop(
      "Illumina annotation columns appear structurally shifted after parsing:\n",
      paste(
        names(numeric_fraction)[numeric_fraction > 0.80],
        sprintf("%.1f%% numeric", 100 * numeric_fraction[numeric_fraction > 0.80]),
        sep = ": ",
        collapse = "\n"
      ),
      "\nDo not continue with this parsed matrix."
    )
  }

  structure(
    list(
      E = expr_mat,
      genes = genes,
      other = other_list,
      targets = data.frame(
        sampleID = sample_ids,
        row.names = sample_ids,
        stringsAsFactors = FALSE
      )
    ),
    class = "EListRaw"
  )
}

detection_pval <- function(raw_data) {
  if (is.null(raw_data$other) || !length(raw_data$other)) {
    return(NULL)
  }
  names <- names(raw_data$other)

  preferred <- c(
    "Detection PValue",
    "Detection PVal",
    "Detection",
    "DetectionPValue",
    "Detection_PValue",
    "gDetectionPValue",
    "PValue"
  )

  idx <- match(tolower(preferred), tolower(names), nomatch=0L)
  idx <- idx[idx>0L]
  if (!length(idx)) {
    idx <- grep(
      "detect.*p.?val|p.?val.*detect|detection",
      names,
      ignore.case = TRUE,
      perl = TRUE
    )
  }
  if (!length(idx)) {
    return(NULL)
  }

  detectionName <- names[idx[1L]]
  detectionName
}

background_correction <- function(raw_data, annotate = T, clean.genes = T, condition_col = "combined_condition", save.view = F, save.dir = NULL) {
  supplier <- get_supplier(raw_data)
  if (isTRUE(attr(raw_data, "preprocessing_done"))) {
    message("Skipping background correction, as data is already preprocessed (.chp)")
    return (annotate_data(raw_data))
  }
  
  if (inherits(raw_data, c("EListRaw", "EList")) &&
      supplier %in% c("Illumina HumanHT-12 v4.0",
                      "Illumina HumanHT-12 v4",
                      "Illumina HumanHT-12 v3.0",
                      "Illumina HumanHT-12 v3",
                      "Illumina MouseRef-8 v2.0",
                      "Illumina MouseRef-8 v2")) {
    detectionName <- detection_pval(raw_data)
    if (!is.null(detectionName)) {
      raw_data_filtered <- filter_illumina(raw_data, detectionName)
      background_corrected <- limma::nec(raw_data_filtered, detection.p = detectionName)
    } else {
      background_corrected <- limma::nec(raw_data)
    }
  }
  else if (inherits(raw_data,"uRNAList")) {
    background_corrected <- raw_data #read function for Agilent MircoRNA Chips already did background correction
    background_corrected <- urna_to_elist(background_corrected)
  }
  else if (inherits(raw_data, c("RGList","EListRaw","EList")) &&
           supplier %in% c("Agilent-084555 026652QM_RCUG_HomoSapiens",
                           "Agilent-074809 SurePrint G3 Mouse GE v2",
                           "Agilent-072363 SurePrint G3 Human GE v3")) {
    background_corrected <- limma::backgroundCorrect(raw_data, method = "normexp", offset=50)
  }
  else if (inherits(raw_data,"AffyBatch")) {
    background_corrected <- affyPLM::preprocess(raw_data, normalize = F, background = T, background.method = "RMA.2")
    background_corrected <- affy::rma(background_corrected, background = F, normalize = F)
  }
  else if (inherits(raw_data, "ExpressionFeatureSet")) {
    background_corrected <- oligo::rma(raw_data, normalize = F, background = T)
  }
  else {
    stop("Background correction not implemented for class ", paste(class(raw_data), collapse=", "))
  }
  if (annotate)
    background_corrected <- annotate_data(background_corrected)
    if (clean.genes) background_corrected <- clean_genes(background_corrected)
    if (save.view) {
      if (is.null(save.dir)) stop("'save.dir' needs to be set.")
      dir.create(save.dir, recursive = TRUE, showWarnings = FALSE)
      if (inherits(background_corrected, c("ExpressionSet", "ExpressionFeatureSet"))) {
        E <- Biobase::exprs(background_corrected)
        metadata <- Biobase::pData(background_corrected)
      } else {
        E <- background_corrected$E
        metadata <- background_corrected$targets
      }
      accession <- unique(trimws(as.character(metadata[["accession"]])))
      accession <- accession[!is.na(accession) & nzchar(accession)]
      conditions <- unique(metadata$combined_condition)
      conditions <- conditions[!is.na(conditions) & nzchar(conditions)]

      for (condition in conditions) {
        sample_indices <- which(metadata$combined_condition == condition)

        E_condition <- E[, sample_indices, drop = FALSE]

        safe_condition <- gsub(
          "[^[:alnum:]_.-]+",
          "_",
          condition
        )
        safe_condition <- gsub("_+", "_", safe_condition)
        safe_condition <- gsub("^_|_$", "", safe_condition)

        output_file <- file.path(
          save.dir,
          paste0(accession, "_", safe_condition, "_background_corr_data.csv")
        )

        utils::write.csv(
          E_condition,
          file = output_file,
          row.names = TRUE,
          quote = FALSE
        )
      }
    }
  background_corrected
}

normalization <- function(raw_data, annotate = T, clean.genes = T, condition_col = "combined_condition", save.view = F, save.dir = NULL) {
  if (isTRUE(attr(raw_data, "preprocessing_done"))) {
    message("Skipping background correction, as data is already preprocessed (.chp)")
    return (annotate_data(raw_data))
  }
  supplier <- get_supplier(raw_data)
  if (inherits(raw_data, c("EListRaw","EList")) &&
      supplier %in% c("Illumina HumanHT-12 v4.0",
                      "Illumina HumanHT-12 v4",
                      "Illumina HumanHT-12 v3.0",
                      "Illumina HumanHT-12 v3",
                      "Illumina MouseRef-8 v2.0",
                      "Illumina MouseRef-8 v2")) {
    detectionName <- detection_pval(raw_data)
    if (!is.null(detectionName)) {
      raw_data_filtered <- filter_illumina(raw_data, detectionName)
      data <- limma::neqc(raw_data_filtered, detection.p = detectionName)
    } else {
      data <- limma::neqc(raw_data)
    }
  }
  else if (inherits(raw_data,"uRNAList")) {
    targets_path <- file.path(raw_data$targets$FileName[1] |> dirname(), "targets.txt")
    data <- AgiMicroRna::tgsMicroRna(raw_data, half = T, makePLOT = F)
    data <- AgiMicroRna::tgsNormalization(data,
                             targets = AgiMicroRna::readTargets(targets_path),
                             makePLOTpre = F, makePLOTpost = F)
    data <- urna_to_elist(data)
  }
  else if (inherits(raw_data, c("RGList","EListRaw","EList")) &&
           supplier %in% c("Agilent-084555 026652QM_RCUG_HomoSapiens",
                           "Agilent-074809 SurePrint G3 Mouse GE v2",
                           "Agilent-072363 SurePrint G3 Human GE v3")) {
    data <- limma::backgroundCorrect(raw_data, method = "normexp", offset=50)
    data <- limma::normalizeBetweenArrays(data)
  }
  else if (inherits(raw_data,"AffyBatch")) {
    data <- affy::rma(raw_data, background = T, normalize = T)
  }
  else if (inherits(raw_data, "ExpressionFeatureSet")) {
    data <- oligo::rma(raw_data, background = T, normalize = T)
  }
  else {
    stop("Normalization not implemented for class ", paste(class(raw_data), collapse=", "))
  }
  if (annotate)
    data <- annotate_data(data)
    if (clean.genes) data <- clean_genes(data)
    if (save.view) {
      if (is.null(save.dir)) stop("'save.dir' needs to be set.")
      dir.create(save.dir, recursive = TRUE, showWarnings = FALSE)
      if (inherits(data, c("ExpressionSet", "ExpressionFeatureSet"))) {
        E <- Biobase::exprs(data)
        metadata <- Biobase::pData(data)
      } else {
        E <- data$E
        metadata <- data$targets
      }
      accession <- unique(trimws(as.character(metadata[["accession"]])))
      accession <- accession[!is.na(accession) & nzchar(accession)]
      conditions <- unique(metadata$combined_condition)
      conditions <- conditions[!is.na(conditions) & nzchar(conditions)]

      for (condition in conditions) {
        sample_indices <- which(metadata$combined_condition == condition)

        E_condition <- E[, sample_indices, drop = FALSE]

        safe_condition <- gsub(
          "[^[:alnum:]_.-]+",
          "_",
          condition
        )
        safe_condition <- gsub("_+", "_", safe_condition)
        safe_condition <- gsub("^_|_$", "", safe_condition)

        output_file <- file.path(
          save.dir,
          paste0(accession, "_", safe_condition, "_norm_data.csv")
        )

        utils::write.csv(
          E_condition,
          file = output_file,
          row.names = TRUE,
          quote = FALSE
        )
      }
    }
  data
}

annotate_data <- function(data, install_missing = F) {
  supplier <- get_supplier(data)

  is_biobase <- inherits(data, c("AffyBatch", "ExpressionSet", "ExpressionFeatureSet"))
  is_limma <- inherits(data, c("EListRaw", "EList", "RGList", "MAList"))

  if (!is_biobase & !is_limma) {
    stop("Unsupported object class for annotation: ", paste(class(data), collapse=", "))
  }
  if (inherits(data, "AffyBatch")) {
    message("Affymetrix microarrays annotate after probe-set summarization (background correction, normalization).")
    return(data)
  }
  get_feature_data <- function(x) {
    if (is_biobase) {
      Biobase::fData(x)
    } else {
      if(is.null(x$genes)) {
        n <- if (!is.null(x$E)) nrow(x$E) else {
          if (!is.null(x$R)) nrow(x$R) else 0L
        }
        data.frame(row.names = seq_len(n), stringsAsFactors = F)
      } else {
        x$genes
      }
    }
  }

  set_feature_data <- function(x, fd) {
    if (is_biobase) {
      rownames(fd) <- Biobase::featureNames(x)
      Biobase::fData(x) <- fd
    } else {
      x$genes <- fd
    }
    x
  }

  feature_data <- get_feature_data(data)

  has_usable_symbol <- function(x, min_valid_frac = 0.1) {
    x <- trimws(as.character(x))
    keep <- !is.na(x) & nzchar(x)
    if (!any(keep)) return(F)
    numeric_only <- grepl(
      "^[+-]?(?:\\d+\\.?\\d*|\\.\\d+)(?:[eE][+-]?\\d+)?$",
      x[keep], perl = T
    )
    symbol_like <- grepl(
      "^[[:alnum:]][[:alnum:].-]*$", x[keep], perl = T
    )

    mean(!numeric_only & symbol_like) >= min_valid_frac
  }

  symbol_col <- intersect(c("Symbol", "SYMBOL", "GeneSymbol", "GENE_SYMBOL"), colnames(feature_data))
  existing_symbol_col <- symbol_col[
    vapply(
      symbol_col,
      function(nm) has_usable_symbol(feature_data[nm]),
      logical(1)
    )
  ]

  if (length(existing_symbol_col)) {
    source_col <- existing_symbol_col[1L]
    feature_data$Symbol <- trimws(as.character(feature_data[[source_col]]))
    message("Using existing annotation column '", source_col, "' as 'Symbol'.")
    return(set_feature_data(data, feature_data))
  }

  if (length(symbol_col)) {
    message("Existing SYMBOL/Symbol field is empty. Reannotate through platform information.")
  }
  
  annotation_map <- list(
    "Affymetrix Human genome U133 Plus 2.0 Array"          = list(pkg = "hgu133plus2.db",                  keytype = "PROBEID"),
    "Affymetrix GeneChip Human Genome U133 Plus 2.0 Array" = list(pkg = "hgu133plus2.db",                  keytype = "PROBEID"),
    "Affymetrix Human Genome U133 Plot 2.0 Array"          = list(pkg = "hgu133plus2.db",                  keytype = "PROBEID"),
    "Affymetrix Human genome U133A 2.0 Array"              = list(pkg = "hgu133a2.db",                     keytype = "PROBEID"),
    "Affymetrix HT Human Genome U133A"                     = list(pkg = "hgu133a.db",                      keytype = "PROBEID"),
    #"GeneChip PrimeView Human Gene Expression Array"       = list(pkg = "primeview.db",                    keytype = "PROBEID"), #https://support.bioconductor.org/p/130727/
    "Affymetrix Rat Genome 230A"                           = list(pkg = "rat2302.db",                      keytype = "PROBEID"),
    "Affymetrix Human Gene 2.0 ST"                         = list(pkg = "hugene20sttranscriptcluster.db",  keytype = "PROBEID"),
    "Affymetrix Mouse Transcriptome Array 1"             = list(pkg = "mta10transcriptcluster.db",       keytype = "PROBEID"),
    "Affymetrix GeneChip Mouse Gene 1.0 ST Array"          = list(pkg = "mogene10sttranscriptcluster.db",  keytype = "PROBEID"),
    "Affymetrix Rat Gene 1.0 ST"                           = list(pkg = "ragene10sttranscriptcluster.db",  keytype = "PROBEID"),
    "Illumina HumanHT-12 v4.0"                             = list(pkg = "illuminaHumanv4.db",              keytype = "PROBEID"),
    "Illumina HumanHT-12 v4"                               = list(pkg = "illuminaHumanv4.db",              keytype = "PROBEID"),  
    "Illumina HumanHT-12 v3"                               = list(pkg = "illuminaHumanv3.db",              keytype = "PROBEID"),
    "Illumina HumanHT-12 v3.0"                             = list(pkg = "illuminaHumanv3.db",              keytype = "PROBEID"),
    "Illumina MouseRef-8 v2"                               = list(pkg = "illuminaMousev2.db",              keytype = "PROBEID"),
    "Illumina MouseRef-8 v2.0"                             = list(pkg = "illuminaMousev2.db",              keytype = "PROBEID")
  )
  args <- annotation_map[[supplier]]
  
  if (!is.null(args)) {
    if (!requireNamespace(args$pkg, quietly = T)){
      if (!install_missing) {
        stop("Require annotation package '", args$pkg, "' is not installed. Install with: \n BiocManager::install('", args$pkg, "')")
      }
      BiocManager::install(args$pkg, ask=F, update=F)
    }

    db <- get(args$pkg, envir = loadNamespace(args$pkg))
    
    probe_col <- intersect(
      c("ProbeName", "PROBE_ID", "ProbeID", "probe_id", "Probe_Id"), colnames(feature_data)
    )

    if (is_biobase) {
      probe_ids <- Biobase::featureNames(data)
      if (length(probe_col)) {
        candidate <- as.character(feature_data[[probe_col[1L]]])
        if (length(candidate) ==length(probe_ids) && mean(!is.na(candidate) & nzchar(candidate)) > 0.9) {
          probe_ids <- candidate
        }
      } 
    } else {
      if (!length(probe_col)) {
        stop("No probe-ID field found. Available columns: ", paste(colnames(feature_data), collapse = ", "))
      }
      probe_ids <- as.character(feature_data[[probe_col[1L]]])
    }
    valid <- !is.na(probe_ids) & nzchar(probe_ids)
    symbols <- rep(NA_character_, length(probe_ids))
    symbols[valid] <- unname(
      AnnotationDbi::mapIds(
        db, keys = probe_ids[valid], column = "SYMBOL", keytype = args$keytype, multiVals = "first"
      )
    )
    feature_data$Symbol <- symbols
    message("Mapped", sum(!is.na(symbols)), " / ", length(symbols), " features to symbols using ", args$pkg, ".")
    return(set_feature_data(data, feature_data))
  }

  if (supplier %in% c(
    "Agilent-084555 026652QM_RCUG_HomoSapiens",
    "Agilent-074809 SurePrint G3 Mouse GE v2",
    "Agilent-072363 SurePrint G3 Human GE v3"
  )) {
    gene_col <- intersect(
      c("GeneName", "GENE_NAME", "GeneSymbol", "SYMBOL", "Symbol"), colnames(feature_data)
    )
    if (!length(gene_col)) {
      stop("No gene name field avaliable.")
    }
    feature_data$Symbol <- trimws(as.character(feature_data[[gene_col[1L]]]))
    return(set_feature_data(data, feature_data))
  } else if (supplier %in% c(
    "Agilent miRNA microarray",
    "Agilent-070155 Mouse miRNA Microarray",
    "Agilent-070156 Human_miRNA_V21.0_Microarray",
    "Agilent-070156 Human miRNA",
    "Agilent-050340 Custom Rat miRNA Microarry",
    "Agilent-019159 Rat miRNA"
  )) {
    id_col <- intersect(
      c("miRNA_ID", "SystematicName", "ProbeName", "ProbeID"), colnames(feature_data)
    )
    if (!length(id_col)) {
      stop("No miRNA identifier field found.")
    }
    feature_data$Symbol <- trimws(as.character(feature_data[[id_col[1L]]]))
    message("Stored platform miRNA IDs from '", id_col[1L], "' in 'Symbol'")
    return(set_feature_data(data, feature_data))
  } else if (supplier == "Affymetrix GeneChip miRNA 3" || supplier == "Affymetrix GeneChip miRNA 3.0") {
    annot_file <- file.path(
      "/usr/local/storage/data_microarray/annotations/", "GPL16384_miRNA-3_1-st-v1.annotations.20140513.csv"
    )
    if (!file.exists(annot_file))
      stop("Affymetrix miRNA 3.0 annotation csv not found.")
    annot <- read.csv(annot_file, comment.char = "#", stringsAsFactors = F, check.names = F)
    probe_ids <- if (is_biobase) {
      Biobase::featureNames(data)
    } else {
      as.character(feature_data$ProbeName)
    }
    feature_data$Symbol <- annot[["Transcript ID(Array Design)"]][
      match(probe_ids, annot[["Probe Set ID"]])
    ]
    return(set_feature_data(data, feature_data))
  } else stop("chip not implemented!")
}

filter_illumina <- function(x, detection_name = "Detection PValue") {
  D <- x$other[[detection_name]]
  if (!identical(dim(x$E), dim(D)) || !identical(colnames(x$E), colnames(D)) || !identical(rownames(x$E), rownames(D))) {
    stop("Expression and expression matrix should have identical dimensions (permutation sensitive)")
  }

  keep <- rowSums(
    is.finite(x$E) & is.finite(D) & D>=0 & D<=1
  ) == ncol(x$E)
  dropped <- which(!keep)

  if (length(dropped)) {
    probe_label <- if (!is.null(x$genes) && "ProbeName" %in% colnames(x$genes)) {
      as.character(x$genes$ProbeName[dropped])
    } else {
      rownames(x$E)[dropped]
    }
    message(length(dropped), " features removed (missing, non-finite or invalid detection pvalue).\nFeatures: ", 
    paste(head(probe_label, 5L), collapse=", "), if(length(probe_label)>5L) " ..." else ""
    )
    x$E <- x$E[keep, , drop=F]
    if (!is.null(x$genes)) {
      x$genes <- x$genes[keep, , drop=F]
    }

    if (!is.null(x$other)) {
      x$other <- lapply(
        x$other, function(m) {
          if (is.matrix(m) && nrow(m) == length(keep)) {
            m[keep, , drop=F]
          } else {
            m
          }
        }
      )
    }
  }
  attr(x, "pre_nec_filter") <- list(
    detection_name = detection_name,
    n_features_in = length(keep),
    n_features_kept = sum(keep),
    n_features_removed = sum(!keep),
    removed_feature_indices = dropped
  )

  x
}

read_agilent_mirna_safe <- function(targets, verbose = FALSE) {
  stopifnot("FileName" %in% colnames(targets))

  files <- as.character(targets$FileName)

  if (!length(files)) {
    stop("No files supplied in targets$FileName.")
  }

  if (!all(file.exists(files))) {
    stop(
      "Missing AFE file(s): ",
      paste(basename(files[!file.exists(files)]), collapse = ", ")
    )
  }

  find_afe_feature_header <- function(file) {
    lines <- readLines(file, warn = FALSE)

    if (!length(lines)) {
      stop("Empty AFE file: ", basename(file))
    }

    split_tab <- function(x) {
      strsplit(x, "\t", fixed = TRUE)[[1L]]
    }

    # Standard full AFE export:
    # a line consisting of "FEATURES" marks the start of the probe-level table.
    features_idx <- which(
      toupper(trimws(lines)) == "FEATURES"
    )

    if (length(features_idx)) {
      start <- features_idx[1L] + 1L

      candidate_idx <- seq.int(
        start,
        min(length(lines), start + 10L)
      )

      for (i in candidate_idx) {
        fields <- trimws(split_tab(lines[i]))

        is_feature_header <- any(
          tolower(fields) %in% c(
            "featurenum",
            "probename",
            "systematicname",
            "controltype",
            "row",
            "col"
          )
        ) &&
          length(fields) >= 3L

        if (is_feature_header) {
          return(i)
        }
      }

      stop(
        "Found a FEATURES section but not its column header in: ",
        basename(file), "\n",
        "Lines around FEATURES:\n",
        paste(
          sprintf(
            "%05d: %s",
            candidate_idx,
            lines[candidate_idx]
          ),
          collapse = "\n"
        )
      )
    }

    # Compact export:
    # locate a tab-separated line containing recognized feature-level fields.
    feature_field_candidates <- c(
      "featurenum",
      "probename",
      "systematicname",
      "controltype",
      "row",
      "col",
      "genename"
    )

    for (i in seq_along(lines)) {
      fields <- tolower(trimws(split_tab(lines[i])))

      n_matches <- sum(fields %in% feature_field_candidates)

      if (length(fields) >= 3L && n_matches >= 2L) {
        return(i)
      }

      # Handles compact files with headers such as:
      # SystematicName  ControlType  gTotalGeneSignal ...
      if (
        length(fields) >= 3L &&
        "systematicname" %in% fields &&
        any(grepl("^g(total|mean|processed)", fields))
      ) {
        return(i)
      }
    }

    preview_n <- min(80L, length(lines))

    stop(
      "Could not find a per-feature AFE table header in: ",
      basename(file), "\n",
      "First ", preview_n, " lines:\n",
      paste(
        sprintf(
          "%05d: %s",
          seq_len(preview_n),
          lines[seq_len(preview_n)]
        ),
        collapse = "\n"
      )
    )
  }

  get_available_columns <- function(file) {
    header_idx <- find_afe_feature_header(file)

    x <- limma::read.columns(
      file,
      required.col = NULL,
      text.to.search = "",
      skip = header_idx - 1L,
      sep = "\t",
      quote = "\"",
      stringsAsFactors = FALSE,
      flush = TRUE
    )
    colnames(x)
  }

  available_by_file <- lapply(files, get_available_columns)

  common_columns <- Reduce(intersect, available_by_file)

  if (!length(common_columns)) {
    stop("No common columns across the supplied Agilent AFE files.")
  }

  select_first <- function(candidates, available, required = TRUE, role = "") {
    hit <- candidates[
      match(tolower(candidates), tolower(available), nomatch = 0L) > 0L
    ]

    if (!length(hit)) {
      if (required) {
        stop(
          "No usable ", role, " field found.\n",
          "Tried: ", paste(candidates, collapse = ", "), "\n",
          "Available common fields: ",
          paste(common_columns, collapse = ", ")
        )
      }
      return(NULL)
    }

    available[
      match(
        tolower(hit[1L]),
        tolower(available)
      )
    ]
  }

  probe_name_col <- select_first(
    candidates = c(
      "ProbeName",
      "SystematicName",
      "FeatureNum",
      "ID"
    ),
    available = common_columns,
    role = "probe identifier"
  )

  gene_name_col <- select_first(
    candidates = c(
      "GeneName",
      "miRNA_ID",
      "miRNA",
      "GeneSymbol"
    ),
    available = common_columns,
    required = FALSE,
    role = "gene/miRNA annotation"
  )

  control_type_col <- select_first(
    candidates = c(
      "ControlType",
      "ControlTypeName"
    ),
    available = common_columns,
    required = FALSE,
    role = "control-probe annotation"
  )

  total_signal_col <- select_first(
    candidates = c(
      "gTotalGeneSignal",
      "gProcessedSignal",
      "gMeanSignal",
      "gMedianSignal",
      "gMeanSignalNoBkg",
      "gTotalProbeSignal"
    ),
    available = common_columns,
    role = "primary expression signal"
  )

  mean_signal_col <- select_first(
    candidates = c(
      "gMeanSignal",
      "gProcessedSignal",
      "gTotalGeneSignal",
      "gMedianSignal"
    ),
    available = common_columns,
    required = FALSE,
    role = "mean signal"
  )

  processed_signal_col <- select_first(
    candidates = c(
      "gProcessedSignal",
      "gTotalGeneSignal",
      "gMeanSignal",
      "gMedianSignal"
    ),
    available = common_columns,
    required = FALSE,
    role = "processed signal"
  )

  detection_col <- select_first(
    candidates = c(
      "gIsGeneDetected",
      "gIsWellAboveBG",
      "gIsWellAboveBGFlag",
      "gDetectionPValue",
      "gPValue"
    ),
    available = common_columns,
    required = FALSE,
    role = "detection field"
  )

  background_col <- select_first(
    candidates = c(
      "gBGMedianSignal",
      "gBGMeanSignal",
      "gBGPValue",
      "gBGUsed"
    ),
    available = common_columns,
    required = FALSE,
    role = "background field"
  )

  expression_columns <- list(
    TGS = total_signal_col
  )

  if (!is.null(mean_signal_col)) {
    expression_columns$meanS <- mean_signal_col
  }

  if (!is.null(processed_signal_col)) {
    expression_columns$procS <- processed_signal_col
  }

  annotation_columns <- c(
    control_type_col,
    probe_name_col,
    gene_name_col
  )

  annotation_columns <- Filter(
      Negate(is.null),
      c(
        control_type_col,
        probe_name_col,
        gene_name_col
      )
    )

  other_columns <- list()

  if (!is.null(detection_col)) {
    other_columns$Detection <- detection_col
  }

  if (!is.null(background_col)) {
    other_columns$Background <- background_col
  }

  message(
    "Agilent miRNA columns selected:\n",
    "  Probe ID: ", probe_name_col, "\n",
    "  Primary signal: ", total_signal_col, "\n",
    "  Annotation: ", paste(annotation_columns, collapse = ", "), "\n",
    "  Other: ",
    if (length(other_columns)) {
      paste(names(other_columns), "=", unlist(other_columns), collapse = ", ")
    } else {
      "<none>"
    }
  )

  AgiMicroRna:::read.agiMicroRna(
    targets = targets,
    columns = expression_columns,
    other.columns = other_columns,
    annotation = annotation_columns,
    verbose = verbose
  )
}

run_dea <- function(
  data,
  group_col,
  group_levels = NULL,
  contrast_str = NULL,
  covariate_cols = NULL,
  p_adjust_method = "BH",
  save.view = F,
  save.dir = NULL
) {
  is_eset <- inherits(data, c("ExpressionSet", "ExpressionFeatureSet"))
  is_limma <- inherits(data, c("EList", "EListRaw", "uRNAList"))

  if (!is_eset && !is_limma) {
    stop(
      "Unsupported class for DEA"
    )
  }
  if (is_eset) {
    expr <- Biobase::exprs(data)
    genes <- Biobase::fData(data)
    meta <- Biobase::pData(data)
  } else {
    expr <- data$E
    genes <- data$genes
    meta <- data$targets
  }

  sample_ids <- colnames(expr)

  missing_metadata <- setdiff(sample_ids, rownames(meta))

  if (length(missing_metadata)) {
    stop(
      "Expression samples are absent from metadata:\n",
      paste(missing_metadata, collapse = ", ")
    )
  }

  meta <- meta[sample_ids, , drop=F]

  if (!identical(rownames(meta), colnames(expr))) {
    stop("metadata is not aligned to expression data.")
  }
  
  missing_covariates <- setdiff(
    covariate_cols %||% character(0), colnames(meta)
  )

  if (length(missing_covariates)) {
    stop("Missing covariate columns: ", paste(missing_covariates, collapse = ", "))
  }

  group_values <- as.character(meta[[group_col]])

  if (is.null(group_levels)) {
    group_levels <- unique(group_values)
  }

  group <- factor(group_values, levels = group_levels)
  if (anyNA(group)) {
    stop(
      "Missing/unrecognized group values in ", group_col, ".\n",
      "Observed: ", paste(unique(group_values), collapse = ", "), "\n",
      "Requested: ", paste(group_levels, collapse = ", ")
    )
  }
  if (nlevels(group)<2L) stop("DEA needs at least 2 levels.")
  
  design_data <- meta
  design_data$group <- group

  design_formula <- if (is.null(covariate_cols) || !length(covariate_cols)) {
    stats::as.formula("~ 0 + group")
  } else {
    stats::as.formula(
      paste0("~ 0 + group + ", paste(covariate_cols, collapse = " + "))
    )
  }
  design <- stats::model.matrix(design_formula, data=design_data)

  colnames(design) <- sub("^group", "", colnames(design))
  colnames(design) <- make.names(colnames(design))

  if (qr(design)$rank < ncol(design)) stop("The DEA design matrix is not full rank.")
  
  fit <- limma::lmFit(object = expr, design = design)

  if (is.null(contrast_str)) {
    if (nlevels(group) != 2L) stop("Require 2 group levels for automatic contrast detection.")
    lev <- make.names(levels(group))
    contrast_str <- paste0(lev[2L], "-", lev[1L])
  }
  contrast_matrix <- limma::makeContrasts(contrasts = contrast_str, levels = design)
  fit2 <- limma::contrasts.fit(fit, contrasts = contrast_matrix)
  fit2 <- limma::eBayes(fit2)

  results <- limma::topTable(
    fit2, coef = 1L, number = Inf, adjust.method = p_adjust_method, sort.by = "P"
  )

  if (!is.null(genes) && nrow(genes) == nrow(expr)) {
    if (is.null(rownames(genes))) rownames(genes) <- rownames(expr)
    genes <- genes[match(rownames(results), rownames(genes)), , drop=F]
    results <- cbind(genes, results)
  }

  if (save.view) {
    if (is.null(save.dir)) stop("'save.dir' must be provided when view should be saved.")
    dir.create(save.dir, recursive = TRUE, showWarnings = FALSE)
    if (inherits(data, c("ExpressionSet", "ExpressionFeatureSet"))) {
      metadata <- Biobase::pData(data)
    } else {
      metadata <- data$targets
    }
    accession <- unique(trimws(as.character(metadata$accession)))
    accession <- accession[!is.na(accession) & nzchar(accession)]

    output_file <- file.path(
      save.dir,
      paste0(paste(accession, sep = "_"), "_dea.csv")
    )
    utils::write.csv(
      results,
      file = output_file,
      row.names = TRUE,
      quote = FALSE
    )
  }

  list(fit = fit2, design = design, metadata = meta, contrast = contrast_str, results = results)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

clean_genes <- function(data, symbol_col = "Symbol", remove_controls = T, verbose = T) {
  is_eset <- inherits(data, c("ExpressionSet", "ExpressionFeatureSet"))
  is_limma <- inherits(data, c("EListRaw", "EList", "RGList"))

  if (!is_eset && !is_limma) stop("Unsupported object")
  
  if (is_eset) {
    E <- Biobase::exprs(data)
    genes <- Biobase::fData(data)
  } else if (is_limma) {
    E <- data$E
    genes <- data$genes
  }

  if (!symbol_col %in% colnames(genes)) stop("Gene-symbol column '", symbol_col, "' is missing.")
  
  symbols <- trimws(as.character(genes[[symbol_col]]))

  supplier <- as.character(get_supplier(data))[1L]
  
  is_mirna_chip <- !is.na(supplier) && nzchar(supplier) && grepl("mirna", supplier, ignore.case = TRUE)
  missing_symbol <- is.na(symbols) | !nzchar(symbols)

  if (!is_mirna_chip) {
    missing_symbol <- missing_symbol |
      symbols %in% c("---", "NA", "N/A", "NULL", "null")
  }

  is_control <- rep(FALSE, nrow(E))

  if (remove_controls) {
    if ("ControlType" %in% colnames(genes)) {
      control_type <- trimws(as.character(genes$ControlType))
      control_type_num <- as.numeric(control_type)

      det_controls <- !is.na(control_type_num) & control_type_num != 0
      det_text <- tolower(control_type) %in% c(
        "pos", "positive", "positive control",
        "neg", "negative", "negative control",
        "control", "not probe", "not_probe", "ignore"
      )
      det_text[is.na(det_text)] <- FALSE
      is_control <- is_control | det_controls | det_text
    }

    if ("Status" %in% colnames(genes)) {
      status <- trimws(tolower(as.character(genes$Status)))
      known_status <- !is.na(status) & nzchar(status)

      status_control <- known_status & status %in% c(
        "negative", "positive", "control", "neg", "pos", "housekeeping", "blank"
      )

      has_gene_status <- any(status == "gene", na.rm = T)

      if (has_gene_status) status_control <- status_control | (known_status & status != "gene")
      
      status_control[is.na(status_control)] <- FALSE
      is_control <- is_control | status_control
    }

    is_affy <- !is.na(supplier)  && grepl("^Affymetrix", supplier, ignore.case = T)
    if (is_affy) {
      id_col <- intersect(c("ProbeName", "PROBE_ID", "ProbeID", "ID_REF"), colnames(genes))
      feature_id <- if (length(id_col)) as.character(genes[[id_col[1L]]]) else rownames(E)
      affx_control <- grepl("^AFFX-", feature_id, ignore.case = TRUE)
      affx_control[is.na(affx_control)] <- FALSE
      is_control <- is_control | grepl("^AFFX-", feature_id, ignore.case = T)
    }

    is_agilent <- !is.na(supplier) && grepl("^Agilent", supplier, ignore.case = T)
    if (is_agilent && !"ControlType" %in% colnames(genes)) {
      id_col <- intersect(c("SystematicName", "ProbeName", "ProbeID", "ID_REF"), colnames(genes))
      feature_id <- if (length(id_col)) trimws(as.character(genes[[id_col[1L]]])) else rownames(E)
            agilent_control <- grepl(
        "^(DarkCorner|SCorner|NC[0-9]+_)",
        feature_id,
        ignore.case = TRUE,
        perl = TRUE
      )
      agilent_control[is.na(agilent_control)] <- FALSE
      is_control <- is_control | grepl("^(DarkCorner|SCorner|NC[0-9]+_)", feature_id, ignore.case = T, perl = T)
    }

    if (!is.na(supplier) && supplier %in% c("Affymetrix GeneChip miRNA 3.0", "Affymetrix GeneChip miRNA 3") &&
      "SequenceType" %in% colnames(genes)) {
      sequence_type <- trimws(tolower(as.character(genes$SequenceType)))
      is_control <- is_control | (!is.na(sequence_type) & nzchar(sequence_type) & sequence_type != "mirna")
    }
  }

  keep <- !missing_symbol & !is_control
  keep[is.na(keep)] <- FALSE
  if (!any(keep)) stop("No features retained after Control removal.")
  
  audit <- data.frame(
    total_features = nrow(E),
    retained_features = sum(keep),
    removed_features = sum(!keep),
    removed_missing_symbol = sum(missing_symbol),
    removed_controls = sum(is_control),
    removed_both = sum(missing_symbol & is_control),
    stringsAsFactors = T
  )

  if (verbose) {
    message(
      "Feature cleaning:\n",
      "Input features: ", nrow(E), "\n",
      "Removed: ", sum(!keep), "\n",
      "no mapped Symbol: ", sum(missing_symbol), "\n",
      "Controls: ",  sum(is_control), "\n",
      "Remaining genes:", sum(keep)
    )
  }

  removed <- data.frame(
    feature_id = rownames(E)[!keep],
    symbol = symbols[!keep],
    missing_symbol = missing_symbol[!keep],
    control = is_control[!keep],
    stringsAsFactors = T
  )

  if (is_eset) {
    cleaned <- data[keep, ]
  } else {
    cleaned <- data
    if (!is.null(cleaned$E) && is.matrix(cleaned$E)) {
      if (nrow(cleaned$E) != length(keep)) {
        stop("`$E` is not aligned with the feature dimension.")
      }

      cleaned$E <- cleaned$E[keep, , drop = FALSE]
    }
    for (channel_name in c("R", "G", "Rb", "Gb", "weights")) {
      channel <- cleaned[[channel_name]]

      if (!is.null(channel) && is.matrix(channel)) {
        if (nrow(channel) != length(keep)) {
          stop(
            "`$", channel_name,
            "` is not aligned with the feature dimension."
          )
        }

        cleaned[[channel_name]] <- channel[keep, , drop = FALSE]
      }
    }
    if (!is.null(cleaned$genes)) cleaned$genes <-cleaned$genes[keep, , drop = F]
    if (!is.null(cleaned$other)) cleaned$other <- lapply(cleaned$other, function(x) {
      if (is.matrix(x) && nrow(x) == length(keep)) x[keep, , drop = F] else x
    })
  }

  attr(cleaned, "feature_cleaning_audit") <- audit
  attr(cleaned, "removed_features") <- removed
  attr(cleaned, "feature_filter") <- list(
    required_symbol_column = symbol_col,
    removed_controls = remove_controls
  )

  cleaned
}

combine_chip_datasets <- function(..., dataset_names = NULL) {
  datasets <- list(...)

  if (length(datasets) < 2L) {
    return(datasets[0L])
  }

  if (is.null(dataset_names)) {
    dataset_names <- paste0("dataset_", seq_along(datasets))
  }

  if (
    (length(dataset_names) != length(datasets)) ||
      (anyDuplicated(dataset_names)) ||
      (anyNA(dataset_names)) ||
      (any(!nzchar(dataset_names)))
  ) {
    stop("'dataset_names' must constain one unique non-empty name per dataset.")
  }

  get_family <- function(x) {
    if (inherits(x, "uRNAList")) return ("uRNAlist")
    if (inherits(x, c("ExpressionSet", "ExpressionFeatureSet"))) return("ExpressionSet")
    if (inherits(x, c("EList", "EListRaw"))) return("EList")
    stop("Unsupported object class. Currently only support for normalized data.")
  }

  families <- vapply(datasets, get_family, character(1))

  if (length(unique(families)) != 1L) {
    stop("At least one dataset has different datatype from the others. Concatination cannot happen.")
  }

  get_meta <- function(x, family) {
    meta <- if (family == "ExpressionSet") Biobase::pData(x) else x$targets

    if (is.null(meta) || !is.data.frame(meta)) stop("Metadata not extracted correctly.")
    meta
  }
  get_chip <- function(x, family) {
    meta <- get_meta(x, family)
    if (!"chip" %in% colnames(meta)) stop("'chip' needs to be present in metadata.")
    chip <- unique(trimws(as.character(meta$chip)))
    chip <- chip[!is.na(chip) & nzchar(chip)]
    if (length(chip) != 1L) stop("More than one non-empty chip within one datasets. Please double check.")
    chip
  }
  chips <- mapply(get_chip, datasets, families, USE.NAMES = F)

  if (length(unique(chips)) != 1L) stop("Datasets not from the same chip!")
  
  extract_standardized <- function(x, family, dataset_name) {
    if (family == "ExpressionSet") {
      E <- Biobase::exprs(x)
      genes <- Biobase::fData(x)
      meta <- Biobase::pData(x)
      feature_ids <- Biobase::featureNames(x)
    } else {
      E <- x$E
      genes <- x$genes
      meta <- x$targets
      feature_ids <- if (!is.null(genes) && "ProbeName" %in% colnames(genes)) {
        trimws(as.character(genes$ProbeName))
      } else {
        rownames(E)
      }
    }

    if (is.null(E)) stop("Expression matrix from dataset '", dataset_name, "' couldn't be extracted.")
    if (
      is.null(feature_ids) || 
        length(feature_ids) != nrow(E) || 
        anyNA(feature_ids) || 
        any(!nzchar(feature_ids))
    ) stop("Invalid feature identifiers in '", dataset_name, "'")
    if (anyDuplicated(feature_ids)) stop("Duplicated feature ids in dataset '", dataset_name, "'")
    if (is.null(genes)) {
      genes <- data.frame(
        ProbeName = feature_ids,
        row.names = feature_ids,
        stringsAsFactors = F
      )
    }

    if (is.null(rownames(genes))) rownames(genes) <- feature_ids

    sample_ids <- colnames(E)

    if (
      is.null(sample_ids) ||
      anyNA(sample_ids) ||
      any(!nzchar(sample_ids))
    ) stop("Invalid sample names in '", dataset_name, "'")
    if (anyDuplicated(sample_ids)) stop("Duplicated sample names in dataset '", dataset_name, "'")
    if (is.null(rownames(meta))) stop("Metadata has no rowname in '", dataset_name, "'")
    
    missing_meta <- setdiff(sample_ids, rownames(meta))

    if (length(missing_meta)) stop("Metadata is missing samples in '", dataset_name, "'")
    
    meta <- meta[sample_ids, , drop = F]

    if (!identical(rownames(meta), sample_ids)) stop("Metadata and expression samples not aligned in '", dataset_name, "'")

    remove_tech_columns <- c("FileName", "Treatment", "GErep")
    present_columns <- intersect(remove_tech_columns,colnames(meta))

    if (length(present_columns)) meta <- meta[, !colnames(meta) %in% present_columns, drop = FALSE]
    
    rownames(E) <- feature_ids

    list(
      E = E,
      genes = genes,
      meta = meta,
      feature_ids = feature_ids
    )
  }

  res <- Map(
    extract_standardized,
    datasets,
    families,
    dataset_names
  )

  ref_features <- res[[1L]]$feature_ids
  max_feature_diff <- 0.001

  feature_check <- lapply(seq_along(res), function(i) {
    curr_features <- res[[i]]$feature_ids
    missing <- setdiff(ref_features, curr_features)
    extra <- setdiff(curr_features, ref_features)

    denominator <- max(length(ref_features), length(curr_features))
    diff <- length(union(missing, extra))

    data.frame(
      dataset = dataset_names[i],
      ref_features = length(ref_features),
      dataset_features = length(curr_features),
      shared = length(intersect(ref_features, curr_features)),
      missing = length(missing),
      extra = length(extra),
      diff_features = diff,
      difference_fraction = diff / denominator,
      stringsAsFactors = F
    )

  })

  feature_audit <- do.call(rbind, feature_check)
  too_different <- feature_audit$difference_fraction > max_feature_diff
  if (any(too_different)) stop("Feature sets differ by more than ", max_feature_diff*100, "% from reference.")
  
  common_features <- Reduce(intersect, lapply(res, '[[', "feature_ids"))
  common_features <- ref_features[ref_features %in% common_features]

  message("Combining datasets using ", length(common_features), " shared probe sets.")

  res <- Map(
    function(x, dataset_name) {
      idx <- match(ref_features, x$feature_ids)
      E <- x$E[idx, , drop=F]
      genes <- x$genes[idx, , drop=F]

      rownames(E) <- ref_features
      rownames(genes) <- ref_features

      org_samples_ids <- colnames(E)
      comb_sample_ids <- paste(dataset_name, org_samples_ids, sep = "__")

      if (anyDuplicated(comb_sample_ids)) stop("Duplicated sample names generated")
      
      colnames(E) <- comb_sample_ids
      meta <- x$meta
      meta$dataset <- dataset_name
      meta$original_sample_id <- org_samples_ids
      rownames(meta) <- comb_sample_ids

      list(E = E, genes = genes, meta = meta)
    },
    res,
    dataset_names
  )

  standardize_meta_colnames <- function(meta) {
    nm <- colnames(meta)

    # Standardize punctuation and whitespace first.
    standardized <- nm |>
      trimws() |>
      gsub("[[:space:]/.]+", "_", x = _) |>
      gsub("_+", "_", x = _) |>
      gsub("^_|_$", "", x = _)

    # Correct known spelling variants after punctuation standardization.
    standardized[standardized == "sample_nr"] <- "sample_nr"

    colnames(meta) <- standardized

    meta
  }

  res <- lapply(
    res,
    function(x) {
      x$meta <- standardize_meta_colnames(x$meta)
      x
    }
  )

  E_comb <- do.call(cbind, lapply(res, '[[', "E"))
  targets_comb <- do.call(rbind, lapply(res, '[[', "meta"))
  genes_comb <- res[[1L]]$genes

  combined <- list(
    E = E_comb,
    genes = genes_comb,
    targets = targets_comb
  )

  class(combined) <- "EList"

  attr(combined, "chip") <- chips[1L]
  attr(combined, "combined_datasets") <- dataset_names
  attr(combined, "combination_level") <- "same_chip/probe_level"
  attr(combined, "preprocessing_status") <- "normalized_separately_then_combined"
  attr(combined, "input_object_family") <- families[1L]
  combined
}

run_pca <- function(data, n_top = 1000L, center = T, scale. = F) {
  is_eset <- inherits(data, c("ExpressionSet", "ExpressionFeatureSet"))
  is_elist <- inherits(data, c("EList", "EListRaw"))

  if (!is_eset && !is_elist) stop("Unsupported object class")
  
  E <- if (is_eset) Biobase::exprs(data) else data$E
  metadata <- if (is_eset) Biobase::pData(data) else data$targets

  missing_meta <- setdiff(colnames(E), rownames(metadata))
  if (length(missing_meta)) stop("Samples absent from metadata")
  
  metadata <- metadata[colnames(E), , drop = F]
  feature_var <- apply(E, 1L, stats::var, na.rm = T)
  keep <- is.finite(feature_var) & feature_var > 0
  E <- E[keep, , drop = F]
  feature_var <- feature_var[keep]
  n_top <- as.integer(n_top)
  n_selected <- min(n_top, nrow(E))
  selected_features <- names(sort(feature_var, decreasing = T))[seq_len(n_selected)]
  sample_by_feature <- t(E[selected_features, , drop = F])
  pca_fit <- stats::prcomp(sample_by_feature, center = center, scale. = scale.)
  variance_explained <- (pca_fit$sdev^2 / sum(pca_fit$sdev^2))*100

  scores <- as.data.frame(pca_fit$x, check.names = F)
  scores$sample_id <- rownames(scores)

  metadata$sample_id <- rownames(metadata)
  score_metadata <- merge(scores, metadata, by = "sample_id", all.x = T, sort = F)
  score_metadata <- score_metadata[match(rownames(pca_fit$x), score_metadata$sample_id), , drop = F]

  list(
    fit = pca_fit,
    scores = score_metadata,
    loadings = pca_fit$rotation,
    variance_explained = variance_explained,
    selected_features = selected_features,
    n_features_in = nrow(E),
    n_features_out = n_selected,
    center = center,
    scale = scale.
  )
}

pca_plot <- function(
  pca_result,
  colour_by,
  shape_by = NULL,
  pc_x = 1L,
  pc_y = 2L,
  point_size = 1.5
) {
  scores <- pca_result$scores
  variance <- pca_result$variance_explained

  x_name <- paste0("PC", pc_x)
  y_name <- paste0("PC", pc_y)

  if (!all(c(x_name, y_name, colour_by) %in% colnames(scores))) {
    stop("Requested PC or metadata column is absent from PCA scores.")
  }

  colour_values <- as.factor(scores[[colour_by]])
  colours <- grDevices::hcl.colors(
    nlevels(colour_values),
    palette = "Dark 3"
  )

  point_col <- colours[as.integer(colour_values)]

  plot(
    scores[[x_name]],
    scores[[y_name]],
    col = point_col,
    pch = 19,
    cex = point_size,
    xlab = sprintf("%s (%.1f%% variance)", x_name, variance[pc_x]),
    ylab = sprintf("%s (%.1f%% variance)", y_name, variance[pc_y]),
    main = "PCA"
  )

  if (!is.null(shape_by)) {
    if (!shape_by %in% colnames(scores)) {
      stop("`shape_by` is absent from PCA scores.")
    }

    shape_values <- as.factor(scores[[shape_by]])

    if (nlevels(shape_values) > 25L) {
      stop("`shape_by` has more than 25 levels; choose another variable.")
    }

    plot(
      scores[[x_name]],
      scores[[y_name]],
      col = point_col,
      pch = as.integer(shape_values) + 14L,
      cex = point_size,
      xlab = sprintf("%s (%.1f%% variance)", x_name, variance[pc_x]),
      ylab = sprintf("%s (%.1f%% variance)", y_name, variance[pc_y]),
      main = "PCA"
    )

    legend(
      "topright",
      legend = levels(shape_values),
      pch = seq_len(nlevels(shape_values)) + 14L,
      title = shape_by,
      bty = "n"
    )
  }

  legend(
    "topleft",
    legend = levels(colour_values),
    col = colours,
    pch = 19,
    title = colour_by,
    bty = "n"
  )
}

urna_to_elist <- function(x) {
  if (!inherits(x, "uRNAList")) stop("'x' must inherit from 'uRNAList'.")
  
  E <- if (!is.null(x$E) && is.matrix(x$E)) {
    x$E
  } else if (!is.null(x$TGS) && is.matrix(x$TGS)) {
    x$TGS
  } else {
    stop("uRNAList has neither a valid `$E` nor `$TGS` matrix.")
  }
  if (is.null(x$genes) || !is.data.frame(x$genes)) {
    stop("uRNAList has no valid `$genes` data.frame.")
  }

  if (nrow(x$genes) != nrow(E)) {
    stop(
      "`$genes` and expression matrix are not aligned: ",
      "nrow(genes) = ", nrow(x$genes),
      "; nrow(E) = ", nrow(E)
    )
  }

  if (is.null(x$targets) || !is.data.frame(x$targets)) {
    stop("uRNAList has no valid `$targets` data.frame.")
  }

  if (is.null(rownames(x$targets))) {
    stop("uRNAList `$targets` must have sample IDs as row names.")
  }

  missing_meta <- setdiff(colnames(E), rownames(x$targets))
  if (length(missing_meta)) {
    stop(
      "uRNAList expression samples missing from `$targets`: ",
      paste(missing_meta, collapse = ", ")
    )
  }

  targets <- x$targets[colnames(E), , drop = FALSE]

  # Ensure a stable, feature-level ID for later annotation/merging.
  if (is.null(rownames(E)) || anyNA(rownames(E)) || any(!nzchar(rownames(E)))) {
    if ("ProbeName" %in% colnames(x$genes)) {
      rownames(E) <- trimws(as.character(x$genes$ProbeName))
    } else if ("SystematicName" %in% colnames(x$genes)) {
      rownames(E) <- trimws(as.character(x$genes$SystematicName))
    } else {
      stop("Expression matrix has no row names and no usable probe-ID column.")
    }
  }

  out <- list(
    E = E,
    genes = x$genes,
    targets = targets
  )

  class(out) <- "EList"
  if (!is.null(x$other)) {
    out$other <- x$other
  }

  # Retain non-expression matrices only if they remain feature-aligned.
  if (!is.null(x$other)) {
    out$other <- x$other
  }

  attr(out, "source_class") <- class(x)
  attr(out, "preprocessing_done") <- TRUE
  attr(out, "preprocessing_status") <- "AgiMicroRna_tgsNormalization"
  attr(out, "expression_assay") <- if (!is.null(x$E) && is.matrix(x$E)) {
    "uRNAList$E"
  } else {
    "uRNAList$TGS"
  }

  out
}