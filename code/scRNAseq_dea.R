suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

usage <- function() {
  cat(
    paste0(
      "\nRun limma-voom pseudobulk differential expression.\n\n",
      "Required:\n",
      "  --counts FILE                 TSV: genes × pseudobulk profiles.\n",
      "  --metadata FILE               TSV: pseudobulk profiles × metadata.\n",
      "  --output FILE                 Output TSV path.\n",
      "  --group-col COLUMN            Metadata column defining condition/group.\n",
      "  --group-levels A,B            Ordered groups. Tests B - A.\n\n",
      "Optional:\n",
      "  --celltype-col COLUMN         Metadata cell-type column.\n",
      "  --celltype VALUE              Analyse only this cell type.\n",
      "  --covariates A,B              Optional metadata covariates.\n",
      "  --n-cells-col COLUMN          Cell-count QC column; default: n_cells.\n",
      "  --min-cells N                 Minimum cells per pseudobulk; default: 20.\n",
      "  --min-samples-per-group N     Minimum independent samples/group; default: 2.\n",
      "  --contrast EXPR               limma contrast, e.g. disease-control.\n",
      "  --help                        Print this help.\n\n",
      "The first column in both TSV files must contain row IDs:\n",
      "  counts: gene IDs\n",
      "  metadata: pseudobulk IDs, matching count-matrix column names\n\n"
    )
  )
}

parse_args <- function(args) {
  options <- list()
  i <- 1L

  while (i <= length(args)) {
    key <- args[[i]]

    if (key %in% c("--help", "-h")) {
      options$help <- TRUE
      i <- i + 1L
      next
    }

    if (!startsWith(key, "--")) {
      stop("Unexpected positional argument: ", key)
    }

    if (i == length(args)) {
      stop("Missing value for argument: ", key)
    }

    value <- args[[i + 1L]]
    name <- sub("^--", "", key)
    options[[name]] <- value
    i <- i + 2L
  }

  options
}

opt <- parse_args(commandArgs(trailingOnly = TRUE))

opt <- parse_args(commandArgs(trailingOnly = TRUE))

if (isTRUE(opt$help)) {
  usage()
  quit(status = 0L)
}

required <- c(
  "counts",
  "metadata",
  "output",
  "group-col",
  "group-levels"
)

missing_required <- required[
  !vapply(
    required,
    function(x) !is.null(opt[[x]]) && nzchar(opt[[x]]),
    logical(1)
  )
]

if (length(missing_required)) {
  usage()
  stop(
    "Missing required argument(s): ",
    paste(paste0("--", missing_required), collapse = ", ")
  )
}

counts_file <- normalizePath(opt[["counts"]], mustWork = TRUE)
metadata_file <- normalizePath(opt[["metadata"]], mustWork = TRUE)
output_file <- opt[["output"]]

group_col <- opt[["group-col"]]

group_levels <- trimws(
  strsplit(opt[["group-levels"]], ",", fixed = TRUE)[[1]]
)

if (length(group_levels) != 2L || any(!nzchar(group_levels))) {
  stop(
    "--group-levels must contain exactly two nonempty comma-separated values, ",
    "e.g. --group-levels control,disease"
  )
}

if (anyDuplicated(group_levels)) {
  stop("--group-levels must contain two distinct values.")
}

celltype_col <- opt[["celltype-col"]] %||% NULL
celltype <- opt[["celltype"]] %||% NULL

if (xor(is.null(celltype_col), is.null(celltype))) {
  stop(
    "Provide both --celltype-col and --celltype, or neither."
  )
}

covariate_cols <- opt[["covariates"]] %||% ""
covariate_cols <- trimws(
  strsplit(covariate_cols, ",", fixed = TRUE)[[1]]
)
covariate_cols <- covariate_cols[nzchar(covariate_cols)]

n_cells_col <- opt[["n-cells-col"]] %||% "n_cells"

min_cells <- suppressWarnings(
  as.integer(opt[["min-cells"]] %||% "20")
)

min_samples_per_group <- suppressWarnings(
  as.integer(opt[["min-samples-per-group"]] %||% "2")
)

if (is.na(min_cells) || min_cells < 1L) {
  stop("--min-cells must be a positive integer.")
}

if (is.na(min_samples_per_group) || min_samples_per_group < 2L) {
  stop("--min-samples-per-group must be an integer >= 2.")
}

contrast_str <- opt[["contrast"]] %||% NULL

counts_df <- utils::read.delim(
  counts_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

counts <- as.matrix(counts_df)
storage.mode(counts) <- "numeric"

if (!nrow(counts) || !ncol(counts)) {
  stop("Count matrix is empty.")
}

if (is.null(rownames(counts)) || any(!nzchar(rownames(counts)))) {
  stop("Counts must contain nonempty gene IDs in the first column.")
}

if (is.null(colnames(counts)) || any(!nzchar(colnames(counts)))) {
  stop("Counts must contain nonempty pseudobulk profile IDs in the header.")
}

if (anyDuplicated(rownames(counts))) {
  stop("Counts contains duplicated gene IDs.")
}

if (anyDuplicated(colnames(counts))) {
  stop("Counts contains duplicated pseudobulk IDs.")
}

if (anyNA(counts) || any(!is.finite(counts))) {
  stop("Counts contains missing or non-finite values.")
}

if (any(counts < 0)) {
  stop("Counts contains negative values; raw pseudobulk counts are required.")
}

if (any(abs(counts - round(counts)) > 1e-8)) {
  stop(
    "Counts are not integer-like. Use summed raw counts, not normalized values."
  )
}

meta <- utils::read.delim(
  metadata_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (is.null(rownames(meta)) || any(!nzchar(rownames(meta)))) {
  stop("Metadata must have nonempty pseudobulk IDs in its first column.")
}

if (anyDuplicated(rownames(meta))) {
  stop("Metadata contains duplicated pseudobulk IDs.")
}

missing_meta <- setdiff(colnames(counts), rownames(meta))

if (length(missing_meta)) {
  stop(
    "Count-matrix profiles absent from metadata:\n",
    paste(missing_meta, collapse = ", ")
  )
}

meta <- meta[colnames(counts), , drop = FALSE]

if (!identical(rownames(meta), colnames(counts))) {
  stop("Metadata row names are not aligned with count-matrix columns.")
}

required_metadata_columns <- c(
  group_col,
  n_cells_col,
  covariate_cols,
  if (!is.null(celltype_col)) celltype_col else character(0)
)

missing_metadata_columns <- setdiff(
  required_metadata_columns,
  colnames(meta)
)

if (length(missing_metadata_columns)) {
  stop(
    "Missing metadata column(s): ",
    paste(missing_metadata_columns, collapse = ", ")
  )
}

if (!is.null(celltype_col)) {
  keep_celltype <- as.character(meta[[celltype_col]]) == celltype

  if (!any(keep_celltype)) {
    stop(
      "No pseudobulk profiles found for ",
      celltype_col, " = '", celltype, "'."
    )
  }

  counts <- counts[, keep_celltype, drop = FALSE]
  meta <- meta[colnames(counts), , drop = FALSE]
}

keep_cells <- !is.na(meta[[n_cells_col]]) &
  suppressWarnings(as.numeric(meta[[n_cells_col]])) >= min_cells

if (!any(keep_cells)) {
  stop(
    "No pseudobulk profiles remain after filtering at min_cells = ",
    min_cells
  )
}

counts <- counts[, keep_cells, drop = FALSE]
meta <- meta[colnames(counts), , drop = FALSE]

group_values <- as.character(meta[[group_col]])

group <- factor(
  group_values,
  levels = group_levels
)

if (anyNA(group)) {
  invalid_groups <- unique(group_values[is.na(group)])

  stop(
    "Missing/unrecognized group values in '", group_col, "'.\n",
    "Expected: ", paste(group_levels, collapse = ", "), "\n",
    "Invalid: ", paste(invalid_groups, collapse = ", ")
  )
}

n_samples_by_group <- table(group)

if (any(n_samples_by_group < min_samples_per_group)) {
  stop(
    "Insufficient independent pseudobulk profiles per group:\n",
    paste(
      paste0(names(n_samples_by_group), "=", n_samples_by_group),
      collapse = ", "
    ),
    "\nRequired per group: ", min_samples_per_group
  )
}

design_data <- meta
design_data$group <- group

design_formula <- if (!length(covariate_cols)) {
  stats::as.formula("~ 0 + group")
} else {
  stats::as.formula(
    paste0("~ 0 + group + ", paste(covariate_cols, collapse = " + "))
  )
}

design <- stats::model.matrix(
  design_formula,
  data = design_data
)

colnames(design) <- sub("^group", "", colnames(design))
colnames(design) <- make.names(colnames(design))

if (qr(design)$rank < ncol(design)) {
  stop(
    "The DEA design matrix is not full rank. ",
    "Condition may be confounded with a covariate."
  )
}

y <- edgeR::DGEList(
  counts = round(counts),
  samples = meta
)

keep_genes <- edgeR::filterByExpr(
  y,
  design = design
)

if (!any(keep_genes)) {
  stop("No genes passed edgeR::filterByExpr().")
}

y <- y[keep_genes, , keep.lib.sizes = FALSE]
y <- edgeR::calcNormFactors(y, method = "TMM")

v <- limma::voom(
  y,
  design = design,
  plot = FALSE
)

fit <- limma::lmFit(
  v,
  design
)

if (is.null(contrast_str)) {
  group_levels_safe <- make.names(group_levels)
  contrast_str <- paste0(
    group_levels_safe[2L],
    "-",
    group_levels_safe[1L]
  )
}

contrast_matrix <- limma::makeContrasts(
  contrasts = contrast_str,
  levels = design
)

fit2 <- limma::contrasts.fit(
  fit,
  contrasts = contrast_matrix
)

fit2 <- limma::eBayes(
  fit2,
  robust = TRUE
)

results <- limma::topTable(
  fit2,
  coef = 1L,
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

results <- cbind(
  gene = rownames(results),
  results
)

results$contrast <- contrast_str
results$n_samples_reference <- unname(
  n_samples_by_group[group_levels[1L]]
)
results$n_samples_test <- unname(
  n_samples_by_group[group_levels[2L]]
)
results$n_genes_tested <- nrow(y)

output_dir <- dirname(output_file)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

utils::write.table(
  results,
  file = output_file,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

message("Completed limma-voom pseudobulk DEA.")
message("Contrast: ", contrast_str)
message(
  "Pseudobulk profiles: ",
  paste(
    paste0(names(n_samples_by_group), "=", n_samples_by_group),
    collapse = ", "
  )
)
message("Genes tested: ", nrow(y))
message("Results: ", normalizePath(output_file, mustWork = FALSE))