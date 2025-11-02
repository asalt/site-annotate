suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("--matrix"), type = "character", help = "Path to TSV matrix"),
  make_option(c("--output"), type = "character", help = "Output PDF path"),
  make_option(c("--title"), type = "character", default = "Gene Heatmap"),
  make_option(c("--metric"), type = "character", default = "log2ratio"),
  make_option(c("--cluster-rows"), type = "logical", default = TRUE),
  make_option(c("--cluster-cols"), type = "logical", default = TRUE),
  make_option(c("--zscore"), type = "logical", default = FALSE),
  make_option(c("--scale-min"), type = "double", default = NA_real_),
  make_option(c("--scale-max"), type = "double", default = NA_real_)
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$matrix) || is.null(opt$output)) {
  stop("--matrix and --output are required")
}

if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  stop("ComplexHeatmap package is required to render heatmaps")
}
if (!requireNamespace("circlize", quietly = TRUE)) {
  stop("circlize package is required by ComplexHeatmap")
}

suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))

mat <- read.delim(opt$matrix, header = TRUE, check.names = FALSE, row.names = 1, sep = "\t")
if (nrow(mat) == 0) {
  stop("Matrix has no rows")
}

num_df <- lapply(mat, function(x) suppressWarnings(as.numeric(as.character(x))))
num_df <- as.data.frame(num_df, check.names = FALSE, stringsAsFactors = FALSE)
mat_num <- as.matrix(num_df)
rownames(mat_num) <- rownames(mat)
colnames(mat_num) <- colnames(mat)

# Optional row-wise z-score
if (isTRUE(opt$zscore)) {
  # compute row-wise z-score, guarding 0 sd
  row_means <- rowMeans(mat_num, na.rm = TRUE)
  row_sds <- apply(mat_num, 1, sd, na.rm = TRUE)
  row_sds[is.na(row_sds) | row_sds == 0] <- 1
  mat_num <- sweep(mat_num, 1, row_means, FUN = "-")
  mat_num <- sweep(mat_num, 1, row_sds, FUN = "/")
}

# Drop rows/columns that are entirely NA (post z-score if applied)
if (nrow(mat_num) > 0) {
  keep_r <- apply(mat_num, 1, function(x) !all(is.na(x)))
  mat_num <- mat_num[keep_r, , drop = FALSE]
}
if (ncol(mat_num) > 0) {
  keep_c <- apply(mat_num, 2, function(x) !all(is.na(x)))
  mat_num <- mat_num[, keep_c, drop = FALSE]
}

if (nrow(mat_num) == 0 || ncol(mat_num) == 0) {
  stop("Matrix is empty after numeric conversion (all-NA rows/cols)")
}

# Replace underscores with spaces in labels for readability
colnames(mat_num) <- gsub("_", " ", colnames(mat_num))
rownames(mat_num) <- gsub("_", " ", rownames(mat_num))

width <- 6 + pmin(ncol(mat_num) * 0.6, 8)
height <- 4 + pmin(nrow(mat_num) * 0.25, 6)

pdf(opt$output, width = width, height = height)
cl_title <- paste0(
  opt$title,
  " (rc", ifelse(isTRUE(opt$`cluster-rows`), "T", "F"),
  "_cc", ifelse(isTRUE(opt$`cluster-cols`), "T", "F"), ")"
)
if (isTRUE(opt$zscore)) {
  cl_title <- paste0(cl_title, " [z]")
}

legend_name <- if (isTRUE(opt$zscore)) paste0(opt$metric, " (z)") else opt$metric

# Build color function using global dataset scale if provided
col_fun <- NULL
if (!is.na(opt$`scale-min`) && !is.na(opt$`scale-max`)) {
  smin <- opt$`scale-min`
  smax <- opt$`scale-max`
  if (smin == smax) {
    smin <- smin - 1
    smax <- smax + 1
  }
  if (smin < 0 && smax > 0) {
    col_fun <- circlize::colorRamp2(c(smin, 0, smax), c("#2c7bb6", "#ffffbf", "#d7191c"))
  } else {
    col_fun <- circlize::colorRamp2(c(smin, smax), c("#2c7bb6", "#d7191c"))
  }
}

ht <- ComplexHeatmap::Heatmap(
  mat_num,
  name = legend_name,
  column_title = cl_title,
  cluster_rows = isTRUE(opt$`cluster-rows`),
  cluster_columns = isTRUE(opt$`cluster-cols`),
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_names_rot = 45,
  heatmap_legend_param = list(direction = "horizontal"),
  col = col_fun
)
ComplexHeatmap::draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()
