logging_tools <- tryCatch(get_tool_env("logging"), error = function(e) NULL)
if (is.null(logging_tools)) {
  log_info <- function(...) invisible(NULL)
  log_debug <- function(...) invisible(NULL)
  log_warn <- function(...) invisible(NULL)
  log_skip <- function(...) invisible(NULL)
  log_error <- function(...) invisible(NULL)
} else {
  log_info <- logging_tools$log_info
  log_debug <- logging_tools$log_debug
  log_warn <- logging_tools$log_warn
  log_skip <- if ("log_skip" %in% ls(envir = logging_tools)) logging_tools$log_skip else logging_tools$log_info
  log_error <- if ("log_error" %in% ls(envir = logging_tools)) logging_tools$log_error else logging_tools$log_warn
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

resolve_sample_order <- function(samples, config, utils_env) {
  config_extra <- tryCatch(config$extra, error = function(e) NULL)
  sample_order_cfg <- NULL
  if (!is.null(config_extra)) {
    sample_order_cfg <- tryCatch(config_extra$sample_order, error = function(e) NULL)
  }
  if (!is.null(sample_order_cfg) && length(sample_order_cfg) > 0) {
    ordered <- intersect(sample_order_cfg, samples)
    if (length(ordered) > 0) {
      return(ordered)
    }
  }
  if (!is.null(utils_env) && "infer_ordered_factor" %in% ls(envir = utils_env)) {
    inferred <- tryCatch(utils_env$infer_ordered_factor(samples), error = function(e) NULL)
    if (!is.null(inferred)) {
      return(levels(inferred))
    }
  }
  samples
}

run_qc_step <- function(label, expr) {
  tryCatch(
    {
      result <- expr()
      log_debug("QC step completed", context = list(step = label))
      if (is.null(result)) {
        return(TRUE)
      }
      result
    },
    error = function(err) {
      log_warn(
        "QC step failed",
        context = list(step = label, error = err$message)
      )
      FALSE
    }
  )
}

# Drop-in replacements for try(…, silent = TRUE)
safe_try <- function(expr) {
  tryCatch(expr, error = function(e) NULL)
}

safe_try_warn <- function(expr, warn_msg = NULL) {
  tryCatch(expr,
    error = function(e) {
      if (!is.null(warn_msg)) log_warn(warn_msg)
      NULL
    }
  )
}

plot_violin_distributions <- function(df_long, sample_order, qcdir) {
  df_long$sample <- factor(df_long$sample, levels = sample_order)
  p_violin <- ggplot2::ggplot(df_long, ggplot2::aes(x = sample, y = value)) +
    ggplot2::geom_violin(fill = "#a6cee3", color = NA) +
    ggplot2::geom_boxplot(width = 0.1, outlier.size = 0.2) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::labs(title = "Sample distributions (normalized)", x = "Sample", y = "Intensity")
  ggplot2::ggsave(file.path(qcdir, "violin_distributions.pdf"), p_violin, width = 8.5, height = 11)
}

plot_nonzero_counts <- function(mat_num, sample_order, qcdir) {
  nz <- colSums(mat_num > 0, na.rm = TRUE)
  df_nz <- data.frame(sample = names(nz), nonzero = as.integer(nz))
  df_nz$sample <- factor(df_nz$sample, levels = sample_order)
  p_nz <- ggplot2::ggplot(df_nz, ggplot2::aes(x = sample, y = nonzero)) +
    ggplot2::geom_col(fill = "#b2df8a") +
    ggplot2::coord_flip() +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::labs(title = "Non-zero site counts per sample", x = "Sample", y = "# nonzero")
  ggplot2::ggsave(file.path(qcdir, "nonzero_counts.pdf"), p_nz, width = 8.5, height = 11)
}

plot_sty_overall <- function(gct, qcdir) {
  if (!"AA" %in% colnames(gct@rdesc)) {
    return(FALSE)
  }
  aa_df <- data.frame(AA = unique(gct@rdesc$AA), stringsAsFactors = FALSE)
  aa_df$AA <- as.character(aa_df$AA)
  aa_df <- aa_df[!is.na(aa_df$AA), , drop = FALSE]
  aa_df <- aa_df[aa_df$AA %in% c("S", "T", "Y"), , drop = FALSE]
  if (nrow(aa_df) == 0) {
    return(FALSE)
  }

  sty_counts <- as.data.frame(table(aa_df$AA), stringsAsFactors = FALSE)
  colnames(sty_counts) <- c("AA", "count")
  sty_counts$AA <- factor(sty_counts$AA, levels = c("Y", "T", "S"), ordered = TRUE)
  sty_counts <- dplyr::arrange(sty_counts, AA)

  p_sty <- ggplot2::ggplot(sty_counts, ggplot2::aes(x = AA, y = count, fill = AA)) +
    ggplot2::geom_col(width = 0.65, show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = c(S = "#8da0cb", T = "#66c2a5", Y = "#fc8d62")) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::labs(title = "Site composition (S/T/Y)", x = "Residue", y = "# sites")
  ggplot2::ggsave(file.path(qcdir, "sty_counts_total.pdf"), p_sty, width = 5.6, height = 4.2)

  y_only <- subset(sty_counts, AA == "Y")
  if (nrow(y_only) == 1) {
    p_y <- ggplot2::ggplot(y_only, ggplot2::aes(x = AA, y = count)) +
      ggplot2::geom_col(fill = "#fc8d62", width = 0.5, show.legend = FALSE) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::scale_x_discrete(labels = c("Y")) +
      ggplot2::labs(title = "Tyrosine-only site count", x = "Residue", y = "# sites")
    ggplot2::ggsave(file.path(qcdir, "sty_counts_Y_only.pdf"), p_y, width = 3.6, height = 4.0)
  }
  TRUE
}

plot_sty_by_sample <- function(gct, mat_num, all_samples, sample_order, qcdir) {
  present_long <- as.data.frame(mat_num)
  present_long$id <- rownames(mat_num)
  present_long <- tidyr::pivot_longer(present_long, cols = -id, names_to = "sample", values_to = "value")
  present_long <- dplyr::filter(present_long, is.finite(value), value > 0)
  aa_map <- data.frame(id = rownames(gct@rdesc), AA = as.character(gct@rdesc$AA), stringsAsFactors = FALSE)
  present_long <- dplyr::left_join(present_long, aa_map, by = "id")
  present_long <- dplyr::filter(present_long, AA %in% c("S", "T", "Y"))
  if (nrow(present_long) == 0) {
    return(FALSE)
  }

  sty_by_sample <- present_long %>%
    dplyr::group_by(sample, AA) %>%
    dplyr::summarize(count = dplyr::n(), .groups = "drop")

  sty_by_sample$AA <- factor(sty_by_sample$AA, levels = c("S", "T", "Y"))
  p_sty_by <- ggplot2::ggplot(
    sty_by_sample,
    ggplot2::aes(y = factor(sample, levels = sample_order), x = count, fill = AA)
  ) +
    ggplot2::geom_col(width = 0.7, position = "stack") +
    ggplot2::scale_fill_manual(values = c(S = "#8da0cb", T = "#66c2a5", Y = "#fc8d62")) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::labs(title = "S/T/Y site counts per sample (non-zero)", x = "# sites", y = "Sample")
  ggplot2::ggsave(file.path(qcdir, "sty_counts_by_sample_STY.pdf"), p_sty_by, width = 7.8, height = 0.3 + 0.28 * length(sample_order))

  y_by_sample <- sty_by_sample %>%
    dplyr::filter(AA == "Y") %>%
    tidyr::complete(sample = all_samples, fill = list(count = 0))
  y_by_sample$sample <- factor(y_by_sample$sample, levels = sample_order)
  p_y_by <- ggplot2::ggplot(y_by_sample, ggplot2::aes(y = sample, x = count)) +
    ggplot2::geom_col(fill = "#fc8d62", width = 0.7) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::labs(title = "Y site counts per sample (non-zero)", x = "# sites (Y)", y = "Sample")
  ggplot2::ggsave(file.path(qcdir, "sty_counts_by_sample_Y.pdf"), p_y_by, width = 6.8, height = 0.3 + 0.28 * length(all_samples))
  TRUE
}

compute_missingness_matrix <- function(mat_num, value_threshold = NULL) {
  missing_matrix <- !is.finite(mat_num)
  if (!is.null(value_threshold) && is.finite(value_threshold)) {
    missing_matrix <- missing_matrix | (mat_num <= value_threshold)
  }
  storage.mode(missing_matrix) <- "logical"
  missing_matrix
}

compute_missingness_stats <- function(missing_matrix) {
  total_sites <- nrow(missing_matrix)
  if (total_sites == 0) {
    return(data.frame(
      sample = character(),
      missing = integer(),
      present = integer(),
      total = integer(),
      fraction = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  counts <- colSums(missing_matrix, na.rm = TRUE)
  data.frame(
    sample = names(counts),
    missing = as.integer(counts),
    present = as.integer(total_sites - counts),
    total = as.integer(total_sites),
    fraction = if (total_sites > 0) counts / total_sites else 0,
    stringsAsFactors = FALSE
  )
}

plot_missingness_counts <- function(missing_stats, sample_order, qcdir) {
  if (nrow(missing_stats) == 0) {
    log_warn("Missingness statistics empty; skipping counts plot")
    return(FALSE)
  }
  missing_stats$sample <- factor(missing_stats$sample, levels = sample_order)
  p_missing_counts <- ggplot2::ggplot(missing_stats, ggplot2::aes(x = sample, y = missing)) +
    ggplot2::geom_col(fill = "#fdbf6f") +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.1f%%", fraction * 100)),
      hjust = -0.1,
      size = 2.7
    ) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::labs(
      title = "Missing values per sample",
      x = "Sample",
      y = "# missing (≤ threshold or non-finite)",
      caption = "Percent labels reflect missing / total sites"
    )
  ggplot2::ggsave(file.path(qcdir, "missing_counts.pdf"), p_missing_counts, width = 8.5, height = 11)
  TRUE
}

plot_missingness_fraction <- function(missing_stats, sample_order, qcdir) {
  if (nrow(missing_stats) == 0) {
    log_warn("Missingness statistics empty; skipping fraction plot")
    return(FALSE)
  }
  missing_stats$sample <- factor(missing_stats$sample, levels = sample_order)
  p_missing_fraction <- ggplot2::ggplot(missing_stats, ggplot2::aes(x = sample, y = fraction)) +
    ggplot2::geom_col(fill = "#fb9a99") +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.1f%%", fraction * 100)),
      hjust = -0.1,
      size = 2.7
    ) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      expand = ggplot2::expansion(mult = c(0, 0.05)),
      labels = function(x) sprintf("%.0f%%", x * 100)
    ) +
    ggplot2::labs(
      title = "Missingness fraction per sample",
      x = "Sample",
      y = "% missing (≤ threshold or non-finite)"
    )
  ggplot2::ggsave(file.path(qcdir, "missing_fraction.pdf"), p_missing_fraction, width = 8.5, height = 11)
  TRUE
}

# Sample correlation heatmap (Pearson by default, Spearman optional)
plot_sample_correlations <- function(mat_num, sample_order, qcdir, method = "pearson", use_zscore = FALSE) {
  if (ncol(mat_num) < 2) {
    log_warn("Not enough samples for correlation heatmap")
    return(FALSE)
  }
  # replace non-finite with NA to avoid stopping cor()
  mat_num[!is.finite(mat_num)] <- NA_real_
  if (isTRUE(use_zscore)) {
    mat_num <- t(scale(t(mat_num)))
  }
  cor_mat <- tryCatch(stats::cor(mat_num, use = "pairwise.complete.obs", method = method), error = function(e) NULL)
  if (is.null(cor_mat)) {
    log_warn("Failed to compute sample correlations")
    return(FALSE)
  }
  cor_df <- reshape2::melt(cor_mat, varnames = c("sample1", "sample2"), value.name = "cor")
  cor_df$sample1 <- factor(cor_df$sample1, levels = sample_order)
  cor_df$sample2 <- factor(cor_df$sample2, levels = sample_order)
  p_cor <- ggplot2::ggplot(cor_df, ggplot2::aes(x = sample1, y = sample2, fill = cor)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_distiller(palette = "RdYlBu", limits = c(-1, 1), oob = scales::squish, direction = 1) +
    ggplot2::theme_bw(base_size = 8) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ggplot2::labs(title = sprintf("Sample correlations (%s)", method), x = "", y = "")
  ggplot2::ggsave(file.path(qcdir, sprintf("cor_heatmap_%s.pdf", method)), p_cor, width = 8.0, height = 8.0)
  TRUE
}

# PCA of samples (on log2-transformed data)
plot_sample_pca <- function(mat_num, sample_order, qcdir, topn = 5000, cdesc = NULL, log_transformed = FALSE) {
  if (ncol(mat_num) < 2) {
    log_warn("Not enough samples for PCA")
    return(FALSE)
  }
  mat_num[!is.finite(mat_num)] <- NA_real_
  # log2 transform with small pseudo-count (skip if already log-transformed upstream)
  if (isTRUE(log_transformed)) {
    mat_log <- mat_num
  } else {
    mat_log <- log2(mat_num + 1e-3)
  }
  # keep top variable features
  vars <- matrixStats::rowVars(mat_log, na.rm = TRUE)
  ord <- order(vars, decreasing = TRUE)
  keep <- ord[seq_len(min(length(ord), topn))]
  mat_use <- mat_log[keep, , drop = FALSE]
  # impute missing for PCA with column medians
  col_meds <- matrixStats::colMedians(mat_use, na.rm = TRUE)
  for (i in seq_len(ncol(mat_use))) {
    nz <- mat_use[, i]
    nz[!is.finite(nz)] <- col_meds[i]
    mat_use[, i] <- nz
  }
  pca <- tryCatch(stats::prcomp(t(mat_use), scale. = TRUE), error = function(e) NULL)
  if (is.null(pca)) {
    log_warn("PCA failed")
    return(FALSE)
  }
  pc_df <- as.data.frame(pca$x[, 1:3, drop = FALSE])
  pc_df$sample <- colnames(mat_num)
  pc_df$sample <- factor(pc_df$sample, levels = sample_order)

  # optional coloring by cdesc fields
  color_fields <- NULL
  if (!is.null(cdesc)) {
    color_fields <- intersect(c("group", "batch", "condition"), colnames(cdesc))
    if (length(color_fields) == 0) color_fields <- NULL
  }

  plots <- list()
  var_exp <- (pca$sdev^2) / sum(pca$sdev^2)
  title_suffix <- sprintf(" (top %d features)", nrow(mat_use))
  for (pair in list(c(1,2), c(1,3), c(2,3))) {
    ax <- pair[1]; ay <- pair[2]
    p <- ggplot2::ggplot(pc_df, ggplot2::aes_string(x = paste0("PC", ax), y = paste0("PC", ay))) +
      ggplot2::geom_point(size = 2, alpha = 0.8) +
      ggplot2::theme_bw(base_size = 9) +
      ggplot2::labs(
        title = sprintf("PCA PC%d vs PC%d%s", ax, ay, title_suffix),
        x = sprintf("PC%d (%.1f%%)", ax, var_exp[ax]*100),
        y = sprintf("PC%d (%.1f%%)", ay, var_exp[ay]*100)
      )
    if (!is.null(color_fields)) {
      f <- color_fields[1]
      if (f %in% colnames(cdesc)) {
        pc_df[[f]] <- cdesc[match(pc_df$sample, rownames(cdesc)), f]
        p <- p + ggplot2::geom_point(ggplot2::aes_string(color = f), size = 2, alpha = 0.8) +
          ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1)))
      }
    }
    plots[[length(plots)+1]] <- p
  }
  for (i in seq_along(plots)) {
    ggplot2::ggsave(file.path(qcdir, sprintf("pca_%d.pdf", i)), plots[[i]], width = 6.0, height = 5.0)
  }
  TRUE
}
plot_missingness_upset <- function(
    missing_matrix,
    sample_order,
    qcdir,
    top_n = 8,
    min_intersections = 2,
    max_features = 5000) {
  if (!requireNamespace("ComplexUpset", quietly = TRUE)) {
    log_warn("ComplexUpset not available; skipping missingness UpSet plot")
    return(FALSE)
  }

  if (is.null(missing_matrix) || ncol(missing_matrix) < 2) {
    log_warn("Not enough samples for missingness UpSet plot")
    return(FALSE)
  }

  missing_counts <- colSums(missing_matrix)
  keep_samples <- names(missing_counts[missing_counts > 0])
  if (length(keep_samples) < 2) {
    log_warn("Insufficient missingness across samples for UpSet plot")
    return(FALSE)
  }

  ordered_by_missing <- names(sort(missing_counts[keep_samples], decreasing = TRUE))
  top_candidates <- head(ordered_by_missing, min(top_n, length(ordered_by_missing)))
  if (!is.null(sample_order)) {
    top_samples <- intersect(sample_order, top_candidates)
  } else {
    top_samples <- top_candidates
  }
  if (length(top_samples) < 2) {
    log_warn("Fewer than two samples selected for missingness UpSet plot")
    return(FALSE)
  }

  missing_top <- missing_matrix[, top_samples, drop = FALSE]
  row_counts <- rowSums(missing_top)
  missing_top <- missing_top[row_counts >= min_intersections, , drop = FALSE]
  if (nrow(missing_top) == 0) {
    log_warn("No features meet missingness intersection threshold for UpSet plot")
    return(FALSE)
  }

  if (!is.null(max_features) && is.finite(max_features) && nrow(missing_top) > max_features) {
    ord <- order(row_counts[row_counts >= min_intersections], decreasing = TRUE)
    missing_top <- missing_top[ord, , drop = FALSE]
    missing_top <- missing_top[seq_len(min(max_features, nrow(missing_top))), , drop = FALSE]
  }

  missing_df <- as.data.frame(missing_top)
  missing_df[] <- lapply(missing_df, as.integer)
  missing_df$feature_id <- rownames(missing_top)

  upset_plot <- ComplexUpset::upset(
    missing_df,
    intersect = colnames(missing_top),
    base_annotations = list(
      "Intersection size" = ComplexUpset::intersection_size(counts = TRUE) +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.1)))
    ),
    set_sizes = ComplexUpset::set_size(counts = TRUE)
  )

  # Rotate set labels for readability when possible
  try({
    upset_plot[[3]] <- upset_plot[[3]] +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
  }, silent = TRUE)

  ggplot2::ggsave(file.path(qcdir, "missing_upset.pdf"), upset_plot, width = 8.5, height = 6.2)
  TRUE
}

generate_qc_reports <- function(gct, config, output_dir) {
  if (is.null(gct) || is.null(gct@mat)) {
    stop("GCT input for QC is invalid.")
  }
  qcdir <- ensure_dir(file.path(output_dir, "qc"))
  mat_num <- gct@mat
  all_samples <- colnames(mat_num)
  utils_env <- tryCatch(get_tool_env("utils"), error = function(e) NULL)
  sample_order_levels <- resolve_sample_order(all_samples, config, utils_env)
  downsample_n <- tryCatch(config$qc$downsample_n, error = function(e) NULL)
  if (is.null(downsample_n)) {
    downsample_n <- 200000
  }
  qc_level <- tryCatch(config$qc$level, error = function(e) NULL)
  qc_level <- suppressWarnings(as.integer(qc_level)[1])
  if (!is.finite(qc_level) || qc_level < 1L) {
    qc_level <- 1L
  }

  missing_cfg <- tryCatch(config$qc$missingness, error = function(e) NULL)
  log_transform_flag <- tryCatch(config$norm$log_transform, error = function(e) NULL)
  value_threshold <- NULL
  if (!is.null(missing_cfg) && !is.null(missing_cfg$value_threshold)) {
    value_threshold <- suppressWarnings(as.numeric(missing_cfg$value_threshold))[1]
    if (!is.finite(value_threshold)) {
      value_threshold <- NULL
    }
  } else if (!is.null(log_transform_flag) && isTRUE(log_transform_flag)) {
    value_threshold <- NULL
  } else {
    value_threshold <- 0
  }
  missing_matrix <- compute_missingness_matrix(mat_num, value_threshold = value_threshold)
  missing_stats <- compute_missingness_stats(missing_matrix)
  log_debug(
    "Missingness statistics prepared",
    context = list(
      threshold = if (is.null(value_threshold)) "non-finite" else value_threshold,
      log_transform = isTRUE(log_transform_flag),
      samples = ncol(missing_matrix),
      features = nrow(missing_matrix)
    )
  )
  upset_top_n <- if (!is.null(missing_cfg) && !is.null(missing_cfg$upset_top)) {
    as.integer(missing_cfg$upset_top)
  } else {
    8L
  }
  if (!is.finite(upset_top_n) || upset_top_n < 2L) {
    upset_top_n <- 8L
  }
  upset_min_intersections <- if (!is.null(missing_cfg) && !is.null(missing_cfg$upset_min_intersections)) {
    as.integer(missing_cfg$upset_min_intersections)
  } else {
    2L
  }
  if (!is.finite(upset_min_intersections) || upset_min_intersections < 1L) {
    upset_min_intersections <- 2L
  }
  upset_max_features <- if (!is.null(missing_cfg) && !is.null(missing_cfg$upset_max_features)) {
    as.integer(missing_cfg$upset_max_features)
  } else {
    5000L
  }
  if (!is.finite(upset_max_features) || upset_max_features <= 0L) {
    upset_max_features <- 5000L
  }

  log_info(
    "Running QC plot generation",
    context = list(samples = length(all_samples), qcdir = qcdir)
  )

  df_long <- as.data.frame(mat_num)
  df_long$id <- rownames(mat_num)
  df_long <- tidyr::pivot_longer(df_long, cols = -id, names_to = "sample", values_to = "value")
  if (nrow(df_long) > downsample_n) {
    df_long <- df_long[sample(nrow(df_long), downsample_n), , drop = FALSE]
  }
  df_long$sample <- as.character(df_long$sample)

  run_qc_step("violin_distributions", function() plot_violin_distributions(df_long, sample_order_levels, qcdir))
  run_qc_step("nonzero_counts", function() plot_nonzero_counts(mat_num, sample_order_levels, qcdir))
  run_qc_step("missing_counts", function() plot_missingness_counts(missing_stats, sample_order_levels, qcdir))
  run_qc_step("missing_fraction", function() plot_missingness_fraction(missing_stats, sample_order_levels, qcdir))

  cor_do <- tryCatch(config$qc$cor_do, error = function(e) NULL)
  if (is.null(cor_do)) cor_do <- (qc_level >= 2L)
  if (isTRUE(cor_do)) {
    cor_method <- tryCatch(config$qc$cor_method, error = function(e) NULL)
    if (is.null(cor_method)) cor_method <- "pearson"
    cor_use_z <- tryCatch(config$qc$cor_use_zscore, error = function(e) NULL)
    if (is.null(cor_use_z)) cor_use_z <- FALSE
    run_qc_step(
      "cor_heatmap",
      function() plot_sample_correlations(mat_num, sample_order_levels, qcdir, method = cor_method, use_zscore = cor_use_z)
    )
  } else {
    log_skip("QC correlation heatmap disabled", context = list(level = qc_level))
  }

  # QC PCA is optional; prefer the main [params.pca] routine for "level 2" QC.
  pca_do <- tryCatch(config$qc$pca_do, error = function(e) NULL)
  if (is.null(pca_do)) {
    pca_do <- FALSE
  }
  if (isTRUE(pca_do) && qc_level >= 2L) {
    pca_topn <- tryCatch(config$qc$pca_topn, error = function(e) NULL)
    if (is.null(pca_topn)) pca_topn <- 5000
    run_qc_step(
      "pca",
      function() plot_sample_pca(mat_num, sample_order_levels, qcdir, topn = pca_topn, cdesc = gct@cdesc, log_transformed = isTRUE(log_transform_flag))
    )
  } else {
    log_skip("QC PCA disabled", context = list(level = qc_level))
  }

  sty_any <- run_qc_step("sty_overall", function() plot_sty_overall(gct, qcdir))
  if (isTRUE(sty_any)) {
    run_qc_step("sty_by_sample", function() plot_sty_by_sample(gct, mat_num, all_samples, sample_order_levels, qcdir))
  }

  upset_do <- tryCatch(config$qc$missingness$upset_do, error = function(e) NULL)
  if (is.null(upset_do)) upset_do <- (qc_level >= 2L)
  if (isTRUE(upset_do)) {
    run_qc_step(
      "missing_upset",
      function() plot_missingness_upset(
        missing_matrix,
        sample_order_levels,
        qcdir,
        top_n = upset_top_n,
        min_intersections = upset_min_intersections,
        max_features = upset_max_features
      )
    )
  } else {
    log_skip("QC missingness UpSet disabled", context = list(level = qc_level))
  }

  log_info("QC plot generation finished", context = list(qcdir = qcdir))
  invisible(qcdir)
}
