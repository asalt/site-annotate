# run.R
# main entry point for all downstream analyses

source(file.path("setup_environment.R"))
source(file.path("lazyloader.R"))

log_tools <- tryCatch(get_tool_env("logging"), error = function(e) NULL)
if (is.null(log_tools)) {
  format_ctx <- function(ctx) {
    if (is.null(ctx) || length(ctx) == 0) {
      return("")
    }
    if (is.null(names(ctx))) {
      names(ctx) <- rep("", length(ctx))
    }
    pieces <- Map(function(name, value) {
      if (is.list(value)) {
        value <- unlist(value, recursive = TRUE, use.names = FALSE)
      }
      if (length(value) > 1) {
        value <- paste(value, collapse = ",")
      }
      value <- encodeString(as.character(value), quote = "")
      if (nzchar(name)) {
        sprintf("%s=%s", name, value)
      } else {
        value
      }
    }, names(ctx), ctx)
    paste(unlist(pieces), collapse = " ")
  }
  base_logger <- function(level) {
    function(msg, ..., context = NULL) {
      if (length(list(...)) > 0) {
        msg <- sprintf(msg, ...)
      }
      ctx <- format_ctx(context)
      prefix <- sprintf("[site-annotate][%s]", level)
      message(trimws(paste(prefix, msg, ctx)))
    }
  }
  log_info <- base_logger("INFO")
  log_debug <- base_logger("DEBUG")
  log_warn <- base_logger("WARN")
  log_skip <- base_logger("SKIP")
  log_error <- base_logger("ERROR")
} else {
  log_info <- log_tools$log_info
  log_debug <- log_tools$log_debug
  log_warn <- log_tools$log_warn
  log_skip <- log_tools$log_skip
  log_error <- log_tools$log_error
}

sanitize_path_component <- function(name) {
  if (is.null(name) || !nzchar(name)) {
    return("unnamed")
  }
  name <- trimws(as.character(name))
  name <- gsub("[[:space:]]+", "_", name)
  name <- gsub("[^A-Za-z0-9_.-]", "_", name)
  name <- gsub("_+", "_", name)
  name <- gsub("^[._]+", "", name)
  name <- gsub("[._]+$", "", name)
  if (!nzchar(name)) name <- "unnamed"
  name
}

cache_tools <- get_tool_env("cache")
`%cached_by%` <- cache_tools$`%cached_by%`



process_gct_file <- function(gct_file, config) {
  util_tools <- get_tool_env("utils")
  log_info("Processing GCT file", context = list(file = gct_file))

  if (!file.exists(gct_file)) {
    log_error("GCT file missing", context = list(file = gct_file))
    stop(paste0("File ", gct_file, " does not exist, exiting"))
  }

  # Step 1: Load or Parse GCT File
  # gct <- load_or_parse_gct(gct_file)
  log_debug("Parsing GCT file", context = list(file = gct_file))
  gct <- parse_gctx(gct_file) %cached_by% rlang::hash_file(gct_file)

  # if ("species" %in% names(config)){
  # }

  # Step 2: Normalize GCT
  norm_config <- config$norm
  # Support both historical 'mednorm' and current 'normalize' flags; default to TRUE for median-normalization
  mednorm_flag <- if (!is.null(norm_config$normalize)) norm_config$normalize else (norm_config$mednorm %||% TRUE)
  # Log-transform is opt-in; default FALSE unless explicitly enabled in config
  log_flag <- norm_config$log_transform %||% FALSE
  log_debug(
    "Normalizing matrix",
    context = list(
      mednorm = mednorm_flag,
      log_transform = log_flag
    )
  )
  normed_gct <- util_tools$normalize_gct(gct,
    mednorm = mednorm_flag,
    log_transform = log_flag
  ) %cached_by%
    rlang::hash(c(gct@mat, mednorm_flag, log_flag))

  # Step 3: Filter Non-Zeros
  orig_dim <- dim(normed_gct@mat)
  filtered_gct <- util_tools$filter_nonzeros(
    normed_gct, config$filter$non_zeros %||% 1.0,
    config$filter$nonzero_subgroup %||% NULL
  ) %cached_by% rlang::hash(c(normed_gct@mat, config$filter$non_zeros %||% 1.0, config$filter$nonzero_subgroup %||% NULL))

  result_dim <- dim(filtered_gct@mat)
  log_info(
    "Filtered zero-only rows",
    context = list(
      before = paste0(as.character(orig_dim), collapse = "x"),
      after = paste0(as.character(result_dim), collapse = "x"),
      threshold = config$filter$non_zeros %||% 1.0
    )
  )


  # Step 4: Apply Batch Correction (if required)
  # need more error checking ehre or inside
  final_gct <- filtered_gct
  if (is.character(config$norm$batch) %||% FALSE) {
    log_info("Applying ComBat batch correction", context = list(batch = config$norm$batch))
    final_gct <- util_tools$do_combat(filtered_gct, by = config$norm$batch) %cached_by% rlang::hash(c(filtered_gct@mat, config$norm$batch))
    # final_gct <- util_tools$do_limmaremovebatch(filtered_gct, by = config$norm$batch) %cached_by% rlang::hash(c(filtered_gct@mat, config$norm$batch))
  }

  return(final_gct)
}


exclude_samples <- function(gct, sample_exclude) {
  if (!is.null(sample_exclude)) {
    to_exclude <- intersect(rownames(gct@cdesc), sample_exclude)
    if (length(to_exclude) > 0) {
      log_info("Excluding samples", context = list(count = length(to_exclude)))
    } else {
      log_debug("No matching samples found for exclusion", context = list(requested = length(sample_exclude)))
    }
    to_keep <- setdiff(rownames(gct@cdesc), to_exclude)
    gct <- gct %>% subset_gct(cid = to_keep)
  }
  return(gct)
}

load_and_validate_config <- function(config_file, data_dir) {
  io_tools <- get_tool_env("io")
  log_info("Loading configuration", context = list(file = config_file, data_dir = data_dir %||% "<NA>"))
  config <- io_tools$load_config(config_file, root_dir = data_dir)
  if (is.null(config)) {
    log_error("Configuration load failed", context = list(file = config_file))
    stop("Config file could not be loaded.")
  }
  return(config)
}


# Main function to handle heatmap generation
generate_global_heatmaps <- function(gct_z, config, outdir, friendly_name, outname = NULL) {
  param_grid <- expand.grid(
    cut_by = unique(c(config$heatmap$cut_by %||% NA, NA)),
    cluster_columns = c(TRUE, FALSE),
    stringsAsFactors = FALSE
  )
  log_info(
    "Generating global heatmaps",
    context = list(
      combinations = nrow(param_grid),
      friendly_name = friendly_name,
      outdir = outdir
    )
  )

  # Step 2: Iterate over parameter grid and generate heatmaps
  param_grid %>%
    purrr::pmap(~ generate_global_heatmap(
      gct_z = gct_z,
      cut_by = ..1,
      cluster_columns = ..2,
      config = config,
      outdir = outdir,
      friendly_name = friendly_name,
      outname = outname
    ))
}

# Helper: Generate a single heatmap
generate_global_heatmap <- function(gct_z, cut_by, cluster_columns, config, outdir, friendly_name, outname = NULL) {
  # Prepare the subset GCT
  if (!"fifteenmer" %in% colnames(gct_z@rdesc)) {
    gct_z@rdesc$fifteenmer <- ""
  }
  rids_to_keep <- gct_z@rdesc %>%
    distinct(fifteenmer, .keep_all = TRUE) %>%
    pull(id)

  gct_z_subset <- gct_z %>% subset_gct(rid = rids_to_keep)
  log_debug(
    "Prepared global heatmap subset",
    context = list(
      cut_by = ifelse(is.na(cut_by), "<NA>", as.character(cut_by)),
      cluster_columns = cluster_columns,
      rows = nrow(gct_z_subset@mat),
      cols = ncol(gct_z_subset@mat),
      friendly_name = friendly_name
    )
  )

  # Generate heatmap code
  heatmap_tools <- get_tool_env("heatmap")
  ht_draw_code <- gct_z_subset %>% heatmap_tools$make_heatmap_fromgct(
    show_row_names = FALSE,
    optimal_size = FALSE,
    cut_by = cut_by,
    cluster_columns = cluster_columns,
    include_row_annotations = FALSE,
    meta_to_include = config$heatmap$legend_include %||% NULL,
    meta_to_exclude = config$heatmap$legend_exclude %||% NULL,
    column_title = friendly_name
  )
  ht_dims <- attr(ht_draw_code, "heatmap_dims")
  n_rows_plot <- if (!is.null(ht_dims) && !is.null(ht_dims$nrow)) ht_dims$nrow else nrow(gct_z_subset@mat)
  n_cols_plot <- if (!is.null(ht_dims) && !is.null(ht_dims$ncol)) ht_dims$ncol else ncol(gct_z_subset@mat)

  # Prepare file paths
  outf <- heatmap_tools$prepare_heatmap_output_path(outdir,
    "heatmap",
    outname,
    prefix = glue::glue("X{n_rows_plot}x{n_cols_plot}"),
    # prefix = list(nrow=nrow(gct_z_subset@mat), ncol(gct_z_subset@mat)),
    cut_by = cut_by,
    cluster_rows = T,
    cluster_columns = cluster_columns
  )

  # Skip if file exists and replace is FALSE
  if (file.exists(outf) && config$advanced$replace == FALSE) {
    log_skip(
      "Global heatmap already exists",
      context = list(outf = outf, replace = config$advanced$replace %||% FALSE)
    )
    return(NULL)
  }

  # Draw and save the heatmap
  log_info("Rendering global heatmap", context = list(outf = outf))
  heatmap_tools$draw_and_save_heatmap(ht_draw_code, outf, gct_z_subset, width = NULL, height = 11)
  log_info("Saved global heatmap", context = list(outf = outf))
}

generate_selection_heatmaps <- function(gct_z, config, outdir, friendly_name = "", outname = "") {
  io_tools <- get_tool_env("io")
  heatmap_tools <- get_tool_env("heatmap")


  param_grid <- expand.grid(
    cut_by = config$heatmap$selection$cut_by %||% config$heatmap$cut_by %||% NA,
    cluster_rows = config$heatmap$selection$cluster_rows %||% config$heatmap$cluster_rows %||% TRUE,
    cluster_columns = config$heatmap$selection$cluster_columns %||% config$heatmap$cluster_columns %||% c(TRUE, FALSE),
    stringsAsFactors = FALSE
  )
  log_info(
    "Generating selection heatmaps",
    context = list(
      combinations = nrow(param_grid),
      friendly_name = friendly_name,
      outdir = outdir
    )
  )

  file_list <- config$heatmap$selection$file_list
  if (is.null(file_list)) {
    log_skip("Selection heatmaps disabled", context = list(reason = "no file list"))
    return(NULL)
  }

  config$heatmap$selection$file_list %>%
    purrr::list_flatten() %>%
    purrr::imap(~ {
      the_file <- .x
      file_name <- .y

      log_info(
        "Preparing selection heatmap",
        context = list(file = the_file, label = file_name)
      )

      # Read and subset GCT
      file_data <- io_tools$read_genelist_file(the_file)
      subgct_z <- gct_z %>% subset_gct(rid = file_data$id %||% file_data[[1]])

      log_debug(
        "Selection subset dimensions",
        context = list(file = the_file, rows = nrow(subgct_z@mat), cols = ncol(subgct_z@mat))
      )

      quantiles <- quantile(
        probs = seq(0, 1, 0.025),
        na.rm = T
      )
      heatmap_colorbar_bounds <- c(quantiles[["2.5%"]], quantiles[["97.5%"]])

      if (nrow(subgct_z@mat) == 0) {
        log_skip("Selection heatmap skipped - no rows", context = list(file = the_file))
        return(NULL)
      }
      # if (nrow(subgct_z@mat) > 499) {
      #   message(paste("Too many rows in ", the_file, ", skipping"))
      #   return(NULL)
      # }

      # Generate heatmaps for this subset
      param_grid %>%
        purrr::pmap(~ generate_single_selection_heatmap(
          subgct_z = subgct_z,
          cut_by = ..1,
          cluster_rows = ..2,
          cluster_columns = ..3,
          outdir = outdir,
          outname = outname,
          friendly_name = friendly_name,
          file_name = file_name,
          config = config,
          heatmap_colorbar_bounds = heatmap_colorbar_bounds
        ))
    })
}


generate_single_selection_heatmap <- function(subgct_z, cut_by, cluster_rows, cluster_columns, outdir, outname, file_name = "", friendly_name = "", config = NULL,
                                              heatmap_colorbar_bounds = NULL) {
  heatmap_tools <- get_tool_env("heatmap")

  outf <- heatmap_tools$prepare_heatmap_output_path(outdir,
    file.path("heatmap", "selection"),
    outname,
    prefix = file_name,
    dimensions = list(nrow = nrow(subgct_z@mat), ncol(subgct_z@mat)),
    cut_by = cut_by,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns
  )

  log_debug(
    "Prepared selection heatmap",
    context = list(
      outf = outf,
      file = file_name,
      cut_by = ifelse(is.na(cut_by), "<NA>", as.character(cut_by)),
      cluster_rows = cluster_rows,
      cluster_columns = cluster_columns,
      rows = nrow(subgct_z@mat),
      cols = ncol(subgct_z@mat)
    )
  )

  if (file.exists(outf) && config$advanced$replace == FALSE) {
    log_skip(
      "Selection heatmap already exists",
      context = list(outf = outf, replace = config$advanced$replace %||% FALSE)
    )
    return(NULL)
  }


  save_func <- purrr::partial(.f = heatmap_tools$draw_and_save_heatmap, outf = outf)

  # Generate heatmap code
  ht_draw_code <- heatmap_tools$make_heatmap_fromgct(
    gct = subgct_z,
    column_title = file_name,
    cut_by = cut_by,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    meta_to_include = config$heatmap$legend_include %||% NULL,
    meta_to_exclude = config$heatmap$legend_exclude %||% NULL,
    show_row_names = TRUE,
    include_row_annotations = TRUE,
    # column_title = friendly_name,
    optimal_size = TRUE,
    heatmap_colorbar_bounds = heatmap_colorbar_bounds,
    save_func = save_func
  )
  # Draw and save the heatmap
  log_info("Rendering selection heatmap", context = list(outf = outf, file = file_name))
}


run <- function(
    data_dir = NULL,
    output_dir = NULL,
    metadata = NULL,
    config_file = NULL,
    gct_file = NULL,
    save_env = FALSE,
    here_dir = NULL,
    ...) {
  log_info(
    "Starting site-annotate run",
    context = list(
      data_dir = data_dir %||% "<NA>",
      output_dir = output_dir %||% "<NA>",
      gct_file = gct_file %||% "<NA>",
      config_file = config_file %||% "<NA>"
    )
  )
  # knitr::opts_chunk$set(echo = TRUE)
  # suppressPackageStartupMessages(library(rlang))
  # suppressPackageStartupMessages(library(tidyverse))
  # suppressPackageStartupMessages(library(magrittr))
  # suppressPackageStartupMessages(library(purrr))
  # suppressPackageStartupMessages(library(ComplexHeatmap))
  # suppressPackageStartupMessages(library(circlize))
  # suppressPackageStartupMessages(library(reactable))
  # suppressPackageStartupMessages(library(RColorBrewer))
  # suppressPackageStartupMessages(library(here))
  # print(paste0("here is defined as: ", here()))

  library(tictoc)

  tic("setting up env")
  setup_environment() # populates the global environment with the tools and libraries
  toc(log = FALSE)


  if (exists("get_tool_env", where = .GlobalEnv, mode = "function")) {
    get_tool_env <- .GlobalEnv$get_tool_env
  } # TRUE

  io_tools <- get_tool_env("io")
  heatmap_tools <- get_tool_env("heatmap")
  # util_tools <- get_tool_env("utils")
  # limma_tools <- get_tool_env("limma")

  OUTDIR <- output_dir %||% data_dir %||% "."

  #
  config_file <- config_file %||% file.path(here(), "config", "base.toml")
  config <- load_and_validate_config(config_file, data_dir)
  # config <- io_tools$load_config(config_file, root_dir = data_dir)
  log_debug("Configuration sections loaded", context = list(keys = paste(names(config), collapse = ",")))

  z_score_flag <- config$heatmap$z_score %||% TRUE
  z_score_by <- config$heatmap$z_score_by %||% config$heatmap$zscore_by %||% NULL


  if (is.null(gct_file)) {
    gct_file <- config$gct_file
  }
  gct <- process_gct_file(gct_file, config)
  gct <- exclude_samples(gct, config$extra$sample_exclude)
  log_info(
    "GCT ready",
    context = list(rows = nrow(gct@mat), cols = ncol(gct@mat))
  )

  # QC plots (optional; defaults on via [params.qc])
  if (config$qc$do %||% TRUE) {
    qc_tools <- get_tool_env("qc")
    qc_tools$generate_qc_reports(
      gct = gct,
      config = config,
      output_dir = OUTDIR
    )
  } else {
    log_skip("QC disabled", context = list(reason = "config$qc$do == FALSE"))
  }

  .outdir <- file.path(OUTDIR, "export")
  if (!dir.exists(.outdir)) dir.create(.outdir, recursive = TRUE)
  export_prefix <- file.path(.outdir, "export")
  overwrite_export <- config$advanced$overwrite_export %||% FALSE
  # Detect any existing export GCT/GCTX written previously with size suffixes
  existing_exports <- list.files(.outdir, pattern = "^export_.*\\.(gct|gctx)$", full.names = TRUE)
  if (length(existing_exports) > 0 && !overwrite_export) {
    log_skip(
      "Export GCT already exists",
      context = list(n_files = length(existing_exports), overwrite = overwrite_export)
    )
  } else {
    log_info(
      "Writing export GCT",
      context = list(outdir = .outdir, overwrite = overwrite_export)
    )
    gct %>% write_gct(export_prefix)
  }

  # ==========================================
  # NAMES SETUP
  suboutdir <- basename(gct_file) %>%
    basename() %>%
    tools::file_path_sans_ext()

  OUTDIR <- file.path(OUTDIR, suboutdir)
  if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)
  OUTNAME <- "" # don't set an actual name for the internal files saved within the main dir

  FRIENDLY_NAME <- suboutdir %>%
    str_replace_all("_", " ") %>%
    str_replace(., "\\d+$", "") %>%
    str_wrap(width = 60)
  # ==========================================

  util_tools <- get_tool_env("utils")
  gct_z <- util_tools$scale_gct(gct,
    group_by = z_score_by,
    apply_z = z_score_flag
  ) # %cached_by% rlang::hash(c(gct@mat, config$heatmap$z_score_by))



  if (config$heatmap$selection$do %||% FALSE) {
    # Heatmap generation for selections
    generate_selection_heatmaps(gct_z, config, outdir = OUTDIR, friendly_name = FRIENDLY_NAME, outname = OUTNAME)
  }


  if (config$heatmap$do %||% FALSE) {
    generate_global_heatmaps(
      gct_z = gct_z,
      config = config,
      outdir = OUTDIR,
      friendly_name = FRIENDLY_NAME,
      outname = NULL # OUTNAME
    )
  }

  if (config$limma$do) {
    limma_tools <- get_tool_env("limma")
    topTables <- limma_tools$run_limma(gct, config)

    # Write LIMMA summary
    limma_summary <- attr(topTables, "limma_summary")
    if (!is.null(limma_summary)) {
      sumf <- file.path(OUTDIR, "limma", "limma_summary.txt")
      dir.create(dirname(sumf), recursive = TRUE, showWarnings = FALSE)
      writeLines(limma_summary, sumf)
      log_info("Wrote LIMMA summary", context = list(file = sumf))
    }

    topTables %>% purrr::imap(~ { # first loop save tables
      .table <- .x %>% arrange(desc(t))
      contrast_label <- .y
      safe_label <- str_replace_all(contrast_label, "-", "minus")

      .outdir <- file.path(OUTDIR, "limma", "tables")
      .outname <- file.path(.outdir, make.names(paste0(OUTNAME, safe_label, ".tsv")))

      if (!dir.exists(.outdir)) dir.create(.outdir, recursive = TRUE)
      if (!file.exists(.outname) || config$advanced$replace == TRUE) {
        .table %>% write_tsv(.outname)
        log_info("Wrote LIMMA table", context = list(file = .outname))
      }

      .outname <- file.path(.outdir, make.names(paste0(OUTNAME, safe_label, "_pdist.png")))
      if (!file.exists(.outname) || config$advanced$replace == TRUE) {
        png(.outname, width = 7, height = 4.5, units = "in", res = 300)
        limma_tools$plot_p_dist(.table, main_title = contrast_label)
        dev.off()
        #

      log_info("Wrote LIMMA plot", context = list(file = .outname))
      }
    })

    # Volcano plots per contrast
    volcdir <- file.path(OUTDIR, "limma", "volcano")
    if (!dir.exists(volcdir)) dir.create(volcdir, recursive = TRUE)
    p_cut <- config$limma$p_cutoff %||% 0.05
    logfc_cut <- config$limma$logfc_cutoff %||% 1.0
    # Allow either a single value or a list of values
    label_top_n_list <- config$limma$label_top_n_list %||% config$limma$label_top_n %||% 15
    label_top_n_list <- unique(as.integer(unlist(label_top_n_list)))
    label_top_n_list <- label_top_n_list[order(label_top_n_list)]
    topTables %>% purrr::imap(~ {
      df <- .x
      lbl <- .y
      if (!all(c("logFC", "adj.P.Val") %in% colnames(df))) return(NULL)
      df$neglog10 <- -log10(pmax(df$adj.P.Val, .Machine$double.xmin))
      fc_filter_on <- is.finite(logfc_cut) && is.numeric(logfc_cut) && (logfc_cut > 0)
      if (fc_filter_on) {
        df$signif <- (df$adj.P.Val <= p_cut) & (abs(df$logFC) >= logfc_cut)
      } else {
        df$signif <- (df$adj.P.Val <= p_cut)
      }
      df$status <- dplyr::case_when(
        df$signif & (df$logFC >= (if (fc_filter_on) logfc_cut else 0)) ~ "up",
        df$signif & (df$logFC <= (if (fc_filter_on) -logfc_cut else 0)) ~ "down",
        TRUE ~ "ns"
      )
      # Choose labels for significant points and rank once
      has_sitename <- "sitename" %in% colnames(df)
      df$label_txt <- if (has_sitename) df$sitename else df$id
      lab_pool <- df[df$signif, , drop = FALSE]
      if (nrow(lab_pool) > 0) {
        ord <- order(-lab_pool$neglog10, -abs(lab_pool$logFC))
        lab_pool <- lab_pool[ord, , drop = FALSE]
      }
      safe_lbl <- str_replace_all(lbl, "-", "minus")
      # Base plot without labels; label count variants below
      total_n <- nrow(df)
      sig_n <- sum(df$signif, na.rm = TRUE)
      cutoff_text <- sprintf("%d / %d at adj.P.Val < %s", sig_n, total_n, format(p_cut, digits = 3))
      if (fc_filter_on) {
        cutoff_text <- paste0(cutoff_text, sprintf(" and |log2FC| >= %s", format(logfc_cut, digits = 3)))
      }
      p_base <- ggplot2::ggplot(df, ggplot2::aes(x = logFC, y = neglog10)) +
        ggplot2::geom_point(ggplot2::aes(color = status), size = 0.8, alpha = 0.75, show.legend = FALSE) +
        ggplot2::geom_hline(yintercept = -log10(p_cut), linetype = "dashed", color = "#999999") +
        ggplot2::geom_vline(xintercept = c(-logfc_cut, logfc_cut), linetype = "dashed", color = "#999999") +
        ggplot2::scale_color_manual(values = c(ns = "#bdbdbd", up = "#d7191c", down = "#2c7bb6")) +
        ggplot2::theme_bw(base_size = 10) +
        ggplot2::theme(plot.caption.position = "plot", plot.caption = ggplot2::element_text(hjust = 1, size = 8)) +
        ggplot2::labs(title = paste0("Volcano: ", lbl), x = "log2 fold-change", y = expression(-log[10](adj.P.Val)), caption = cutoff_text)

      # Render one plot per requested label count; if fewer sig than requested, actual label n is used in filename
      for (n_req in label_top_n_list) {
        n_actual <- min(n_req, nrow(lab_pool))
        outf <- file.path(volcdir, paste0("volcano_", make.names(safe_lbl), "_n", n_actual, ".pdf"))
        if (file.exists(outf) && config$advanced$replace == FALSE) {
          log_skip("Volcano exists", context = list(file = outf))
          next
        }
        p <- p_base
        if (n_actual > 0) {
          lab_df <- utils::head(lab_pool, n_actual)
          # dynamic label sizing: shrink as label count grows, bounded by label_size_min
          label_size_base <- as.numeric(config$limma$label_size_base %||% 2.6)
          label_size_min <- as.numeric(config$limma$label_size_min %||% 1.6)
          label_size <- label_size_base
          if (n_actual > 30) {
            label_size <- max(label_size_min, label_size_base * sqrt(30 / n_actual))
          }

          if (requireNamespace("ggrepel", quietly = TRUE)) {
            # ggrepel tuning knobs from config
            seg_size <- as.numeric(config$limma$label_segment_size %||% 0.25)
            seg_alpha <- as.numeric(config$limma$label_segment_alpha %||% 0.6)
            box_pad <- as.numeric(config$limma$label_box_padding %||% 0.25)
            point_pad <- as.numeric(config$limma$label_point_padding %||% 0.25)
            repel_force <- as.numeric(config$limma$label_force %||% 1.0)
            max_time <- as.numeric(config$limma$label_max_time %||% 2.0)
            max_overlaps <- config$limma$label_max_overlaps %||% Inf
            # sanitize label direction
            label_dir_raw <- config$limma$label_direction %||% "both"
            valid_dirs <- c("both", "x", "y")
            label_dir <- if (is.character(label_dir_raw) && length(label_dir_raw) >= 1 && label_dir_raw[1] %in% valid_dirs) label_dir_raw[1] else "both"
            p <- p + ggrepel::geom_text_repel(
              data = lab_df,
              ggplot2::aes(label = label_txt),
              size = label_size,
              max.overlaps = max_overlaps,
              min.segment.length = 0.05,
              box.padding = box_pad,
              point.padding = point_pad,
              segment.size = seg_size,
              segment.alpha = seg_alpha,
              force = repel_force,
              max.time = max_time,
              direction = label_dir
            )
          } else {
            p <- p + ggplot2::geom_text(
              data = lab_df,
              ggplot2::aes(label = label_txt),
              size = label_size, vjust = -0.2
            )
          }
        }
        ggplot2::ggsave(outf, p, width = 6.8, height = 5.4)
        log_info("Wrote volcano plot", context = list(file = outf, n_labels = n_actual))
      }
      NULL
    })

    limma_topn <- config$limma$topn %||% c(50, 100)
    param_grid <- expand.grid(
      cut_by = unique(c(config$heatmap$cut_by, NA)),
      cluster_rows = unique(c(TRUE, FALSE) %||% config$limma$cluster_rows %||% config$heatmap$cluster_rows), # unique(c(config$heatmap$cluster_rows, FALSE)) %||% c(TRUE, FALSE),
      cluster_columns = unique(config$limma$cluster_columns %||%
        config$heatmap$cluster_columns %||%
        c(FALSE, TRUE)),
      top_n = config$limma$top_n %||% c(50, 100)
      # file_data = file_data
    )

    gct_z <- util_tools$scale_gct(
      gct,
      group_by = z_score_by,
      apply_z = z_score_flag
    )

    topTables %>% purrr::imap(~ { # second loop plot topdiff
      .table <- .x %>% arrange(desc(abs(t)))
      contrast_label <- .y
      contrast_expr <- unique(.table$contrast_expression)
      contrast_expr <- contrast_expr[!is.na(contrast_expr)]
      if (length(contrast_expr) == 0) contrast_expr <- contrast_label
      contrast_expr <- contrast_expr[1]
      safe_label_file <- str_replace_all(contrast_label, "-", "minus")
      safe_label_dir <- sanitize_path_component(contrast_label)

      .outdir <- file.path(OUTDIR, "limma", "tables")

      column_title <- paste0(OUTNAME, "\ncontrast: ", contrast_label)

      param_grid %>% purrr::pmap(~ {
        cut_by_val <- ..1
        cluster_rows_val <- ..2
        cluster_columns_val <- ..3
        top_n_val <- ..4


        .top_n_val <- top_n_val %/% 2 # Integer divide for balanced selection
        .top_tvals <- .table %>%
          dplyr::distinct(t) %>%
          arrange(desc(t)) %>%
          head(.top_n_val) %>%
          bind_rows(
            .table %>%
              dplyr::distinct(t) %>%
              arrange(t) %>%
              head(.top_n_val)
          )

        .subsel <- .table %>%
          dplyr::filter(t %in% .top_tvals$t) %>%
          arrange(-logFC)

        .gct_z <- gct_z %>% subset_gct(rid = .subsel$id)
        ht_draw_code <- .gct_z %>% heatmap_tools$make_heatmap_fromgct(
          column_title = column_title,
          cut_by = cut_by_val,
          cluster_rows = cluster_rows_val,
          cluster_columns = cluster_columns_val,
          meta_to_include = config$heatmap$legend_include %||% NULL,
          meta_to_exclude = config$heatmap$legend_exclude %||% NULL,
          show_row_names = TRUE,
          optimal_size = TRUE
        )

        .nrow <- nrow(.gct_z@mat)
        .ncol <- ncol(.gct_z@mat)
        height <- 6 + (.nrow * .20)
        width <- 12 + (.ncol * .26)

          .outdir <- file.path(OUTDIR, "limma", "all", safe_label_dir)
        if (!dir.exists(.outdir)) dir.create(.outdir, recursive = TRUE)
        .outf <- file.path(
          .outdir,
          make.names(
            paste0(
                OUTNAME,
                safe_label_file, "_",
                .nrow, "x", .ncol,
                "_cut_", cut_by_val,
                "cr_", cluster_rows_val,
                "cc_", cluster_columns_val,
                ".pdf"
            )
          )
        )

        log_info("Rendering LIMMA heatmap", context = list(outf = .outf))
        if (file.exists(.outf) && config$advanced$replace == FALSE) {
          log_skip(
            "LIMMA heatmap already exists",
            context = list(outf = .outf, replace = config$advanced$replace %||% FALSE)
          )
          return(NULL)
        }
        tryCatch(
          {
            cairo_pdf(.outf, width = width, height = height)
            ht_draw_code()
            dev.off()
          },
          error = function(e) print(e)
        )

        # ===============================

        # Extract the two values from the contrast string
        # this all needs to be in a trycatch
        tryCatch(
          {
            column_names <- names(gct@cdesc)

            values <- stringr::str_extract_all(contrast_expr, "(?<=- )[A-Za-z0-9_]+")[[1]]
            matched_column <- column_names[purrr::map_lgl(column_names, ~ str_detect(contrast_expr, .x))]
            if (length(matched_column) == 0) {
              warning(paste0("No column found matching contrast: ", contrast_expr))
              return()
            }
            col_name <- matched_column[1] # Assuming there's one match; adjust if multiple matches needed
            cleaned_contrast <- str_replace_all(contrast_expr, col_name, "")
            group_values <- str_split(cleaned_contrast, "\\s?-\\s?") %>%
              unlist() %>%
              str_trim()

            .cdesc <- gct@cdesc %>% dplyr::filter(!!sym(col_name) %in% group_values)
            sub_gct_z <- gct %>%
              subset_gct(cid = rownames(.cdesc)) %>%
              util_tools$scale_gct() # %>%
            .sub_gct_z <- sub_gct_z %>% subset_gct(rid = .subsel$id)

            # Skip subset heatmap if the selected columns match the full set (duplicate of 'all')
            if (identical(sort(colnames(.sub_gct_z@mat)), sort(colnames(gct@mat)))) {
              log_skip("Subset equals full set; skipping duplicate", context = list(contrast = contrast_label))
              return(NULL)
            }


            .nrow <- nrow(.top_tvals)
            .ncol <- ncol(.sub_gct_z@mat)

            .outdir <- file.path(OUTDIR, "limma", "subsets", safe_label_dir)
        if (!dir.exists(.outdir)) dir.create(.outdir, recursive = TRUE)
        .outf <- file.path(
          .outdir,
          make.names(
            paste0(
                  OUTNAME,
                  safe_label_file, "_",
                  .nrow, "x", .ncol,
                  "_cut_", cut_by_val,
                  "cr_", cluster_rows_val,
                  "cc_", cluster_columns_val,
                  ".pdf"
                )
              )
            )

            ht_draw_code <- .sub_gct_z %>% heatmap_tools$make_heatmap_fromgct(
              column_title = column_title,
              cut_by = cut_by_val,
              cluster_rows = cluster_rows_val,
              cluster_columns = cluster_columns_val,
              meta_to_include = config$heatmap$legend_include %||% NULL,
              meta_to_exclude = config$heatmap$legend_exclude %||% NULL,
              show_row_names = TRUE,
              optimal_size = TRUE
            )

            .nrow <- nrow(.sub_gct_z@mat)
            .ncol <- ncol(.sub_gct_z@mat)
            height <- 6 + (.nrow * .20)
            width <- 12 + (.ncol * .26)

            log_info("Rendering LIMMA subset heatmap", context = list(outf = .outf))
            if (file.exists(.outf) && config$advanced$replace == FALSE) {
              log_skip(
                "LIMMA subset heatmap already exists",
                context = list(outf = .outf, replace = config$advanced$replace %||% FALSE)
              )
              return(NULL)
            }
            tryCatch(
              {
                cairo_pdf(.outf, width = width, height = height)
                ht_draw_code()
                dev.off()
              },
              error = function(e) print(e)
            )
          },
          error = function(e) print(e)
        )
      }) # end of pmap of param grid
    }) # end of imap of topTables
  } # end of if limma


  if (config$pca$do) {
    pca_tools <- get_tool_env("pca")

    param_grid <- expand.grid(
      scale = c(TRUE),
      colby = config$pca$colby %||% config$pca$color %||% NA,
      top_loadings = config$pca$top_loadings %||% 5,
      markby = config$pca$markby %||% config$pca$marker %||% NA,
      #

      stringsAsFactors = FALSE
    )

    param_grid %>%
      purrr::pmap(~ {
        do_scale <- ..1
        colby <- ..2
        ntopLoadings <- ..3
        markby <- ..4
        rids_to_keep <- gct@rdesc %>%
          distinct(fifteenmer, .keep_all = TRUE) %>%
          pull(id)
        .gct <- gct %>% subset_gct(rid = rids_to_keep)

        pca_obj <- pca_tools$do_pca(.gct, scale = do_scale)

        .outdir <- file.path(OUTDIR, "pca")
        if (!dir.exists(.outdir)) dir.create(.outdir, recursive = TRUE)

        .outf <- file.path(
          .outdir,
          make.names(
            paste0(
              "pca",
              "_scale_", do_scale,
              "_col_", colby,
              ifelse(!(is.null(markby) && is.na(markby)), paste0("_markby_", markby), ""),
              ifelse(!is.null(ntopLoadings) && ntopLoadings > 0, paste0("_nload_", ntopLoadings), ""),
              ".pdf"
            )
          )
        )

        plts <- pca_tools$plot_biplot(pca_obj,
          top_pc = config$pca$max_pc %||% 3,
          showLoadings = T,
          labSize = config$pca$label_size %||% 1.8,
          pointSize = config$pca$point_size %||% 3,
          sizeLoadingsNames = 2,
          colby = colby,
          shape = markby,
          encircle = !is.null(colby),
          title = FRIENDLY_NAME,
          ntopLoadings = ntopLoadings,
          fig_width = config$pca$width %||% 8.4,
          fig_height = config$pca$height %||% 7.6,
          basename = .outf
        )


        # not working yet
        # pca_tools$make_heatmap_from_loadings(gct_obj=.gct, pca_object=pca_obj)
      }) # end of map
    # browser()

    if (config$totalprotein$do) {
      # this is not functional yet
      # sitediff_dir <- file.path(OUTDIR, "limma", "tables")
      # total_protein_dir <- file.path(OUTDIR, "limma", "tables", "totalprotein")
      # if (!dir.exists(total_protein)) dir.create(total_protein, recursive = TRUE)
      # comparisons <- io_tools$match_comparisons(sitediff_dir, total_protein_dir)
    }
  }

  log_info("site-annotate run complete", context = list(outdir = OUTDIR))
}
# end of run
