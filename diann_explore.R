#!/usr/bin/env Rscript

# Minimal, dependency-light R exploration of DIA-NN outputs
# - Uses tidyverse/ggplot2 style for summaries
# - Reads DIA-NN TSV and (optionally) Parquet (via arrow if installed)
# - Writes figures to results/explore/{phos,prof}

suppressPackageStartupMessages({
  ok_tidyverse <- requireNamespace("tidyverse", quietly = TRUE)
  if (ok_tidyverse) {
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(readr)
    library(stringr)
    library(purrr)
    library(forcats)
  } else {
    stop("Please install tidyverse: install.packages('tidyverse')")
  }
})

have_arrow <- requireNamespace("arrow", quietly = TRUE)
have_complexheatmap <- requireNamespace("ComplexHeatmap", quietly = TRUE)

# Logging -------------------------------------------------------------------

.now <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
.log <- function(level, fmt, ...) {
  msg <- sprintf(fmt, ...)
  cat(sprintf("[%s] [%s] %s\n", .now(), level, msg))
}
log_info  <- function(fmt, ...) .log("INFO", fmt, ...)
log_warn  <- function(fmt, ...) .log("WARN", fmt, ...)
log_error <- function(fmt, ...) .log("ERROR", fmt, ...)

# I/O helpers ---------------------------------------------------------------

dir_exists <- function(path) {
  isTRUE(file.info(path)$isdir)
}

ensure_dir <- function(path) {
  if (!dir_exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

list_files <- function(root, pattern) {
  list.files(root, pattern = pattern, recursive = TRUE, full.names = TRUE)
}

# Paths (override with env vars if desired) --------------------------------

# PHOS_DIR   <- Sys.getenv("DIANN_PHOS_DIR", unset = "data/processed/phos")
PROF_DIR   <- Sys.getenv("DIANN_PROF_DIR", unset = "data/processed/prof")
PHOS_DIR   <- Sys.getenv("DIANN_PHOS_DIR", unset = "/mnt/amms06/06262024/MSPC001482/data/processed/")
#PROF_DIR   <- Sys.getenv("DIANN_PHOS_DIR", unset = "/mnt/amms06/06262024/MSPC001482/data/processed/prof")
OUT_DIR    <- Sys.getenv("DIANN_OUT_DIR",  unset = "results/explorephos")

# Derive PROC_ROOT if not explicitly set
PROC_ROOT  <- Sys.getenv("DIANN_PROC_ROOT", unset = NA_character_)
if (is.na(PROC_ROOT) || !nzchar(PROC_ROOT)) {
  if (nzchar(PHOS_DIR) && grepl("/phos/?$", PHOS_DIR)) {
    PROC_ROOT <- dirname(PHOS_DIR)
  } else if (nzchar(PROF_DIR) && grepl("/prof/?$", PROF_DIR)) {
    PROC_ROOT <- dirname(PROF_DIR)
  } else if (nzchar(PHOS_DIR)) {
    PROC_ROOT <- PHOS_DIR
  } else if (nzchar(PROF_DIR)) {
    PROC_ROOT <- PROF_DIR
  } else {
    PROC_ROOT <- "data/processed"
  }
}
RAW_SUMMARY_PATH <- Sys.getenv("DIANN_RAW_SUMMARY", unset = "results/raw/raw_summary.csv")

ensure_dir(file.path(OUT_DIR, "phos"))
ensure_dir(file.path(OUT_DIR, "prof"))

if (have_arrow) {
  # Print brief Arrow build info to help debug codec issues (e.g., zstd not built)
  info <- tryCatch(arrow::arrow_info(), error = function(e) NULL)
  if (!is.null(info)) {
    log_info("arrow version: %s; features: %s", as.character(info$version), paste(names(info$capabilities)[unlist(info$capabilities)], collapse = ","))
  } else {
    log_warn("arrow installed but arrow_info() failed; Parquet will be best-effort")
  }
} else {
  log_info("arrow not installed; Parquet will be skipped")
}

# Inference helpers ---------------------------------------------------------

infer_method <- function(spectrum_file) {
  # Extract run/method indicator from names like: 57013_1_TOF_... -> method 1
  m <- stringr::str_match(spectrum_file, "_([0-9])_TOF")
  out <- ifelse(is.na(m[,2]), NA_character_, paste0("meth", m[,2]))
  out
}

infer_meta_experiment <- function(path_or_file) {
  if (is.null(path_or_file) || !length(path_or_file)) return(rep(NA_character_, length(path_or_file)))
  m <- stringr::str_match(path_or_file, "(meth[0-9]+)")
  out <- ifelse(is.na(m[,2]), NA_character_, m[,2])
  out
}

infer_ms_method <- function(spectrum_file) {
  # From SpectrumFile pattern, e.g., 57013_1_TOF -> meth1
  infer_method(spectrum_file)
}

infer_sample_id <- function(spectrum_file) {
  # Leading digits prior to first underscore, e.g., '57013' or '57014'
  m <- stringr::str_match(spectrum_file, "^([0-9]+)_")
  ifelse(is.na(m[,2]), NA_character_, m[,2])
}

normalize_meta <- function(x) {
  out <- as.character(x)
  out[is.na(out) | !nzchar(out)] <- "meta_unknown"
  out
}

extract_meta_from_path <- function(path) {
  if (is.null(path) || !length(path)) return(rep(NA_character_, length(path)))
  bn <- basename(path)
  match <- grepl("^meth[0-9]+$", bn, ignore.case = TRUE)
  out <- bn
  out[!match] <- NA_character_
  out
}

sanitize_token <- function(x) {
  out <- gsub("[^A-Za-z0-9]+", "-", x)
  out <- gsub("(^-+|-+$)", "", out)
  ifelse(out == "", "NA", out)
}

calc_plot_width <- function(n_meta, base = 8, extra_per = 1.5, max_width = 18) {
  if (is.na(n_meta) || n_meta < 1) n_meta <- 1
  width <- base + (n_meta - 1) * extra_per
  min(width, max_width)
}

order_meta_levels <- function(meta_vec) {
  vals <- normalize_meta(meta_vec)
  uniq_vals <- unique(vals)
  uniq_vals <- uniq_vals[!is.na(uniq_vals)]
  if (!length(uniq_vals)) return(NULL)
  numeric_part <- suppressWarnings(as.numeric(stringr::str_extract(uniq_vals, "(?<=meth)\\d+")))
  is_meth <- !is.na(numeric_part)
  is_unknown <- uniq_vals == "meta_unknown"
  order_df <- tibble::tibble(
    meta = uniq_vals,
    is_meth = is_meth,
    num = numeric_part,
    is_unknown = is_unknown
  ) %>%
    arrange(is_unknown, dplyr::desc(is_meth), num, meta)
  order_df$meta
}

format_count <- function(x) {
  ifelse(is.na(x), "0", format(round(as.numeric(x)), big.mark = ",", trim = TRUE))
}

filter_phospho_psms <- function(df) {
  if (!nrow(df) || !"SequenceModi" %in% names(df)) return(tibble::tibble())
  pat <- "(S|T|Y)\\(pho\\)|(S|T|Y)\\(UniMod:21\\)"
  df %>%
    dplyr::filter(!is.na(SequenceModi)) %>%
    dplyr::filter(stringr::str_detect(SequenceModi, pat))
}

dedupe_psms <- function(df) {
  if (!nrow(df) || !"SequenceModi" %in% names(df) || !"Charge" %in% names(df)) return(df)
  time_cols <- intersect(c("RTmin", "RT", "Retention.Time"), names(df))
  time_col <- if (length(time_cols)) time_cols[1] else NULL
  key_cols <- c("SpectrumFile", "SequenceModi", "Charge")
  if (!is.null(time_col)) key_cols <- c(key_cols, time_col)
  df %>% dplyr::distinct(dplyr::across(any_of(key_cols)), .keep_all = TRUE)
}

expand_phospho_residues <- function(df) {
  if (!nrow(df) || !"SequenceModi" %in% names(df)) return(tibble::tibble())
  residues <- stringr::str_extract_all(df$SequenceModi, "[STY](?=\\(pho\\)|\\(UniMod:21\\))")
  res_tbl <- tibble::tibble(row_id = seq_len(nrow(df)), Residue = residues)
  res_tbl <- tidyr::unnest(res_tbl, Residue)
  if (!nrow(res_tbl)) return(tibble::tibble())
  expanded <- df[res_tbl$row_id, , drop = FALSE]
  expanded$Residue <- res_tbl$Residue
  expanded
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

log10p1 <- function(x) log10(pmax(as.numeric(x), 0) + 1)

# Reader functions ----------------------------------------------------------

read_diann_reports_parquet <- function(root) {
  files <- list_files(root, pattern = "report.*\\.parquet$")
  if (!length(files)) return(NULL)
  if (!have_arrow) {
    log_info("Skipping Parquet in %s (arrow not installed)", root)
    return(NULL)
  }
  log_info("Found %d Parquet report file(s) under %s", length(files), root)
  res <- tryCatch({
    ds <- arrow::open_dataset(files)
    df <- as.data.frame(ds)
    tibble::as_tibble(df) %>% mutate(source_file = NA_character_)
  }, error = function(e) {
    msg <- conditionMessage(e)
    if (grepl("zstd|codec", msg, ignore.case = TRUE)) {
      log_warn("Parquet read failed: %s", msg)
      log_warn("Your Arrow build may lack ZSTD support; skipping Parquet. Consider reinstalling arrow with zstd.")
    } else {
      log_warn("Parquet read failed: %s", msg)
    }
    NULL
  })
  res
}

read_diann_psms_quant <- function(meta_roots) {
  if (is.null(meta_roots) || !nrow(meta_roots)) return(tibble::tibble())

  files_tbl <- purrr::map2_dfr(meta_roots$path, meta_roots$meta_experiment, function(path, meta) {
    files <- list_files(path, pattern = "psms_.*QUANT\\.tsv$")
    if (!length(files)) {
      log_warn("  No QUANT files under %s (meta_experiment=%s)", path, meta)
      return(tibble::tibble())
    }
    tibble::tibble(
      file = files,
      meta_experiment_root = meta,
      meta_root_path = path,
      file_parent = dirname(files)
    )
  })

  if (!nrow(files_tbl)) {
    log_warn("No QUANT files found under meta roots: %s", paste(unique(meta_roots$path), collapse = ", "))
    return(tibble::tibble())
  }

  log_info("Reading %d QUANT file(s) across %d meta roots", nrow(files_tbl), nrow(meta_roots))
  max_show <- 20
  log_file_listing <- function(tbl_slice) {
    purrr::pwalk(tbl_slice, function(file, meta_experiment_root, meta_root_path, file_parent) {
      label <- meta_experiment_root
      if (is.null(label) || !nzchar(label)) label <- "NA"
      log_info("  QUANT[%s] %s", label, file)
    })
  }
  if (nrow(files_tbl) <= max_show) {
    log_file_listing(files_tbl)
  } else {
    log_file_listing(files_tbl[seq_len(max_show), , drop = FALSE])
    log_info("  ... %d more QUANT files omitted from log ...", nrow(files_tbl) - max_show)
  }
  safe_read <- purrr::safely(function(f) {
    suppressMessages(readr::read_tsv(f, show_col_types = FALSE, guess_max = 100000))
  })

  out <- purrr::map(files_tbl$file, safe_read)
  errs <- purrr::keep(out, ~ !is.null(.x$error))
  if (length(errs)) {
    purrr::iwalk(errs, function(x, i) log_warn("Failed to read QUANT: %s -> %s", files_tbl$file[[i]], conditionMessage(x$error)))
  }

  ok_idx <- which(vapply(out, function(x) !is.null(x$result), logical(1)))
  if (!length(ok_idx)) return(tibble::tibble())

  dfs <- purrr::map(ok_idx, function(i) out[[i]]$result)
  files_ok <- files_tbl[ok_idx, , drop = FALSE]
  nrows <- vapply(dfs, nrow, integer(1))
  if (!sum(nrows)) return(tibble::tibble())

  outdf <- dplyr::bind_rows(dfs, .id = NULL)
  outdf$source_file <- rep(files_ok$file, times = nrows)
  outdf$meta_root <- rep(files_ok$meta_experiment_root, times = nrows)
  outdf$meta_root_path <- rep(files_ok$meta_root_path, times = nrows)
  outdf$file_parent <- rep(files_ok$file_parent, times = nrows)
  outdf$meta_from_parent <- extract_meta_from_path(outdf$file_parent)

  outdf %>%
    mutate(
      SpectrumFile = dplyr::coalesce(.data$SpectrumFile, .data$`SpectrumFile`),
      ms_method = infer_ms_method(SpectrumFile),
      meta_experiment = normalize_meta(dplyr::coalesce(meta_from_parent, meta_root, infer_meta_experiment(source_file), ms_method)),
      sample_id = infer_sample_id(SpectrumFile)
    )
}

read_diann_psms_qual <- function(meta_roots) {
  if (is.null(meta_roots) || !nrow(meta_roots)) return(tibble::tibble())

  files_tbl <- purrr::map2_dfr(meta_roots$path, meta_roots$meta_experiment, function(path, meta) {
    files <- list_files(path, pattern = "psms_.*QUAL\\.tsv$")
    if (!length(files)) {
      log_warn("  No QUAL files under %s (meta_experiment=%s)", path, meta)
      return(tibble::tibble())
    }
    tibble::tibble(
      file = files,
      meta_experiment_root = meta,
      meta_root_path = path,
      file_parent = dirname(files)
    )
  })

  if (!nrow(files_tbl)) {
    log_warn("No QUAL files found under meta roots: %s", paste(unique(meta_roots$path), collapse = ", "))
    return(tibble::tibble())
  }

  log_info("Reading %d QUAL file(s) across %d meta roots", nrow(files_tbl), nrow(meta_roots))
  max_show <- 20
  log_file_listing <- function(tbl_slice) {
    purrr::pwalk(tbl_slice, function(file, meta_experiment_root, meta_root_path, file_parent) {
      label <- meta_experiment_root
      if (is.null(label) || !nzchar(label)) label <- "NA"
      log_info("  QUAL[%s] %s", label, file)
    })
  }
  if (nrow(files_tbl) <= max_show) {
    log_file_listing(files_tbl)
  } else {
    log_file_listing(files_tbl[seq_len(max_show), , drop = FALSE])
    log_info("  ... %d more QUAL files omitted from log ...", nrow(files_tbl) - max_show)
  }
  safe_read <- purrr::safely(function(f) {
    suppressMessages(readr::read_tsv(f, show_col_types = FALSE, guess_max = 100000))
  })

  out <- purrr::map(files_tbl$file, safe_read)
  errs <- purrr::keep(out, ~ !is.null(.x$error))
  if (length(errs)) {
    purrr::iwalk(errs, function(x, i) log_warn("Failed to read QUAL: %s -> %s", files_tbl$file[[i]], conditionMessage(x$error)))
  }

  ok_idx <- which(vapply(out, function(x) !is.null(x$result), logical(1)))
  if (!length(ok_idx)) return(tibble::tibble())

  dfs <- purrr::map(ok_idx, function(i) out[[i]]$result)
  files_ok <- files_tbl[ok_idx, , drop = FALSE]
  nrows <- vapply(dfs, nrow, integer(1))
  if (!sum(nrows)) return(tibble::tibble())

  outdf <- dplyr::bind_rows(dfs, .id = NULL)
  outdf$source_file <- rep(files_ok$file, times = nrows)
  outdf$meta_root <- rep(files_ok$meta_experiment_root, times = nrows)
  outdf$meta_root_path <- rep(files_ok$meta_root_path, times = nrows)
  outdf$file_parent <- rep(files_ok$file_parent, times = nrows)
  outdf$meta_from_parent <- extract_meta_from_path(outdf$file_parent)

  outdf %>%
    mutate(
      SpectrumFile = dplyr::coalesce(.data$SpectrumFile, .data$`SpectrumFile`),
      ms_method = infer_ms_method(SpectrumFile),
      meta_experiment = normalize_meta(dplyr::coalesce(meta_from_parent, meta_root, infer_meta_experiment(source_file), ms_method)),
      sample_id = infer_sample_id(SpectrumFile)
    )
}

# Plot helpers --------------------------------------------------------------

ggsave_quiet <- function(filename, plot, width = 7, height = 5) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)
  tryCatch({
    suppressMessages(ggplot2::ggsave(filename, plot, width = width, height = height, limitsize = FALSE))
    log_info("Wrote %s", filename)
  }, error = function(e) {
    log_warn("Failed to save %s: %s", filename, conditionMessage(e))
  })
}

theme_min <- function() theme_bw(base_size = 12) + theme(panel.grid = element_blank())

# Counts per run / file -----------------------------------------------------

plot_psm_counts_per_run <- function(qt, title = "PSMs per run (deduplicated)") {
  if (!nrow(qt)) return(NULL)
  df <- qt %>%
    mutate(meta_experiment = normalize_meta(meta_experiment)) %>%
    group_by(meta_experiment, SpectrumFile) %>%
    summarise(n_psms = n(), .groups = "drop") %>%
    arrange(meta_experiment, desc(n_psms)) %>%
    mutate(SpectrumFile = forcats::fct_reorder(SpectrumFile, n_psms))
  meta_levels <- order_meta_levels(df$meta_experiment)
  if (!is.null(meta_levels)) {
    df <- df %>% mutate(meta_experiment = factor(meta_experiment, levels = meta_levels))
  }
  threshold <- max(df$n_psms) * 0.15
  inside_df <- df %>% dplyr::filter(n_psms >= threshold)
  outside_df <- df %>% dplyr::filter(n_psms < threshold)
  threshold <- max(df$n_psms) * 0.15
  inside_df <- df %>% dplyr::filter(n_psms >= threshold)
  outside_df <- df %>% dplyr::filter(n_psms < threshold)
  p <- ggplot(df, aes(x = SpectrumFile, y = n_psms, fill = SpectrumFile)) +
    geom_col(show.legend = FALSE) +
    coord_flip() +
    labs(title = title, x = "Run (SpectrumFile)", y = "PSMs (rows in QUANT)") +
    geom_text(data = inside_df, aes(y = n_psms / 2, label = format_count(n_psms)),
              color = "white", size = 3, fontface = "bold") +
    geom_text(data = outside_df, aes(y = n_psms * 1.05 + 1, label = format_count(n_psms)),
              color = "#222222", size = 3, fontface = "bold") +
    theme_min()
  if (!is.null(meta_levels) && length(meta_levels) > 1) {
    p <- p + facet_grid(meta_experiment ~ ., scales = "free_y", space = "free_y") +
      theme(strip.background = element_blank(), strip.text.y = element_text(face = "bold"))
  }
  p
}

plot_peptide_counts_per_run <- function(qt, title = "Unique peptides per run (deduplicated)") {
  if (!nrow(qt) || !"SequenceModi" %in% names(qt)) return(NULL)
  df <- qt %>%
    mutate(meta_experiment = normalize_meta(meta_experiment)) %>%
    group_by(meta_experiment, SpectrumFile) %>%
    summarise(n_peptides = n_distinct(SequenceModi), .groups = "drop") %>%
    arrange(meta_experiment, desc(n_peptides)) %>%
    mutate(SpectrumFile = forcats::fct_reorder(SpectrumFile, n_peptides))
  meta_levels <- order_meta_levels(df$meta_experiment)
  if (!is.null(meta_levels)) {
    df <- df %>% mutate(meta_experiment = factor(meta_experiment, levels = meta_levels))
  }
  threshold <- max(df$n_peptides) * 0.15
  inside_df <- df %>% dplyr::filter(n_peptides >= threshold)
  outside_df <- df %>% dplyr::filter(n_peptides < threshold)
  p <- ggplot(df, aes(x = SpectrumFile, y = n_peptides, fill = SpectrumFile)) +
    geom_col(show.legend = FALSE) +
    coord_flip() +
    labs(title = title, x = "Run (SpectrumFile)", y = "Unique peptides (SequenceModi)") +
    geom_text(data = inside_df, aes(y = n_peptides / 2, label = format_count(n_peptides)),
              color = "white", size = 3, fontface = "bold") +
    geom_text(data = outside_df, aes(y = n_peptides * 1.05 + 1, label = format_count(n_peptides)),
              color = "#222222", size = 3, fontface = "bold") +
    theme_min()
  if (!is.null(meta_levels) && length(meta_levels) > 1) {
    p <- p + facet_grid(meta_experiment ~ ., scales = "free_y", space = "free_y") +
      theme(strip.background = element_blank(), strip.text.y = element_text(face = "bold"))
  }
  p
}

plot_gene_counts_per_run <- function(qt, title = "Unique genes per run (deduplicated)") {
  gene_col <- dplyr::case_when(
    "GeneID" %in% names(qt) ~ "GeneID",
    "Genes" %in% names(qt) ~ "Genes",
    TRUE ~ NA_character_
  )
  if (is.na(gene_col)) return(NULL)
  df <- qt %>%
    mutate(meta_experiment = normalize_meta(meta_experiment)) %>%
    group_by(meta_experiment, SpectrumFile) %>%
    summarise(n_genes = n_distinct(.data[[gene_col]]), .groups = "drop") %>%
    arrange(meta_experiment, desc(n_genes)) %>%
    mutate(SpectrumFile = forcats::fct_reorder(SpectrumFile, n_genes))
  meta_levels <- order_meta_levels(df$meta_experiment)
  if (!is.null(meta_levels)) {
    df <- df %>% mutate(meta_experiment = factor(meta_experiment, levels = meta_levels))
  }
  threshold <- max(df$n_genes) * 0.15
  inside_df <- df %>% dplyr::filter(n_genes >= threshold)
  outside_df <- df %>% dplyr::filter(n_genes < threshold)
  p <- ggplot(df, aes(x = SpectrumFile, y = n_genes, fill = SpectrumFile)) +
    geom_col(show.legend = FALSE) +
    coord_flip() +
    labs(title = title, x = "Run (SpectrumFile)", y = "Unique genes") +
    geom_text(data = inside_df, aes(y = n_genes / 2, label = format_count(n_genes)),
              color = "white", size = 3, fontface = "bold") +
    geom_text(data = outside_df, aes(y = n_genes * 1.05 + 1, label = format_count(n_genes)),
              color = "#222222", size = 3, fontface = "bold") +
    theme_min()
  if (!is.null(meta_levels) && length(meta_levels) > 1) {
    p <- p + facet_grid(meta_experiment ~ ., scales = "free_y", space = "free_y") +
      theme(strip.background = element_blank(), strip.text.y = element_text(face = "bold"))
  }
  p
}

# Distributions -------------------------------------------------------------

plot_intensity_hist <- function(qt, title = "Intensity distribution", suffix = "") {
  cand <- c("SequenceArea", "PrecursorArea_dstrAdj", "PrecursorArea", "ReporterIntensity")
  col <- cand[cand %in% names(qt)][1]
  if (is.na(col)) return(NULL)
  df <- qt %>% mutate(meta_experiment = normalize_meta(meta_experiment), val = log10p1(.data[[col]]))
  meta_levels <- order_meta_levels(df$meta_experiment)
  if (!is.null(meta_levels)) {
    df <- df %>% mutate(meta_experiment = factor(meta_experiment, levels = meta_levels))
  }
  title_full <- if (nzchar(suffix)) paste0(title, " - ", suffix) else title
  p <- ggplot(df, aes(x = val)) +
    geom_histogram(bins = 60, color = "white") +
    labs(title = paste0(title_full, " (", col, ")"), x = paste0("log10(1+", col, ")"), y = "Count") +
    theme_min()
  if (!is.null(meta_levels) && length(meta_levels) > 1) {
    p <- p + facet_wrap(~meta_experiment)
  }
  p
}

plot_charge_bar <- function(qual, title = "Charge state counts", suffix = "") {
  if (!nrow(qual) || !"Charge" %in% names(qual)) return(NULL)
  df <- qual %>%
    mutate(meta_experiment = normalize_meta(meta_experiment), Charge = factor(Charge)) %>%
    count(meta_experiment, Charge, name = "n")
  meta_levels <- order_meta_levels(df$meta_experiment)
  if (!is.null(meta_levels)) {
    df <- df %>% mutate(meta_experiment = factor(meta_experiment, levels = meta_levels))
  }
  title_full <- if (nzchar(suffix)) paste0(title, " - ", suffix) else title
  p <- ggplot(df, aes(x = Charge, y = n, fill = Charge)) +
    geom_col(show.legend = FALSE) +
    labs(title = title_full, x = "Charge state", y = "Count") +
    theme_min()
  if (!is.null(meta_levels) && length(meta_levels) > 1) {
    p <- p + facet_wrap(~meta_experiment)
  }
  p
}

plot_rt_hist <- function(qual, title = "Retention time distribution", suffix = "") {
  cand <- c("RTmin", "RT", "Retention.Time")
  col <- cand[cand %in% names(qual)][1]
  if (is.na(col)) return(NULL)
  df <- qual %>% mutate(meta_experiment = normalize_meta(meta_experiment))
  meta_levels <- order_meta_levels(df$meta_experiment)
  if (!is.null(meta_levels)) {
    df <- df %>% mutate(meta_experiment = factor(meta_experiment, levels = meta_levels))
  }
  title_full <- if (nzchar(suffix)) paste0(title, " - ", suffix) else title
  p <- ggplot(df, aes(x = .data[[col]])) +
    geom_histogram(bins = 60, color = "white") +
    labs(title = paste0(title_full, " (", col, ")"), x = col, y = "Count") +
    theme_min()
  if (!is.null(meta_levels) && length(meta_levels) > 1) {
    p <- p + facet_wrap(~meta_experiment)
  }
  p
}

plot_im_hist <- function(qual, title = "Ion mobility distribution", suffix = "") {
  cand <- c("IM", "Ion.Mobility")
  col <- cand[cand %in% names(qual)][1]
  if (is.na(col)) return(NULL)
  df <- qual %>% mutate(meta_experiment = normalize_meta(meta_experiment))
  meta_levels <- order_meta_levels(df$meta_experiment)
  if (!is.null(meta_levels)) {
    df <- df %>% mutate(meta_experiment = factor(meta_experiment, levels = meta_levels))
  }
  title_full <- if (nzchar(suffix)) paste0(title, " - ", suffix) else title
  p <- ggplot(df, aes(x = .data[[col]])) +
    geom_histogram(bins = 60, color = "white") +
    labs(title = paste0(title_full, " (", col, ")"), x = col, y = "Count") +
    theme_min()
  if (!is.null(meta_levels) && length(meta_levels) > 1) {
    p <- p + facet_wrap(~meta_experiment)
  }
  p
}

plot_charge_bar_residue <- function(df, title = "Charge state counts by residue") {
  if (!nrow(df) || !"Charge" %in% names(df) || !"Residue" %in% names(df)) return(NULL)
  df <- df %>%
    mutate(meta_experiment = normalize_meta(meta_experiment), Charge = factor(Charge), Residue = factor(Residue, levels = c("S", "T", "Y"))) %>%
    count(meta_experiment, Residue, Charge, name = "n")
  meta_levels <- order_meta_levels(df$meta_experiment)
  if (!is.null(meta_levels)) {
    df <- df %>% mutate(meta_experiment = factor(meta_experiment, levels = meta_levels))
  }
  ggplot(df, aes(x = Charge, y = n, fill = Charge)) +
    geom_col(show.legend = FALSE) +
    labs(title = title, x = "Charge state", y = "Count") +
    facet_grid(meta_experiment ~ Residue, scales = "free_y") +
    theme_min()
}

plot_hist_residue <- function(df, value_col, title, x_lab = NULL, bins = 60) {
  if (!nrow(df) || !value_col %in% names(df) || !"Residue" %in% names(df)) return(NULL)
  df <- df %>% mutate(meta_experiment = normalize_meta(meta_experiment), Residue = factor(Residue, levels = c("S", "T", "Y")))
  meta_levels <- order_meta_levels(df$meta_experiment)
  if (!is.null(meta_levels)) {
    df <- df %>% mutate(meta_experiment = factor(meta_experiment, levels = meta_levels))
  }
  x_label <- if (is.null(x_lab)) value_col else x_lab
  ggplot(df, aes(x = .data[[value_col]])) +
    geom_histogram(bins = bins, color = "white") +
    labs(title = title, x = x_label, y = "Count") +
    facet_grid(meta_experiment ~ Residue, scales = "free_y") +
    theme_min()
}

# Aggregate counts per meta experiment --------------------------------------

plot_counts_per_meta <- function(qt, what = c("psm", "peptide", "gene"), title = NULL) {
  what <- match.arg(what)
  if (!nrow(qt)) return(NULL)
  qt <- qt %>% mutate(meta_experiment = normalize_meta(meta_experiment))
  gene_col <- dplyr::case_when(
    what == "gene" && "GeneID" %in% names(qt) ~ "GeneID",
    what == "gene" && "Genes" %in% names(qt) ~ "Genes",
    TRUE ~ NA_character_
  )
  df <- switch(what,
    psm = qt %>% group_by(meta_experiment) %>% summarise(value = n(), .groups = "drop"),
    peptide = if ("SequenceModi" %in% names(qt)) qt %>% group_by(meta_experiment) %>% summarise(value = n_distinct(SequenceModi), .groups = "drop") else NULL,
    gene = if (!is.na(gene_col)) qt %>% group_by(meta_experiment) %>% summarise(value = n_distinct(.data[[gene_col]]), .groups = "drop") else NULL
  )
  if (is.null(df)) return(NULL)
  meta_levels <- order_meta_levels(df$meta_experiment)
  if (!is.null(meta_levels)) {
    df <- df %>% mutate(meta_experiment = factor(meta_experiment, levels = meta_levels)) %>% arrange(meta_experiment)
  } else {
    df <- df %>% arrange(meta_experiment)
  }
  y_lab <- switch(what,
    psm = "PSMs (rows)",
    peptide = "Unique peptides",
    gene = "Unique genes"
  )
  if (is.null(title)) {
    title <- paste(y_lab, "by meta experiment")
  }
  ggplot(df, aes(x = meta_experiment, y = value, fill = meta_experiment)) +
    geom_col(show.legend = FALSE) +
    geom_text(aes(y = value / 2, label = format_count(value)), color = "white", size = 3.2) +
    labs(title = title, x = "meta experiment (methX)", y = y_lab) +
    theme_min()
}

count_sty <- function(seqmodi) {
  # Count S/T/Y phospho occurrences in a vector of SequenceModi strings
  s <- stringr::str_count(seqmodi, "S\\(pho\\)|S\\(UniMod:21\\)")
  t <- stringr::str_count(seqmodi, "T\\(pho\\)|T\\(UniMod:21\\)")
  y <- stringr::str_count(seqmodi, "Y\\(pho\\)|Y\\(UniMod:21\\)")
  tibble::tibble(S = s, T = t, Y = y)
}

summarise_sty_counts <- function(qt) {
  if (!nrow(qt) || !"SequenceModi" %in% names(qt)) return(tibble::tibble())
  qt <- qt %>% mutate(meta_experiment = normalize_meta(meta_experiment))
  sty_counts <- count_sty(qt$SequenceModi)
  sty <- dplyr::bind_cols(qt %>% dplyr::select(meta_experiment), sty_counts) %>%
    group_by(meta_experiment) %>%
    summarise(across(dplyr::all_of(c("S", "T", "Y")), ~ sum(.x, na.rm = TRUE)), .groups = "drop")
  tidyr::pivot_longer(sty, cols = c("S", "T", "Y"), names_to = "Residue", values_to = "Count")
}

plot_sty_bar <- function(sty_df, title = "Phosphosite counts (S/T/Y)", y_max = NULL) {
  if (!nrow(sty_df)) return(NULL)
  meta_levels <- if ("meta_experiment" %in% names(sty_df)) order_meta_levels(sty_df$meta_experiment) else NULL
  if (!is.null(meta_levels)) {
    sty_df <- sty_df %>% mutate(meta_experiment = factor(meta_experiment, levels = meta_levels))
  }
  p <- ggplot(sty_df, aes(x = Residue, y = Count, fill = Residue)) +
    geom_col(show.legend = FALSE) +
    geom_text(aes(y = Count / 2, label = format_count(Count)), color = "white", size = 3) +
    labs(title = title, x = "Residue", y = "Counts in SequenceModi") +
    theme_min()
  if (!is.null(meta_levels) && length(meta_levels) > 1) {
    p <- p + facet_wrap(~meta_experiment)
  }
  if (!is.null(y_max) && is.finite(y_max)) {
    p <- p + coord_cartesian(ylim = c(0, y_max * 1.05))
  }
  p
}

plot_ty_bar <- function(sty_df, title = "Phosphosite counts (T/Y only)", y_max = NULL) {
  if (!nrow(sty_df)) return(NULL)
  ty <- sty_df %>% dplyr::filter(Residue %in% c("T", "Y"))
  if (!nrow(ty)) return(NULL)
  meta_levels <- if ("meta_experiment" %in% names(ty)) order_meta_levels(ty$meta_experiment) else NULL
  if (!is.null(meta_levels)) {
    ty <- ty %>% mutate(meta_experiment = factor(meta_experiment, levels = meta_levels))
  }
  p <- ggplot(ty, aes(x = Residue, y = Count, fill = Residue)) +
    geom_col(show.legend = FALSE) +
    geom_text(aes(y = Count / 2, label = format_count(Count)), color = "white", size = 3) +
    labs(title = title, x = "Residue", y = "Counts in SequenceModi") +
    theme_min()
  if (!is.null(meta_levels) && length(meta_levels) > 1) {
    p <- p + facet_wrap(~meta_experiment)
  }
  if (!is.null(y_max) && is.finite(y_max)) {
    p <- p + coord_cartesian(ylim = c(0, y_max * 1.05))
  }
  p
}

# Heatmap (optional) --------------------------------------------------------

plot_intensity_heatmap <- function(qt, top_n = 150) {
  if (!have_complexheatmap) return(NULL)
  cand <- c("PrecursorArea_dstrAdj", "PrecursorArea", "SequenceArea")
  col <- cand[cand %in% names(qt)][1]
  if (is.na(col) || !"SequenceModi" %in% names(qt) || !"SpectrumFile" %in% names(qt)) return(NULL)

  df <- qt %>%
    mutate(val = as.numeric(.data[[col]])) %>%
    group_by(SequenceModi, SpectrumFile) %>%
    summarise(val = max(val, na.rm = TRUE), .groups = "drop") %>%
    group_by(SequenceModi) %>%
    mutate(present = sum(!is.na(val))) %>%
    ungroup() %>%
    group_by(SequenceModi) %>%
    summarise(total = sum(val, na.rm = TRUE), present = max(present), .groups = "drop") %>%
    arrange(desc(present), desc(total)) %>%
    slice_head(n = top_n) %>%
    select(SequenceModi)

  wide <- qt %>%
    mutate(val = as.numeric(.data[[col]])) %>%
    semi_join(df, by = "SequenceModi") %>%
    group_by(SequenceModi, SpectrumFile) %>%
    summarise(val = max(val, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = SpectrumFile, values_from = val)

  mat <- as.matrix(wide[,-1])
  rownames(mat) <- wide$SequenceModi
  mat <- log10(mat + 1)
  mat[is.infinite(mat)] <- NA

  ComplexHeatmap::Heatmap(mat, name = "log10+1", cluster_columns = TRUE, cluster_rows = TRUE)
}

# Main runners --------------------------------------------------------------

detect_meta_roots <- function(which) {
  # If PROC_ROOT contains meth*/<which>, use those; else return the base <which> dir
  meth_tags <- paste0("meth", 1:20)
  typed_candidates <- file.path(PROC_ROOT, meth_tags, which)
  typed_dirs <- typed_candidates[dir_exists(typed_candidates)]
  if (length(typed_dirs)) {
    return(tibble::tibble(
      path = typed_dirs,
      meta_experiment = basename(dirname(typed_dirs))
    ))
  }

  raw_candidates <- file.path(PROC_ROOT, meth_tags)
  raw_dirs <- raw_candidates[dir_exists(raw_candidates)]
  if (length(raw_dirs)) {
    return(tibble::tibble(
      path = raw_dirs,
      meta_experiment = basename(raw_dirs)
    ))
  }

  fallback <- switch(which, phos = PHOS_DIR, prof = PROF_DIR)
  fallback <- unique(fallback[nzchar(fallback)])
  if (!length(fallback)) {
    return(tibble::tibble(path = character(), meta_experiment = character()))
  }
  tibble::tibble(
    path = fallback,
    meta_experiment = normalize_meta(infer_meta_experiment(fallback))
  )
}

analyze_space <- function(root, which = c("phos", "prof")) {
  which <- match.arg(which)
  out_dir <- file.path(OUT_DIR, which)
  ensure_dir(out_dir)

  meta_roots <- detect_meta_roots(which)
  if (!nrow(meta_roots)) {
    log_warn("No meta roots found under %s for %s; using %s directly", PROC_ROOT, which, root)
    meta_roots <- tibble::tibble(path = root, meta_experiment = infer_meta_experiment(root))
  }

  if (nrow(meta_roots) == 1) {
    log_info("Loading %s data from %s (meta_experiment=%s)", which, meta_roots$path[[1]], meta_roots$meta_experiment[[1]])
  } else {
    log_info("Loading %s data from %d meta folders under %s", which, nrow(meta_roots), PROC_ROOT)
    purrr::pwalk(meta_roots, function(path, meta_experiment) {
      log_info(" - %s (meta_experiment=%s)", path, meta_experiment)
    })
  }

  qt_raw <- read_diann_psms_quant(meta_roots)
  qt <- qt_raw
  if (nrow(qt_raw)) {
    qt <- dedupe_psms(qt_raw)
    removed <- nrow(qt_raw) - nrow(qt)
    if (removed > 0) log_info("Deduplicated %d rows (SequenceModi+Charge+RT key) from QUANT", removed)
  }
  ql <- read_diann_psms_qual(meta_roots)
  rp <- read_diann_reports_parquet(root)
  log_info("Loaded: QUANT=%d rows, QUAL=%d rows, Parquet=%s", nrow(qt), nrow(ql), ifelse(is.null(rp), "none", "present"))
  if ("meta_experiment" %in% names(qt)) {
    exps <- sort(unique(normalize_meta(qt$meta_experiment)))
    if (length(exps)) log_info("Meta experiments detected (QUANT): %s", paste(exps, collapse = ", "))
    meta_summary <- qt %>%
      mutate(meta_experiment = normalize_meta(meta_experiment)) %>%
      count(meta_experiment, name = "n_rows")
    purrr::pwalk(meta_summary, function(meta_experiment, n_rows) {
      log_info("  meta_experiment=%s -> %d QUANT rows", meta_experiment, n_rows)
    })
    if ("sample_id" %in% names(qt)) {
      sample_summary <- qt %>%
        mutate(
          meta_experiment = normalize_meta(meta_experiment),
          sample_id = dplyr::coalesce(sample_id, "unknown")
        ) %>%
        count(meta_experiment, sample_id, name = "n_rows")
      purrr::pwalk(sample_summary, function(meta_experiment, sample_id, n_rows) {
        log_info("    sample_id=%s (meta=%s) -> %d QUANT rows", sample_id, meta_experiment, n_rows)
      })
    }
  }
  if (nrow(ql) && "meta_experiment" %in% names(ql)) {
    ql_summary <- ql %>%
      mutate(meta_experiment = normalize_meta(meta_experiment)) %>%
      count(meta_experiment, name = "n_rows")
    purrr::pwalk(ql_summary, function(meta_experiment, n_rows) {
      log_info("  meta_experiment=%s -> %d QUAL rows", meta_experiment, n_rows)
    })
  }

  n_meta_quant <- if (nrow(qt)) length(unique(normalize_meta(qt$meta_experiment))) else 0
  n_meta_qual <- if (nrow(ql)) length(unique(normalize_meta(ql$meta_experiment))) else 0
  width_quant <- calc_plot_width(n_meta_quant, base = 9, extra_per = 1.8)
  width_qual <- calc_plot_width(n_meta_qual, base = 8, extra_per = 1.6)
  qt_phos <- if (which == "phos") filter_phospho_psms(qt) else tibble::tibble()
  qt_phos_res <- if (nrow(qt_phos)) expand_phospho_residues(qt_phos) else tibble::tibble()
  if (nrow(qt_phos)) {
    log_info("Phospho-filtered QUANT rows: %d (SequenceModi with S/T/Y phospho)", nrow(qt_phos))
    phos_summary <- qt_phos %>%
      mutate(meta_experiment = normalize_meta(meta_experiment)) %>%
      count(meta_experiment, name = "n_rows")
    purrr::pwalk(phos_summary, function(meta_experiment, n_rows) {
      log_info("  meta_experiment=%s -> %d phospho-filtered rows", meta_experiment, n_rows)
    })
    if (nrow(qt_phos_res)) {
      res_summary <- qt_phos_res %>%
        mutate(meta_experiment = normalize_meta(meta_experiment)) %>%
        count(meta_experiment, Residue, name = "n_rows")
      purrr::pwalk(res_summary, function(meta_experiment, Residue, n_rows) {
        log_info("    %s phospho -> %d rows (meta=%s)", Residue, n_rows, meta_experiment)
      })
    }
  }

  # Plots: counts
  tryCatch({
    p1 <- plot_psm_counts_per_run(qt, "PSMs per run (deduplicated)")
    if (!is.null(p1)) ggsave_quiet(file.path(out_dir, "counts_psms_per_run_by_meta.png"), p1, width_quant, 4.8)
  }, error = function(e) log_warn("counts_psms_per_run failed: %s", conditionMessage(e)))
  tryCatch({
    p2 <- plot_peptide_counts_per_run(qt, "Unique peptides per run (deduplicated)")
    if (!is.null(p2)) ggsave_quiet(file.path(out_dir, "counts_peptides_per_run_by_meta.png"), p2, width_quant, 4.8)
  }, error = function(e) log_warn("counts_peptides_per_run failed: %s", conditionMessage(e)))
  tryCatch({
    p3 <- plot_gene_counts_per_run(qt, "Unique genes per run")
    if (!is.null(p3)) ggsave_quiet(file.path(out_dir, "counts_genes_per_run_by_meta.png"), p3, width_quant, 4.8)
  }, error = function(e) log_warn("counts_genes_per_run failed: %s", conditionMessage(e)))

  if (nrow(qt_phos)) {
    tryCatch({
      p1p <- plot_psm_counts_per_run(qt_phos, "PSMs per run (S/T/Y phospho, deduplicated)")
      if (!is.null(p1p)) ggsave_quiet(file.path(out_dir, "counts_psms_per_run_phos_by_meta.png"), p1p, width_quant, 4.8)
    }, error = function(e) log_warn("counts_psms_per_run_phos failed: %s", conditionMessage(e)))
    tryCatch({
      p2p <- plot_peptide_counts_per_run(qt_phos, "Unique peptides per run (S/T/Y phospho, deduplicated)")
      if (!is.null(p2p)) ggsave_quiet(file.path(out_dir, "counts_peptides_per_run_phos_by_meta.png"), p2p, width_quant, 4.8)
    }, error = function(e) log_warn("counts_peptides_per_run_phos failed: %s", conditionMessage(e)))
  }

  # Aggregates per meta
  tryCatch({
    pm <- plot_counts_per_meta(qt, "psm", "PSMs by meta experiment (deduplicated)")
    if (!is.null(pm)) ggsave_quiet(file.path(out_dir, "counts_psms_per_meta.png"), pm, 7, 4.2)
  }, error = function(e) log_warn("counts_psms_per_meta failed: %s", conditionMessage(e)))
  tryCatch({
    pm <- plot_counts_per_meta(qt, "peptide", "Unique peptides by meta experiment (deduplicated)")
    if (!is.null(pm)) ggsave_quiet(file.path(out_dir, "counts_peptides_per_meta.png"), pm, 7, 4.2)
  }, error = function(e) log_warn("counts_peptides_per_meta failed: %s", conditionMessage(e)))
  tryCatch({
    pm <- plot_counts_per_meta(qt, "gene", "Unique genes by meta experiment")
    if (!is.null(pm)) ggsave_quiet(file.path(out_dir, "counts_genes_per_meta.png"), pm, 7, 4.2)
  }, error = function(e) log_warn("counts_genes_per_meta failed: %s", conditionMessage(e)))

  if (nrow(qt_phos)) {
    tryCatch({
      pm <- plot_counts_per_meta(qt_phos, "psm", "PSMs by meta experiment (S/T/Y phospho)")
      if (!is.null(pm)) ggsave_quiet(file.path(out_dir, "counts_psms_per_meta_phos.png"), pm, 7, 4.2)
    }, error = function(e) log_warn("counts_psms_per_meta_phos failed: %s", conditionMessage(e)))
    tryCatch({
      pm <- plot_counts_per_meta(qt_phos, "peptide", "Unique peptides by meta experiment (S/T/Y phospho)")
      if (!is.null(pm)) ggsave_quiet(file.path(out_dir, "counts_peptides_per_meta_phos.png"), pm, 7, 4.2)
    }, error = function(e) log_warn("counts_peptides_per_meta_phos failed: %s", conditionMessage(e)))
  }

  # Plots: distributions
  tryCatch({
    p4 <- plot_intensity_hist(qt, "Intensity distribution", "SequenceArea (deduplicated)")
    if (!is.null(p4)) ggsave_quiet(file.path(out_dir, "hist_intensity_by_meta.png"), p4, width_quant, 4.8)
  }, error = function(e) log_warn("hist_intensity failed: %s", conditionMessage(e)))
  tryCatch({
    p5 <- plot_charge_bar(ql, "Charge state counts", "all PSMs")
    if (!is.null(p5)) ggsave_quiet(file.path(out_dir, "bar_charge_states_by_meta.png"), p5, width_qual, 4.2)
  }, error = function(e) log_warn("bar_charge_states failed: %s", conditionMessage(e)))
  tryCatch({
    p6 <- plot_rt_hist(ql, "Retention time distribution", "all PSMs")
    if (!is.null(p6)) ggsave_quiet(file.path(out_dir, "hist_RT_by_meta.png"), p6, width_qual, 4.5)
  }, error = function(e) log_warn("hist_RT failed: %s", conditionMessage(e)))
  tryCatch({
    p7 <- plot_im_hist(ql, "Ion mobility distribution", "all PSMs")
    if (!is.null(p7)) ggsave_quiet(file.path(out_dir, "hist_IM_by_meta.png"), p7, width_qual, 4.5)
  }, error = function(e) log_warn("hist_IM failed: %s", conditionMessage(e)))

  # Phospho-specific
  if (which == "phos") {
    tryCatch({
      sty_df <- summarise_sty_counts(qt)
      if (nrow(sty_df)) {
        sty_width <- calc_plot_width(length(unique(sty_df$meta_experiment)), base = 7, extra_per = 1.3)
        sty_max <- max(sty_df$Count, na.rm = TRUE)
        ty_max <- max(sty_df %>% dplyr::filter(Residue %in% c("T", "Y")) %>% dplyr::pull(Count), na.rm = TRUE)
        p8 <- plot_sty_bar(sty_df, "Phosphosite counts (S/T/Y)", y_max = sty_max)
        if (!is.null(p8)) ggsave_quiet(file.path(out_dir, "bar_STY_counts_by_meta.png"), p8, sty_width, 4.2)
        ty_plot <- plot_ty_bar(sty_df, "Phosphosite counts (T/Y only)", y_max = ty_max)
        if (!is.null(ty_plot)) ggsave_quiet(file.path(out_dir, "bar_TY_counts_by_meta.png"), ty_plot, sty_width, 4.2)
        meta_levels <- unique(sty_df$meta_experiment)
        purrr::walk(meta_levels, function(meta_label) {
          sub_df <- dplyr::filter(sty_df, meta_experiment == meta_label)
          if (!nrow(sub_df)) return(NULL)
          slug <- sanitize_token(meta_label)
          single_plot <- plot_sty_bar(sub_df, sprintf("Phosphosite counts (S/T/Y) - %s", meta_label), y_max = sty_max)
          if (!is.null(single_plot)) ggsave_quiet(file.path(out_dir, sprintf("bar_STY_counts_%s.png", slug)), single_plot, 5.5, 4)
          ty_single <- plot_ty_bar(sub_df, sprintf("Phosphosite counts (T/Y) - %s", meta_label), y_max = ty_max)
          if (!is.null(ty_single)) ggsave_quiet(file.path(out_dir, sprintf("bar_TY_counts_%s.png", slug)), ty_single, 5.5, 4)
        })
      }
    }, error = function(e) log_warn("phospho residue plots failed: %s", conditionMessage(e)))
    tryCatch({
      ip <- plot_intensity_hist(qt_phos, "Intensity distribution", "SequenceArea (S/T/Y phospho)")
      if (!is.null(ip)) ggsave_quiet(file.path(out_dir, "hist_intensity_phos_by_meta.png"), ip, width_quant, 4.8)
    }, error = function(e) log_warn("hist_intensity_phos failed: %s", conditionMessage(e)))
    if (nrow(qt_phos_res)) {
      tryCatch({
        charge_ty <- plot_charge_bar(qt_phos_res %>% dplyr::filter(Residue %in% c("T", "Y")),
                                     "Charge state counts (T/Y phospho)", "filtered")
        if (!is.null(charge_ty)) ggsave_quiet(file.path(out_dir, "bar_charge_states_TY_phos_by_meta.png"), charge_ty, width_qual, 4.2)
      }, error = function(e) log_warn("bar_charge_states_TY_phos failed: %s", conditionMessage(e)))
      tryCatch({
        charge_res <- plot_charge_bar_residue(qt_phos_res, "Charge state counts by residue (S/T/Y phospho)")
        if (!is.null(charge_res)) ggsave_quiet(file.path(out_dir, "bar_charge_states_by_residue_phos.png"), charge_res, 7.5, 4.5)
      }, error = function(e) log_warn("bar_charge_states_residue_phos failed: %s", conditionMessage(e)))
      seq_col <- "SequenceArea"
      if (seq_col %in% names(qt_phos_res)) {
        hist_res <- plot_hist_residue(qt_phos_res %>% mutate(log_seq_area = log10p1(.data[[seq_col]])),
                                      "log_seq_area", "Intensity by residue (log10(1+SequenceArea))")
        if (!is.null(hist_res)) ggsave_quiet(file.path(out_dir, "hist_intensity_by_residue_phos.png"), hist_res, 8, 4.8)
      }
      rt_cols <- intersect(c("RTmin", "RT", "Retention.Time"), names(qt_phos_res))
      if (length(rt_cols)) {
        hist_rt_res <- plot_hist_residue(qt_phos_res, rt_cols[1], paste0("Retention time by residue (", rt_cols[1], ")"))
        if (!is.null(hist_rt_res)) ggsave_quiet(file.path(out_dir, "hist_RT_by_residue_phos.png"), hist_rt_res, 8, 4.5)
      }
      im_cols <- intersect(c("IM", "Ion.Mobility"), names(qt_phos_res))
      if (length(im_cols)) {
        hist_im_res <- plot_hist_residue(qt_phos_res, im_cols[1], paste0("Ion mobility by residue (", im_cols[1], ")"))
        if (!is.null(hist_im_res)) ggsave_quiet(file.path(out_dir, "hist_IM_by_residue_phos.png"), hist_im_res, 8, 4.5)
      }
    }
  }

  if (which == "phos" && file.exists(RAW_SUMMARY_PATH)) {
    tryCatch({
      raw_summary <- readr::read_csv(RAW_SUMMARY_PATH, show_col_types = FALSE) %>%
        mutate(
          meta_experiment = normalize_meta(stringr::str_extract(raw_path, "meth\\d+")),
          raw_file = basename(raw_path),
          msms_fraction = dplyr::if_else(total_frames > 0, msms_frames / total_frames, NA_real_),
          rt_span_min = (rt_end - rt_start) / 60,
          rt_start_min = rt_start / 60,
          rt_end_min = rt_end / 60
        )

      raw_width <- calc_plot_width(length(unique(raw_summary$meta_experiment)), base = 8, extra_per = 1.0)

      p_raw_frames <- raw_summary %>%
        ggplot(aes(x = forcats::fct_reorder(raw_file, total_frames), y = total_frames, fill = meta_experiment)) +
        geom_col(show.legend = FALSE) +
        coord_flip() +
        labs(title = "Raw frames per acquisition", x = "Raw file", y = "Frames") +
        theme_min()
      ggsave_quiet(file.path(out_dir, "raw_frames_per_run.png"), p_raw_frames, raw_width, 4.5)

      p_raw_rt <- raw_summary %>%
        ggplot(aes(x = forcats::fct_reorder(raw_file, rt_start_min), ymin = rt_start_min, ymax = rt_end_min, color = meta_experiment)) +
        geom_linerange(linewidth = 1.1) +
        coord_flip() +
        labs(title = "Raw retention time coverage", x = "Raw file", y = "Retention time (min)") +
        theme_min()
      ggsave_quiet(file.path(out_dir, "raw_rt_range.png"), p_raw_rt, raw_width, 4.5)

      raw_dia <- raw_summary %>% dplyr::filter(!is.na(dia_total_windows))
      if (nrow(raw_dia)) {
        p_raw_dia <- raw_dia %>%
          ggplot(aes(x = forcats::fct_reorder(raw_file, dia_total_windows), y = dia_total_windows, fill = meta_experiment)) +
          geom_col(show.legend = FALSE) +
          coord_flip() +
          labs(title = "DIA windows per acquisition", x = "Raw file", y = "Total DIA windows") +
          theme_min()
        ggsave_quiet(file.path(out_dir, "raw_dia_windows.png"), p_raw_dia, raw_width, 4.2)
      }
    }, error = function(e) log_warn("raw summary integration failed: %s", conditionMessage(e)))
  }

  # Optional heatmap
  if (have_complexheatmap) {
    tryCatch({
      hm <- plot_intensity_heatmap(qt, top_n = 120)
      if (!is.null(hm)) {
        f <- file.path(out_dir, "heatmap_peptide_intensity.pdf")
      pdf(f, width = 8, height = 10)
      grid::grid.newpage(); ComplexHeatmap::draw(hm); dev.off()
      log_info("Wrote %s", f)
      }
    }, error = function(e) log_warn("heatmap failed: %s", conditionMessage(e)))
  }

  invisible(list(qt = qt, ql = ql, rp = rp))
}

main <- function() {
  raw_args <- commandArgs(trailingOnly = TRUE)
  args <- tolower(raw_args)
  modes <- c("phos", "prof", "all")
  mode <- ifelse(length(args) >= 1 && args[1] %in% modes, args[1], "all")

  if (mode %in% c("phos", "all")) {
    log_info("Analyzing phospho at %s", PHOS_DIR)
    tryCatch({
      analyze_space(PHOS_DIR, "phos")
      log_info("Phospho plots in %s", file.path(OUT_DIR, "phos"))
    }, error = function(e) {
      log_error("Phospho analysis failed: %s", conditionMessage(e))
    })
  }
  if (mode %in% c("prof", "all")) {
    log_info("Analyzing proteome at %s", PROF_DIR)
    tryCatch({
      analyze_space(PROF_DIR, "prof")
      log_info("Proteome plots in %s", file.path(OUT_DIR, "prof"))
    }, error = function(e) {
      log_error("Proteome analysis failed: %s", conditionMessage(e))
    })
  }
  log_info("Done.")
}

if (identical(environment(), globalenv())) {
  # Run and never hard-crash; log errors and exit 0 to not block pipelines
  tryCatch(main(), error = function(e) {
    log_error("Unexpected error: %s", conditionMessage(e))
  })
}
