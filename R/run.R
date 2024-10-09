# main entry point for all downstream analyses

run <- function(
    data_dir = NULL,
    output_dir = NULL,
    metadata = NULL,
    config_file = NULL,
    gct_file = NULL,
    save_env = FALSE,
    here_dir = NULL,
    ...) {
  knitr::opts_chunk$set(echo = TRUE)
  suppressPackageStartupMessages(library(rlang))
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(magrittr))
  suppressPackageStartupMessages(library(purrr))
  suppressPackageStartupMessages(library(ComplexHeatmap))
  suppressPackageStartupMessages(library(circlize))
  suppressPackageStartupMessages(library(reactable))
  suppressPackageStartupMessages(library(RColorBrewer))
  suppressPackageStartupMessages(library(here))
  library(dendsort)
  print(paste0("here is defined as: ", here()))

  source(file.path("lazyloader.R"))


  io_tools <- get_tool_env("io")
  heatmap_tools <- get_tool_env("heatmap")
  util_tools <- get_tool_env("utils")
  limma_tools <- get_tool_env("limma")

  OUTDIR <- output_dir %||% data_dir %||% "."
  config_file <- config %||% file.path(here(), "config", "base.toml")
  config <- io_tools$load_config(config_file)

  print(config)


  if (!is.null(gct_file)) {
    gct <- cmapR::parse_gctx(gct_file)
  } else {
    cat("no data passed, exiting")
    stop()
  }


  OUTNAME <- basename(gct_file) %>%
    basename() %>%
    gsub("_n\\d+x\\d+\\.gct", "", .)

  if (config$limma$do) {
    topTables <- limma_tools$run_limma(gct, config)
    limma_topn <- config$limma$topn %||% c(50, 100)


    param_grid <- expand.grid(
      cut_by = config$heatmap$cut_by %||% NULL,
      cluster_rows = config$heatmap$cluster_rows %||% TRUE,
      cluster_columns = config$heatmap$cluster_columns %||% FALSE,
      top_n = config$limma$top_n %||% c(50, 100)
      # file_data = file_data
    )

    gct_z <- util_tools$scale_gct(gct, group_by = config$heatmap$zscore_by)


    topTables %>% purrr::imap(~ {
      .table <- .x %>% arrange(desc(abs(t)))
      .contrast <- .y

      .outdir <- file.path(OUTDIR, "limma", "tables")
      .outname <- file.path(.outdir, make.names(paste0(OUTNAME, .contrast, ".tsv")))
      if (!dir.exists(.outdir)) dir.create(.outdir, recursive = TRUE)
      if (!file.exists(.outname) || config$advanced$replace == TRUE) {
        .table %>% write_tsv(.outname)
      }

      column_title <- paste0(OUTNAME, "\ncontrast: ", .contrast)

      param_grid %>% purrr::pmap(~ {
        cut_by_val <- ..1
        cluster_rows_val <- ..2
        cluster_columns_val <- ..3
        top_n_val <- ..4

        .top_tvals <- .table %>%
          distinct(t) %>%
          head(n = top_n_val) # %>% pull(t)
        .subsel <- .table %>% filter(t %in% .top_tvals$t)

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

        .nrow <- nrow(.top_tvals)
        .ncol <- ncol(.gct_z@mat)
        height <- 6 + (.nrow * .20)
        width <- 12 + (.ncol * .26)

        .outdir <- file.path(OUTDIR, "limma", "all")
        if (!dir.exists(.outdir)) dir.create(.outdir, recursive = TRUE)
        .outf <- file.path(
          .outdir,
          make.names(
            paste0(
              OUTNAME,
              .contrast, "_",
              .nrow, "x", .ncol,
              "_cut_", cut_by_val,
              "cr_", cluster_rows_val,
              "cc_", cluster_columns_val,
              ".pdf"
            )
          )
        )

        message(paste0("drawing ", .outf))
        if (file.exists(.outf) && config$advanced$replace == FALSE) {
          message(paste0(.outf, " already exists, skipping"))
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
        column_names <- names(gct@cdesc)

        values <- stringr::str_extract_all(.contrast, "(?<=- )[A-Za-z0-9_]+")[[1]]
        matched_column <- column_names[purrr::map_lgl(column_names, ~ str_detect(.contrast, .x))]
        if (length(matched_column) == 0) {
          warning(paste0("No column found matching contrast: ", .contrast))
          return()
        }
        col_name <- matched_column[1] # Assuming there's one match; adjust if multiple matches needed
        cleaned_contrast <- str_replace_all(.contrast, col_name, "")
        group_values <- str_split(cleaned_contrast, "\\s?-\\s?") %>%
          unlist() %>%
          str_trim()

        .cdesc <- gct@cdesc %>% dplyr::filter(!!sym(col_name) %in% group_values)
        sub_gct_z <- gct %>%
          subset_gct(cid = rownames(.cdesc)) %>%
          util_tools$scale_gct() # %>%
        .sub_gct_z <- sub_gct_z %>% subset_gct(rid = .subsel$id)


        .nrow <- nrow(.top_tvals)
        .ncol <- ncol(.sub_gct_z@mat)

        .outdir <- file.path(OUTDIR, "limma", "subsets")
        if (!dir.exists(.outdir)) dir.create(.outdir, recursive = TRUE)
        .outf <- file.path(
          .outdir,
          make.names(
            paste0(
              OUTNAME,
              .contrast, "_",
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

        message(paste0("drawing ", .outf))
        if (file.exists(.outf) && config$advanced$replace == FALSE) {
          message(paste0(.outf, " already exists, skipping"))
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
      }) # end of pmap of param grid
    }) # end of imap of topTables
  } # end of if limma
}
