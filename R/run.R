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
    tools::file_path_sans_ext()

  OUTDIR <- file.path(OUTDIR, OUTNAME)
  if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)


  # heatmap of everything
  gct_z <- util_tools$scale_gct(gct, group_by = config$heatmap$zscore_by)


  param_grid <- expand.grid(
    cut_by = unique(c(config$heatmap$cut_by, NA)),
    cluster_columns = c(TRUE, FALSE),
    stringsAsFactors = FALSE
  )

  if (config$heatmap$do) {
    param_grid %>% purrr::pmap(~ {
      .cut_by <- ..1
      .cluster_columns <- ..2
      ht_draw_code <- gct_z %>% heatmap_tools$make_heatmap_fromgct(
        show_row_names = FALSE,
        optimal_size = FALSE,
        cut_by = .cut_by,
        cluster_columns = .cluster_columns
      )

      .nrow <- nrow(gct_z@mat)
      .ncol <- ncol(gct_z@mat)
      width <- 6 + (.ncol * .26)
      height <- 12


      .outdir <- file.path(OUTDIR, "heatmap")
      if (!dir.exists(.outdir)) dir.create(.outdir, recursive = TRUE)
      .outf <- file.path(
        .outdir,
        make.names(
          paste0(
            OUTNAME,
            .nrow, "x", .ncol,
            "_cut_", .cut_by,
            "cc_", .cluster_columns,
            ".pdf"
          )
        )
      )

      if (file.exists(.outf) && config$advanced$replace == FALSE) {
        message(paste0(.outf, " already exists, skipping"))
        return(NULL)
      }
      message(paste0("drawing ", .outf))
      tryCatch(
        {
          cairo_pdf(.outf, width = width, height = height)
          ht_draw_code()
          dev.off()
        },
        error = function(e) print(e)
      )
    }) # end of pmap of param grid
  } # end of if heatmap


  if (config$limma$do) {
    topTables <- limma_tools$run_limma(gct, config)

    topTables %>% purrr::imap(~ { # first loop save tables
      .table <- .x %>% arrange(desc(t))
      .contrast <- .y

      .outdir <- file.path(OUTDIR, "limma", "tables")
      .outname <- file.path(.outdir, make.names(paste0(OUTNAME, .contrast, ".tsv")))

      if (!dir.exists(.outdir)) dir.create(.outdir, recursive = TRUE)
      if (!file.exists(.outname) || config$advanced$replace == TRUE) {
        .table %>% write_tsv(.outname)
      }

      .outname <- file.path(.outdir, make.names(paste0(OUTNAME, .contrast, "_pdist.png")))
      png(.outname, width = 7, height = 4.5, units = "in", res = 300)
      limma_tools$plot_p_dist(.table, main_title = .contrast)
      dev.off()
    })

    limma_topn <- config$limma$topn %||% c(50, 100)
    param_grid <- expand.grid(
      cut_by = unique(c(config$heatmap$cut_by, NULL)),
      cluster_rows = c(TRUE, FALSE), # unique(c(config$heatmap$cluster_rows, FALSE)) %||% c(TRUE, FALSE),
      cluster_columns = config$heatmap$cluster_columns %||% FALSE,
      top_n = config$limma$top_n %||% c(50, 100)
      # file_data = file_data
    )

    gct_z <- util_tools$scale_gct(gct, group_by = config$heatmap$zscore_by)

    topTables %>% purrr::imap(~ { # second loop plot topdiff
      .table <- .x %>% arrange(desc(abs(t)))
      .contrast <- .y

      .outdir <- file.path(OUTDIR, "limma", "tables")



      column_title <- paste0(OUTNAME, "\ncontrast: ", .contrast)

      param_grid %>% purrr::pmap(~ {
        cut_by_val <- ..1
        cluster_rows_val <- ..2
        cluster_columns_val <- ..3
        top_n_val <- ..4

        .top_tvals <- .table %>%
          distinct(t) %>%
          head(n = top_n_val) # %>% pull(t)
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


  if (config$pca$do) {
    pca_tools <- get_tool_env("pca")

    param_grid <- expand.grid(
      scale = c(TRUE, FALSE),
      colby = config$pca$colby %||% config$pca$color %||% NA,
      stringsAsFactors = FALSE
    )

    param_grid %>%
     purrr::pmap( ~ {
      do_scale <- ..1
      colby <- ..2
      pca_obj <- pca_tools$do_pca(gct, scale = do_scale)

      .outdir <- file.path(OUTDIR, "pca")

      .outf <- file.path(
        .outdir,
        make.names(
          paste0(
            "pca",
            "_scale_", do_scale,
            ".pdf"
          )
        )
      )

      plts <- pca_tools$plot_biplot(pca_obj,
        top_pc = 3,
        showLoadings = T,
        labSize = 1.8, pointSize = 3,
        sizeLoadingsNames = 2,
        colby = colby,
        shape = NULL,
        encircle = !is.null(colby),
        title = "",
        basename = .outf
      )

    }) # end of map
  # browser()
}
} # end of run
