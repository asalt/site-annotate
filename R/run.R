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


  if (exists("get_tool_env", where = .GlobalEnv, mode = "function")){
    get_tool_env <- .GlobalEnv$get_tool_env
  }  # TRUE

  io_tools <- get_tool_env("io")
  heatmap_tools <- get_tool_env("heatmap")
  util_tools <- get_tool_env("utils")
  limma_tools <- get_tool_env("limma")

  OUTDIR <- output_dir %||% data_dir %||% "."
  config_file <- config %||% file.path(here(), "config", "base.toml")
  config <- io_tools$load_config(config_file, root_dir = data_dir)

  print(config)

  # gct_file <- file.path(data_dir %||% '.', config$gct_file) %||% gct_file %||% NULL
  gct_file <- config$gct_file %||% gct_file %||% NULL

  if (!is.null(gct_file)) {
    gct <- cmapR::parse_gctx(gct_file)
  } else {
    cat("no data passed, exiting")
    stop()
  }

  if (!is.null(config$extra$sample_exclude)) {
    to_exclude <- config$extra$sample_exclude %||% NULL
    to_exclude <- intersect(rownames(gct@cdesc), to_exclude)
    message(paste0("excluding ", length(to_exclude), " samples"))
    to_keep <- setdiff(rownames(gct@cdesc), to_exclude)
    gct <- gct %>% subset_gct(cid = to_keep)
  }

  suboutdir <- basename(gct_file) %>%
    basename() %>%
    tools::file_path_sans_ext()

  OUTDIR <- file.path(OUTDIR, suboutdir)
  if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)
  OUTNAME <- "" # don't set an actual name for the internal files saved within the main dir

  FRIENDLY_NAME <- suboutdir %>% str_replace_all("_", " ") %>% str_replace(., "\\d+$", "") %>% str_wrap(width = 60)

  # heatmap of everything
  if (config$norm$normalize){
    gct <- util_tools$normalize_gct(gct, log_transform=config$norm$log_transform %||% FALSE)
    }
  gct_z <- util_tools$scale_gct(gct, group_by = config$heatmap$zscore_by)




  if (config$heatmap$selection$do %||% FALSE) {

    param_grid <- expand.grid(
      cut_by = config$heatmap$selection$cut_by %||% config$heatmap$cut_by %||% NA,
      cluster_rows = config$heatmap$selection$cluster_rows %||% config$heatmap$cluster_rows %||% TRUE,
      cluster_columns = config$heatmap$selection$cluster_columns %||% config$heatmap$cluster_columns %||% c(TRUE, FALSE),
      stringsAsFactors = FALSE
    )


    config$heatmap$selection$file_list %>%
      flatten %>%
      imap(~{
        the_file <- .x
        file_name <- .y

        file_data <- io_tools$read_genelist_file(the_file)

        subgct_z <- gct_z %>% subset_gct(rid = file_data$id)
        print(paste('nrow submat gct: ', nrow(subgct_z@mat)))
        if (nrow(subgct_z@mat) == 0) {
          message(paste0("No rows found in ", the_file, ", skipping"))
          return()
        }
        if (nrow(subgct_z@mat) > 199) {
          message(paste("too many rows", the_file, ", skipping"))
          return()
        }
        .nrow <- nrow(subgct_z@mat)
        .ncol <- ncol(subgct_z@mat)
        #
        param_grid %>% purrr::pmap(~ {
          cut_by_val <- ..1
          cluster_rows_val <- ..2
          cluster_columns_val <- ..3

          .outdir <- file.path(OUTDIR, "heatmap", "selection", paste0("cut_", make.names(cut_by_val)))
          if (!dir.exists(.outdir)) dir.create(.outdir, recursive = TRUE)
          .outf <- file.path(
            .outdir,
            make.names(
              paste0(
                OUTNAME,
                file_name,
                "_",
                .nrow, "x", .ncol,
                "_cut_", cut_by_val,
                "cr_", cluster_rows_val,
                "cc_", cluster_columns_val,
                ".pdf"
              )
            )
          )

          ht_draw_code <- heatmap_tools$make_heatmap_fromgct(
            gct = subgct_z,
            column_title = file_name,
            cut_by = cut_by_val,
            cluster_rows = cluster_rows_val,
            cluster_columns = cluster_columns_val,
            meta_to_include = config$heatmap$legend_include %||% NULL,
            meta_to_exclude = config$heatmap$legend_exclude %||% NULL,
            show_row_names = TRUE,
            include_row_annotations = TRUE,
            optimal_size = TRUE
            )
          height <- 6 + (.nrow * .20)
          width <- 12 + (.ncol * .26)

          message(paste0('drawing ', .outf))
          if (file.exists(.outf) && config$advanced$replace == FALSE){
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
        })
      }) # end of imap of file list
  } # end of if selection


  param_grid <- expand.grid(
    cut_by = unique(c(config$heatmap$cut_by, NA)),
    cluster_columns = c(TRUE, FALSE),
    stringsAsFactors = FALSE
  )

  if (config$heatmap$do %||% FALSE) {
    param_grid %>% purrr::pmap(~ {
      .cut_by <- ..1
      .cluster_columns <- ..2
      ht_draw_code <- gct_z %>% heatmap_tools$make_heatmap_fromgct(
        show_row_names = FALSE,
        optimal_size = FALSE,
        cut_by = .cut_by,
        cluster_columns = .cluster_columns,
        include_row_annotations = FALSE,
        column_title=FRIENDLY_NAME
      )

      .nrow <- nrow(gct_z@mat)
      .ncol <- ncol(gct_z@mat)
      width <- 6 + (.ncol * .26)
      height <- 7.8


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
      .contrast <- str_replace_all(.contrast, "-", "minus")

      .outdir <- file.path(OUTDIR, "limma", "tables")
      .outname <- file.path(.outdir, make.names(paste0(OUTNAME, .contrast, ".tsv")))

      if (!dir.exists(.outdir)) dir.create(.outdir, recursive = TRUE)
      if (!file.exists(.outname) || config$advanced$replace == TRUE) {
        .table %>% write_tsv(.outname)
        message(paste0("wrote ", .outname))
      }

      .outname <- file.path(.outdir, make.names(paste0(OUTNAME, .contrast, "_pdist.png")))
      if (!file.exists(.outname) || config$advanced$replace == TRUE) {
        png(.outname, width = 7, height = 4.5, units = "in", res = 300)
        limma_tools$plot_p_dist(.table, main_title = .contrast)
        dev.off()
        #

        message(paste0("wrote ", .outname))
      }
    })

    limma_topn <- config$limma$topn %||% c(50, 100)
    param_grid <- expand.grid(
      cut_by = unique(c(config$heatmap$cut_by, NA)),
      cluster_rows = c(TRUE, FALSE), # unique(c(config$heatmap$cluster_rows, FALSE)) %||% c(TRUE, FALSE),
      cluster_columns = unique(config$heatmap$cluster_columns %||% FALSE),
      top_n = config$limma$top_n %||% c(50, 100)
      # file_data = file_data
    )

    gct_z <- util_tools$scale_gct(gct, group_by = config$heatmap$zscore_by)

    topTables %>% purrr::imap(~ { # second loop plot topdiff
      .table <- .x %>% arrange(desc(abs(t)))
      .contrast <- .y
      .contrast <- str_replace_all(.contrast, "-", "minus")

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

        .nrow <- nrow(.gct_z@mat)
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
        # this all needs to be in a trycatch
        tryCatch({
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
        }, error = function(e) print(e))
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
     purrr::pmap( ~ {
      do_scale <- ..1
      colby <- ..2
      ntopLoadings <- ..3
      markby <- ..4
      rids_to_keep <- gct@rdesc %>% distinct(fifteenmer, .keep_all = TRUE) %>% pull(id)
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
            ifelse(!is.null(ntopLoadings) && ntopLoadings > 0, paste0('_nload_', ntopLoadings), ""),
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


      # not work yet
      # pca_obj %>% pca_tools$make_heatmap_from_loadings(gct_obj=gct_obj,)

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
} # end of run
