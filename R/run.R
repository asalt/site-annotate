# run.R
# main entry point for all downstream analyses

source(file.path("setup_environment.R"))
source(file.path("lazyloader.R"))

cache_tools <- get_tool_env("cache")
`%cached_by%` <- cache_tools$`%cached_by%`



process_gct_file <- function(gct_file, config) {
  util_tools <- get_tool_env("utils")

  if (!file.exists(gct_file)) {
    stop(paste0("File ", gct_file, " does not exist, exiting"))
  }

  # Step 1: Load or Parse GCT File
  # gct <- load_or_parse_gct(gct_file)
  gct <- parse_gctx(gct_file) %cached_by% rlang::hash_file(gct_file)

  # if ("species" %in% names(config)){
  # }

  # Step 2: Normalize GCT
  norm_config <- config$norm
  normed_gct <- util_tools$normalize_gct(gct,
    mednorm = norm_config$mednorm %||% TRUE,
    log_transform = norm_config$log_transform %||% TRUE
  ) %cached_by%
    rlang::hash(c(gct@mat, norm_config))

  # Step 3: Filter Non-Zeros
  orig_dim <- dim(normed_gct@mat)
  filtered_gct <- util_tools$filter_nonzeros(
    normed_gct, config$filter$non_zeros %||% 1.0,
    config$filter$nonzero_subgroup %||% NULL
  ) %cached_by% rlang::hash(c(normed_gct@mat, config$filter$non_zeros %||% 1.0, config$filter$nonzero_subgroup %||% NULL))

  result_dim <- dim(filtered_gct@mat)
  message(paste0("dim ", paste0(as.character(orig_dim), collapse = "x"), " -> ", paste0(as.character(result_dim), collapse = "x")))


  # Step 4: Apply Batch Correction (if required)
  # need more error checking ehre or inside
  final_gct <- filtered_gct
  if (is.character(config$norm$batch) %||% FALSE) {
    final_gct <- util_tools$do_combat(filtered_gct, by = config$norm$batch) %cached_by% rlang::hash(c(filtered_gct@mat, config$norm$batch))
    # final_gct <- util_tools$do_limmaremovebatch(filtered_gct, by = config$norm$batch) %cached_by% rlang::hash(c(filtered_gct@mat, config$norm$batch))
  }

  return(final_gct)
}


exclude_samples <- function(gct, sample_exclude) {
  if (!is.null(sample_exclude)) {
    to_exclude <- intersect(rownames(gct@cdesc), sample_exclude)
    message(paste0("Excluding ", length(to_exclude), " samples"))
    to_keep <- setdiff(rownames(gct@cdesc), to_exclude)
    gct <- gct %>% subset_gct(cid = to_keep)
  }
  return(gct)
}

load_and_validate_config <- function(config_file, data_dir) {
  io_tools <- get_tool_env("io")
  config <- io_tools$load_config(config_file, root_dir = data_dir)
  if (is.null(config)) stop("Config file could not be loaded.")
  return(config)
}


# Main function to handle heatmap generation
generate_global_heatmaps <- function(gct_z, config, outdir, friendly_name, outname = NULL) {
  param_grid <- expand.grid(
    cut_by = unique(c(config$heatmap$cut_by %||% NA, NA)),
    cluster_columns = c(TRUE, FALSE),
    stringsAsFactors = FALSE
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

  # Prepare file paths
  outf <- heatmap_tools$prepare_heatmap_output_path(outdir,
    "heatmap",
    outname,
    prefix = glue::glue("X{nrow(gct_z_subset@mat)}x{ncol(gct_z_subset@mat)}"),
    # prefix = list(nrow=nrow(gct_z_subset@mat), ncol(gct_z_subset@mat)),
    cut_by = cut_by,
    cluster_rows = T,
    cluster_columns = cluster_columns
  )

  # Skip if file exists and replace is FALSE
  if (file.exists(outf) && config$advanced$replace == FALSE) {
    message(paste0(outf, " already exists, skipping"))
    return(NULL)
  }

  # Draw and save the heatmap
  heatmap_tools$draw_and_save_heatmap(ht_draw_code, outf, gct_z_subset, width = NULL, height = 11)
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

  file_list <- config$heatmap$selection$file_list
  if (is.null(file_list)) {
    message("No file list provided for selection heatmaps, skipping")
    return(NULL)
  }

  config$heatmap$selection$file_list %>%
    purrr::list_flatten() %>%
    purrr::imap(~ {
      the_file <- .x
      file_name <- .y

      # Read and subset GCT
      file_data <- io_tools$read_genelist_file(the_file)
      subgct_z <- gct_z %>% subset_gct(rid = file_data$id %||% file_data[[1]])

      quantiles <- quantile(
        probs = seq(0, 1, 0.025),
        na.rm = T
      )
      heatmap_colorbar_bounds <- c(quantiles[["2.5%"]], quantiles[["97.5%"]])

      if (nrow(subgct_z@mat) == 0) {
        message(paste0("No rows found in ", the_file, ", skipping"))
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

  if (file.exists(outf) && config$advanced$replace == FALSE) {
    message(paste0(outf, " already exists, skipping"))
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

  print(config)

  z_score_flag <- config$heatmap$z_score %||% TRUE
  z_score_by <- config$heatmap$z_score_by %||% config$heatmap$zscore_by %||% NULL


  if (is.null(gct_file)) {
    gct_file <- config$gct_file
  }
  gct <- process_gct_file(gct_file, config)
  gct <- exclude_samples(gct, config$extra$sample_exclude)

  .outdir <- file.path(OUTDIR, "export")
  if (!dir.exists(.outdir)) dir.create(.outdir, recursive = TRUE)
  gct %>% write_gct(file.path(.outdir, "export"))

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
      .contrast <- .y
      .contrast <- str_replace_all(.contrast, "-", "minus")

      .outdir <- file.path(OUTDIR, "limma", "tables")

      column_title <- paste0(OUTNAME, "\ncontrast: ", .contrast)

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
        tryCatch(
          {
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
} # end of run
