# heatmap.R
# A collection of utility functions for data processing and visualization
# Author: Alex Saltzman
# Date: 2024-

# Load necessary libraries within the functions or ensure they are loaded in the main script
# to avoid loading them multiple times.

# this file has gotten quite large with a variety of helper functions

library(here)

source(file.path(here(), "R", "lazyloader.R"))
util_tools <- get_tool_env("utils")

hide_zero <- function(x) {
  sapply(x, function(y) ifelse(y == 0, "", as.character(y)))
}


# Function: sort_by_number_after_underscore ( is this being used? )
# Description: Sorts a vector of strings based on the numeric part that follows the first underscore.
# Parameters:
#   strings (character vector): Vector of strings to be sorted.
# Returns:
#   A character vector sorted based on the numeric value after the first underscore.
#   If no numeric part is found after the underscore, those strings are placed at the end.
sort_by_number_after_underscore <- function(strings) {
  if (!is.character(strings)) {
    stop("Input 'strings' must be a character vector.")
  }

  # Extract the numeric part after the first underscore
  numeric_part <- as.numeric(gsub(".*?_(\\d+).*", "\\1", strings))

  # Identify strings with valid numeric parts
  valid_numeric <- !is.na(numeric_part)

  # Sort based on numeric_part where available
  sorted_indices <- order(
    ifelse(valid_numeric, numeric_part, Inf)
  )

  sorted_strings <- strings[sorted_indices]

  return(sorted_strings)
}

# Function: condense_groups
# Description: Condenses groups by aggregating non-numeric columns and summarizing numeric columns.
# Parameters:
#   .mat (data.frame or tibble): Data frame containing both numeric and non-numeric columns.
#     Expected Columns:
#       - 'sitename' (character): Name of the site.
#       - 'ms_lit', 'lt_lit' (numeric): Specific numeric columns to exclude from general numeric processing.
#       - Other numeric and non-numeric columns as required.
# Returns:
#   A condensed data frame with aggregated non-numeric columns and summarized numeric columns.
#   The 'sitename' column is formatted for better readability.
#   If expected columns are missing, the function will throw an error.
condense_groups <- function(.mat) {
  # Load required libraries
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required but not installed.")
  }
  if (!requireNamespace("stringr", quietly = TRUE)) {
    stop("Package 'stringr' is required but not installed.")
  }

  library(dplyr)
  library(stringr)

  # Check if .mat is a data frame
  if (!is.data.frame(.mat)) {
    stop("Input '.mat' must be a data frame or tibble.")
  }

  # Check for required columns
  required_non_numeric_exclusions <- c("ms_lit", "lt_lit")
  missing_exclusions <- setdiff(required_non_numeric_exclusions, names(.mat))
  if (length(missing_exclusions) > 0) {
    warning(paste("Missing required columns:", paste(missing_exclusions, collapse = ", ")))
  }

  # Identify numeric and non-numeric columns excluding specific columns
  numeric_cols <- sapply(.mat, is.numeric) & !names(.mat) %in% required_non_numeric_exclusions
  non_numeric_cols <- !numeric_cols

  # Select columns accordingly
  selected_numeric <- names(.mat)[numeric_cols]
  selected_non_numeric <- names(.mat)[non_numeric_cols]

  # Group by all numeric columns and summarize non-numeric columns
  condensed_data <- .mat %>%
    select(all_of(c(selected_non_numeric, selected_numeric))) %>%
    group_by(across(all_of(selected_numeric))) %>%
    summarize(
      across(
        all_of(selected_non_numeric),
        ~ paste(sort(unique(.)), collapse = "/"),
        .names = "aggregated_{col}"
      ),
      .groups = 'drop'
    )

  # Format 'sitename' for better readability, if present
  if ("aggregated_sitename" %in% names(condensed_data)) {
    condensed_data <- condensed_data %>%
      mutate(
        aggregated_sitename = str_replace_all(aggregated_sitename, '/', ' / ') %>%
          str_wrap(width = 48)
      )
  }

  return(condensed_data)
}

# Function: do_plot
# Description: Generates and saves heatmaps based on differential expression results.
# Parameters:
#   df (data.frame or tibble): Data frame containing differential expression results.
#     Expected Columns:
#       - 'P.Value' (numeric): P-values from the analysis.
#       - 'site_id' (character): Unique identifier for each site.
#       - Other necessary columns as used within the function.
#   comp_name (character): Name of the comparison for labeling the plot.
#   n (integer, default = 60): Number of top sites to include in the plot.
#   datal (data.frame or tibble): Data frame containing normalized data with 'site_id', 'name', 'zscore', etc.
#   meta (data.frame or tibble): Metadata data frame containing sample information with 'name', 'group', etc.
#   OUTDIR (character): Output directory path where the plots will be saved.
#   util_tools (list): A list containing utility functions, e.g., 'create_named_color_list'.
# Returns:
#   Saves heatmap PDFs to the specified output directory.
#   Throws errors if required packages are missing or input data is invalid.

.do_plot_example <- function(df, comp_name, n = 60, datal, meta, OUTDIR, util_tools) {
    # this is to be refactored. df is a dataframe of stat results. we should have a gct object to make the data loading easier
  # Load required libraries
  required_packages <- c("dplyr", "stringr", "ComplexHeatmap", "grid", "tidyr", "purrr", "dendsort")
  missing_packages <- setdiff(required_packages, rownames(installed.packages()))
  if (length(missing_packages) > 0) {
    stop(paste("Missing required packages:", paste(missing_packages, collapse = ", ")))
  }

  library(dplyr)
  library(stringr)
  library(ComplexHeatmap)
  library(grid)
  library(tidyr)
  library(purrr)
  library(dendsort)

  # Error checking for inputs
  if (!("P.Value" %in% names(df))) {
    stop("Input 'df' must contain a 'P.Value' column.")
  }
  if (!("site_id" %in% names(df))) {
    stop("Input 'df' must contain a 'site_id' column.")
  }
  if (!is.character(comp_name) || length(comp_name) != 1) {
    stop("Parameter 'comp_name' must be a single character string.")
  }
  if (!is.numeric(n) || length(n) != 1 || n <= 0) {
    stop("Parameter 'n' must be a single positive integer.")
  }
  if (!("site_id" %in% names(datal)) || !("name" %in% names(datal)) || !("zscore" %in% names(datal))) {
    stop("Data frame 'datal' must contain 'site_id', 'name', and 'zscore' columns.")
  }
  if (!("name" %in% names(meta)) || !("group" %in% names(meta))) {
    stop("Data frame 'meta' must contain 'name' and 'group' columns.")
  }
  if (!is.list(util_tools) || !("create_named_color_list" %in% names(util_tools))) {
    stop("Parameter 'util_tools' must be a list containing the 'create_named_color_list' function.")
  }

  # Step 1: Select top n sites based on P.Value
  .sel <- df %>%
    arrange(P.Value) %>%
    head(n = as.integer(100 * n))  # Assuming 'n' needs to be scaled by 100 for some reason

  # Step 2: Filter 'datal' for selected 'site_id' and exclude 'Refmix' entries
  .data <- datal %>%
    filter(site_id %in% .sel$site_id) %>%
    filter(!str_detect(name, "Refmix"))

  # Step 3: Pivot to wide format and join with specific columns
  .mat <- .data %>%
    pivot_wider(id_cols = site_id, names_from = name, values_from = zscore) %>%
    left_join(
      datal %>%
        distinct(site_id, ms_lit, lt_lit, sitename),
      by = "site_id"
    ) %>%
    mutate(
      ms_lit = replace_na(ms_lit, 0),
      lt_lit = replace_na(lt_lit, 0)
    )

  # Step 4: Condense groups and select top 'n' entries
  .mat <- condense_groups(.mat) %>%
    head(n = n)

  # Step 5: Align metadata with the matrix columns
  .cids <- intersect(colnames(.mat), meta$name)
  .meta <- meta %>%
    filter(name %in% .cids)

  # Check if all column names in .mat are present in .meta
  if (!all(.cids %in% .meta$name)) {
    warning("Some column names in the matrix do not have corresponding entries in metadata.")
  }

  # Step 6: Create color list for annotations
  .colors <- util_tools$create_named_color_list(.meta, c("group"))

  # Step 7: Create column annotations
  ca <- columnAnnotation(
    group = .meta$group,
    annotation_name_side = "left",
    col = .colors
  )

  # Function to hide zero values in annotations
  hide_zero <- function(x) {
    sapply(x, function(y) ifelse(y == 0, "", as.character(y)))
  }

  # Determine maximum number of characters for annotation width
  nchar_max <- max(
    nchar(.mat$lt_lit),
    nchar(.mat$ms_lit),
    na.rm = TRUE
  )

  .annot_font_size <- 7
  anno_width <- unit((nchar_max + 1) * (.annot_font_size / 12), "picas")

  # Create row annotations
  row_annotations <- HeatmapAnnotation(
    lt_lit = anno_text(
      hide_zero(.mat$lt_lit) %>% str_wrap(width = 10, whitespace_only = FALSE),
      show_name = TRUE,
      gp = gpar(fontsize = .annot_font_size),
      location = unit(0.5, "npc"),
      just = "centre"
    ),
    ms_lit = anno_text(
      hide_zero(.mat$ms_lit) %>% str_wrap(width = 10, whitespace_only = FALSE),
      show_name = TRUE,
      gp = gpar(fontsize = .annot_font_size),
      location = unit(0.5, "npc"),
      just = "center"
    ),
    which = "row",
    gp = gpar(fontsize = 2, size = 2),
    width = anno_width,
    annotation_name_side = "top",
    border = TRUE,
    show_annotation_name = TRUE
  )

  # Prepare the heatmap matrix
  .ht_mat <- .mat %>%
    select(all_of(.meta$name)) %>%
    as.matrix()

  dist_no_na <- util_tools$dist_no_na
  # Generate the heatmap
  ht <- Heatmap(
    .ht_mat,
    width = unit(ncol(.ht_mat) * 0.2, "in"),
    height = unit(nrow(.ht_mat) * 0.24, "in"),
    column_split = .meta$group,
    row_labels = .mat$sitename,
    row_names_gp = gpar(fontsize = 6),
    column_names_side = "top",
    heatmap_legend_param = list(title = "zscore"),
    column_title = comp_name,
    right_annotation = row_annotations,
    top_annotation = ca,
    cluster_rows = dendsort(hclust(dist_no_na(.ht_mat), method = "ward.D2")),
    clustering_method_columns = "ward.D2",


  )

  # Define output file paths
  .outdir <- file.path(OUTDIR, "results", "modi")

  # Create the output directory if it doesn't exist
  if (!dir.exists(.outdir)) {
    dir.create(.outdir, recursive = TRUE)
  }

  .outf <- str_replace_all(comp_name, " - ", "_minus_")
  .fulloutf <- file.path(.outdir, paste0(make.names(.outf), "_top_", nrow(.ht_mat), "x", ncol(.ht_mat), ".pdf"))

  # Draw and save the heatmap
  draw(ht)
  .height <- 7 + (nrow(.ht_mat) * 0.2)
  .width <- 9 + (ncol(.ht_mat) * 0.2)

  pdf(.fulloutf, height = .height, width = .width)
  print(ht)
  dev.off()

  # Generate a sub-selected heatmap for two comparisons
  .comp <- str_extract_all(comp_name, "(?<=group)[A-Za-z0-9]+") %>%
    unlist()

  .cids <- datal %>%
    filter(group %in% .comp) %>%
    pull(name) %>%
    unique()

  .meta2 <- meta %>%
    filter(name %in% .cids)

  .to_remove <- intersect(setdiff(colnames(.mat), .meta2$name), .meta$name)
  .ids <- setdiff(colnames(.mat), .to_remove)

  .mat2 <- .mat %>%
    select(all_of(.ids)) %>%
    as.data.frame()

  # Ensure 'sitename' consistency
  if (!all(.mat2$sitename == .mat$sitename)) {
    warning("Mismatch in 'sitename' between .mat and .mat2.")
  }

  .ht_mat2 <- .mat2 %>%
    select(all_of(.meta2$name)) %>%
    as.matrix()

  # Recreate annotations for the sub-selected heatmap
  ca <- columnAnnotation(
    group = .meta2$group,
    annotation_name_side = "left",
    col = .colors
  )

  row_annotations <- HeatmapAnnotation(
    lt_lit = anno_text(
      hide_zero(.mat2$lt_lit) %>% str_wrap(width = 10, whitespace_only = FALSE),
      show_name = TRUE,
      gp = gpar(fontsize = .annot_font_size),
      location = unit(0.5, "npc"),
      just = "centre"
    ),
    ms_lit = anno_text(
      hide_zero(.mat2$ms_lit) %>% str_wrap(width = 10, whitespace_only = FALSE),
      show_name = TRUE,
      gp = gpar(fontsize = .annot_font_size),
      location = unit(0.5, "npc"),
      just = "center"
    ),
    which = "row",
    gp = gpar(fontsize = 2, size = 2),
    width = anno_width,
    annotation_name_side = "top",
    border = TRUE,
    show_annotation_name = TRUE
  )

  # Generate the sub-selected heatmap
  ht2 <- Heatmap(
    .ht_mat2,
    column_split = .meta2$group,
    width = unit(ncol(.ht_mat2) * 0.2, "in"),
    height = unit(nrow(.ht_mat2) * 0.24, "in"),
    row_labels = .mat2$sitename,
    row_names_gp = gpar(fontsize = 6),
    column_names_side = "top",
    heatmap_legend_param = list(title = "zscore"),
    column_title = comp_name,
    top_annotation = ca,
    right_annotation = row_annotations,
    cluster_rows = dendsort(hclust(dist(.ht_mat2), method = "ward.D2")),
    clustering_method_columns = "ward.D2"
  )

  # Define output file paths for the sub-selected heatmap
  .fulloutf2 <- file.path(
    .outdir,
    paste0(make.names(.outf), "_top_subsel_", nrow(.ht_mat2), "x", ncol(.ht_mat2), ".pdf")
  )

  # Draw and save the sub-selected heatmap
  draw(ht2)
  .height2 <- 7 + (nrow(.ht_mat2) * 0.2)
  .width2 <- 8 + (ncol(.ht_mat2) * 0.2)

  pdf(.fulloutf2, height = .height2, width = .width2)
  print(ht2)
  dev.off()
}

process_cut_by <- function(cut_by, cdesc) {
  print("***")
  print(cut_by)
  # Return NULL immediately if cut_by is NULL
  if (is.null(cut_by)) {
    return(NULL)
  }

  # If cut_by is a single string containing ':', split it into a vector
  if (is.character(cut_by) && length(cut_by) == 1 && grepl(":", cut_by)) {
    cut_by <- strsplit(cut_by, ":")[[1]]
  }

  # Ensure cut_by is now a character vector
  if (!is.character(cut_by)) {
    #warning("cut_by should be a character string or vector.")
    #return(NULL)
    # this is fine
    cut_by <- as.character(cut_by)
  }

  # Check if all elements in cut_by are valid column names
  invalid_cols <- setdiff(cut_by, colnames(cdesc))
  if (length(invalid_cols) > 0) {
    warning(
      "The following cut_by elements are not column names in cdesc: ",
      paste(invalid_cols, collapse = ", ")
    )
    return(NULL)
  }

  # Subset the relevant columns and create the interaction factor
  cut_by_factor <- interaction(cdesc[, cut_by, drop = FALSE], drop = TRUE)

  print("***")
  print(cut_by_factor)

  return(cut_by_factor)
}


# this is an all encompassing heatmap function that can be used to generate heatmaps from gct objects with lit count cols
make_heatmap_fromgct <- function(
    gct,
    row_title = "",
    column_title = "",
    save_func = NULL,
    cluster_rows = T,
    cluster_columns = F,
    color_mapper = NULL,
    sample_order = NULL,
    show_row_names = FALSE,
    cut_by = NULL,
    scale = T,
    optimal_size = FALSE,
    include_row_annotations = TRUE,
    meta_to_include = NULL, # if NULL uses all
    meta_to_exclude = NULL, # if NULL filters nothing out
    ...) {

  # gct <- subgct
  # gct@cdesc$treat <-
  #   factor(.gct@cdesc$treat , levels = c("untreated", "carboplatin", "IMT", "carboplatin_IMT"), ordered = T)

  default_meta_to_exclude <- c("recno", "runno", "searchno", "label", "expr_col", "expr_file")
  # If user specifies additional meta_to_exclude, merge with defaults
  if (!is.null(meta_to_exclude)) {
    meta_to_exclude <- union(default_meta_to_exclude, meta_to_exclude)
  } else {
    meta_to_exclude <- default_meta_to_exclude
  }



  if (!is.null(sample_order)) { # maybe put this order checking later
    gct <- cmapR::subset_gct(gct, cid = sample_order)
  }

  for (col in colnames(gct@cdesc)) {
    .col <- gct@cdesc[[col]]
    .col <- replace(.col, is.na(.col), "NA")
    if (!is.factor(.col)) {
      .col <- util_tools$infer_ordered_factor(.col)
    }
    gct@cdesc[[col]] <- .col
  }


  rdesc <- gct %>% cmapR::melt_gct(keep_cdesc = FALSE, suffixes = c('.x', '.y')) %>% distinct(id.x, .keep_all = TRUE)
  rdesc <- rdesc %>%
    dplyr::select(id.x, ms_lit, lt_lit, sitename) %>% #distinct() %>%
    mutate(
      ms_lit = replace_na(ms_lit, 0),
      lt_lit = replace_na(lt_lit, 0)
    )
  # note the default call of data.frame has check.names = TRUE
  mat_data <- as.data.frame(gct@mat)[rdesc$id.x, ] %>% rownames_to_column(var = "id.x") %>% left_join(rdesc, by = "id.x")
  mat_data <- condense_groups(mat_data) %>% rename(sitename = aggregated_sitename, ms_lit = aggregated_ms_lit, lt_lit = aggregated_lt_lit)
  row_labels <- mat_data$sitename

  heatmap_legend_param <- list(title = "zscore")
  heatmap_matrix_width  <- NULL
  heatmap_matrix_height <- NULL
  if (optimal_size == TRUE){
    heatmap_matrix_width <- unit(ncol(mat_data) * .2, "in")
    heatmap_matrix_height <- unit(nrow(mat_data) * .2, "in")
  }

  cut_by <- process_cut_by(cut_by, gct@cdesc)
  print("***")
  print("***")
  print("***")
  print("***")
  print("***")
  print(cut_by)

  # if (!is.null(cut_by) && cut_by %in% colnames(gct@cdesc)) {
  #   cut_by <- gct@cdesc[[cut_by]]
  #   cut_by <- factor(cut_by, levels = unique(cut_by))
  # } else {
  #   cut_by <- NULL
  # }


  # Determine maximum number of characters for annotation width
  nchar_max <- max(
    nchar(mat_data$lt_lit),
    nchar(mat_data$ms_lit),
    na.rm = TRUE
  )

  .annot_font_size <- 7
  anno_width <- unit((nchar_max + 1) * (.annot_font_size / 12), "picas")


  row_annotations <- NULL
  if (include_row_annotations) {
    row_annotations <- HeatmapAnnotation(
      lt_lit = anno_text(
        hide_zero(mat_data$lt_lit) %>% str_wrap(width = 10, whitespace_only = FALSE),
        show_name = TRUE,
        gp = gpar(fontsize = .annot_font_size),
        location = unit(0.5, "npc"),
        just = "centre"
      ),
      ms_lit = anno_text(
        hide_zero(mat_data$ms_lit) %>% str_wrap(width = 10, whitespace_only = FALSE),
        show_name = TRUE,
        gp = gpar(fontsize = .annot_font_size),
        location = unit(0.5, "npc"),
        just = "center"
      ),
      which = "row",
      gp = gpar(fontsize = 2.2, size = 2.2),
      width = anno_width,
      annotation_name_side = "bottom",
      border = TRUE,
      show_annotation_name = TRUE,
      annotation_name_gp = gpar(fontsize = 8, size = 8)
    )
  }


  column_annotations <- NULL

  cmeta <- gct@cdesc
  if (!is.null(meta_to_include)){
    message('subsetting cmeta by inclusion list')
    cmeta <- cmeta %>% dplyr::select(any_of(meta_to_include))
  }
  if (!is.null(meta_to_exclude)){
    message('subsetting cmeta by exclusion list')
    cmeta <- cmeta %>% dplyr::select(-any_of(meta_to_exclude))

  }
  #.colors <- util_tools$create_named_color_list(gct@cdesc, colnames(gct@cdesc))
  message(meta_to_include)
  message(colnames(cmeta))
  .colors <- util_tools$create_named_color_list(cmeta, colnames(cmeta))
  # browser()
  column_annotations <- ComplexHeatmap::columnAnnotation(
      df = cmeta,
      col = .colors,
      annotation_name_side = "right"
    )

  .legend_width <- unit(6, "cm")
  .cbar_title <- "zscore"
  heatmap_legend_param <- list(
    title = .cbar_title,
    direction = "horizontal",
    just = "bottom",
    legend_width = .legend_width
  )

  dist_no_na <- util_tools$dist_no_na

  ht <- ComplexHeatmap::Heatmap(
    mat_data[ , colnames(gct@mat)],
    width = heatmap_matrix_width,
    height = heatmap_matrix_height,
    # TODO use z_score_withna or other custom func for handling nas when scaling
    row_labels = row_labels,
    show_row_names = show_row_names,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    #

    row_title = row_title,
    # column_title = column_title,
    column_title_gp = gpar(fontsize = 9, fontface = "bold", just = "left"),

    column_labels = gct@cdesc$id, # id should always be here

    column_split = cut_by,


    top_annotation = column_annotations,
    right_annotation = row_annotations,

    heatmap_legend_param = heatmap_legend_param,
    row_names_gp = grid::gpar(fontsize = 7.4),
    column_names_gp = grid::gpar(fontsize = 7.4),
    cluster_column_slices = FALSE,
    column_names_side = "top",

    clustering_distance_rows = dist_no_na,
    clustering_distance_columns = dist_no_na
  )



  # log_msg(debug = paste0("defining draw func"))
  do_draw <- function() {
    draw(ht,
      heatmap_legend_side = "bottom",
      column_title = column_title,
      column_title_gp = gpar(fontsize = 13, fontface = "bold", just = "left"),
      padding = unit(c(2, 24, 2, 24), "mm"), # top, left, bottom, right
    )
  }

  # log_msg(debug = paste0("save func: ", class(save_func) %>% as.character()))
  # log_msg(debug = paste0("is null save func: ", is.null(save_func)))

  height <- 4 + (nrow(mat_data) * .20)
  width <- 8 + (ncol(mat_data) * .26)

  if (!is.null(save_func)) {

    # log_msg(debug = paste0("save func attrs before: "))
    # log_msg(debug = paste0(names(get_args(save_func)), "-", get_args(save_func)))

    save_func <- make_partial(save_func, height = height, width = width)

    # log_msg(debug = paste0("save func attrs after: "))
    # log_msg(debug = paste0(names(get_args(save_func)), "-", get_args(save_func)))

    save_func(plot_code = do_draw)
  }


  return(do_draw)
}

