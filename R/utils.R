suppressPackageStartupMessages(library(tidyr))


normalize <- function(datal, id_col = "site_id") {
  # Ensure that 'id_col' exists in 'datal'
  if (!id_col %in% names(datal)) {
    stop(paste("Column", id_col, "not found in 'datal'"))
  }

  # Compute the minimum positive 'value'
  minval <- datal %>%
    dplyr::filter(value > 0) %>%
    dplyr::pull(value) %>%
    min()

  # Normalize 'value' and compute 'log2_val' and 'zscore'
  datal <- datal %>%
    dplyr::filter(value > 0) %>%
    dplyr::group_by(name, recno, runno, searchno, label, taxon) %>%
    dplyr::mutate(
      val_mednorm = (value + (minval / 2)) / median(value, na.rm = TRUE)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      log2_val = log2(val_mednorm)
    ) %>%
    dplyr::group_by(!!rlang::sym(id_col)) %>%
    dplyr::mutate(
      zscore = myzscore(log2_val)
    ) %>%
    dplyr::ungroup()

  # Count non-NA 'log2_val' per 'id_col'
  counts <- datal %>%
    dplyr::filter(!is.na(log2_val)) %>%
    dplyr::group_by(!!rlang::sym(id_col)) %>%
    dplyr::summarize(n = n()) %>%
    dplyr::arrange(desc(n))

  # Determine 'maxn' based on 'counts' and 'sample names'
  maxn <- min(max(counts$n), length(unique(datal[['name']])))

  # Identify groups to keep
  tokeep <- counts %>%
    dplyr::filter(n == maxn) %>%
    dplyr::pull(!!rlang::sym(id_col))

  # Filter 'datal' to include only groups in 'tokeep'
  datal <- datal %>%
    dplyr::filter(!!rlang::sym(id_col) %in% tokeep)

  return(datal)
}



create_named_color_list <- function(df, columns) {
  # Matplotlib default colors (from memory)
  matplotlib_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
                         "#9467bd", "#8c564b", "#e377c2",
                         "#bcbd22", "#17becf")

  # Initialize an empty list to store the result
  color_list <- list()

  # Iterate over each column
  for (col_name in columns) {
    df[[col_name]] <-  df[[col_name]] %>% replace_na("NA")
    unique_vals <- unique(df[[col_name]])
    # message(paste0('unique values are ', unique_vals))
    n_vals <- length(unique_vals)

    # Assign colors based on the number of unique values, recycling matplotlib colors if needed
    colors_assigned <- rep(matplotlib_colors, length.out = n_vals)
    color_map <- setNames(colors_assigned, unique_vals)


    if ("NA" %in% unique_vals) {
      color_map["NA"] <- "#646464"
    }
    if ("NULL" %in% unique_vals) {
      color_map["NULL"] <- "#646464"
    }


    # Create a named vector for the unique values with corresponding colors
    color_list[[col_name]] <- setNames(colors_assigned, unique_vals)
  }

  return(color_list)
}


myzscore <- function(value, minval = NA, remask = TRUE) {
  mask <- is.na(value)
  if (is.na(minval)) minval <- min(value, na.rm = TRUE)

  if (minval == Inf) {
    minval <- 0
  }

  value[is.na(value)] <- minval
  # todo make smaller than min val
  out <- scale(value, center = TRUE, scale = TRUE) # these are defaults, just being explicit

  # if all NA:
  if (sum(!is.finite(out)) == length(out)) {
    out[, 1] <- 0
  }

  if (remask == TRUE) {
    out[mask] <- NA
  }
  return(out)
}

dist_no_na <- function(mat) {
  mat[is.na(mat)] <- min(mat, na.rm = TRUE)
  edist <- dist(mat)
  return(edist)
}

normalize_gct <- function(gct, log_transform=FALSE) {
  # log_msg(msg="zscoring gct file by row")
  # log_msg(msg=paste0("group_by is set to: ", group_by))
  res <- gct %>%
    melt_gct() %>%
    group_by(id.y) %>%
    dplyr::mutate(value_norm = value/median(value, na.rm=T)) %>%
    dplyr::ungroup() 

    if (log_transform) {
        res %<>% mutate(value_norm = log2(value_norm))
    } 

  # make a new gct and return
  res <- res %>%
    dplyr::select(id.x, id.y, value_norm) %>%
    tidyr::pivot_wider(names_from = id.y, values_from = value_norm) %>%
    as.data.frame()
  rownames(res) <- res$id.x
  res$id.x <- NULL
  res <- as.matrix(res)
  .minval <- res[ res > -Inf] %>% min(na.rm=T)
  res[res == -Inf] = NA
  res <- res + abs(.minval*1.1)
  rdesc <- gct@rdesc
  cdesc <- gct@cdesc

  if (length(colnames(gct@rdesc)) == 1) {
    rdesc["dummy"] <- "X"
  }

  if (length(colnames(res)) == 1) {
    rdesc["dummy"] <- "X"
  }

  newgct <- new("GCT",
    mat = res,
    rid = rownames(res),
    cid = colnames(res),
    rdesc = rdesc[rownames(res), ],
    cdesc = cdesc[colnames(res), ]
  )

  return(newgct)
}

scale_gct <- function(gct, group_by = NULL) {
  # log_msg(msg="zscoring gct file by row")
  # log_msg(msg=paste0("group_by is set to: ", group_by))
  if (!is.null(group_by) && group_by == FALSE) group_by <- NULL
  res <- gct %>%
    melt_gct() %>%
    {
      # Conditionally add group_by
      if (!is.null(group_by) && group_by != FALSE) {
        group_by(., id.x, !!sym(group_by))
      } else {
        group_by(., id.x)
      }
    } %>%
    dplyr::mutate(zscore = myzscore(value)) %>%
    dplyr::ungroup()

  # make a new gct and return
  res <- res %>%
    dplyr::select(id.x, id.y, zscore) %>%
    tidyr::pivot_wider(id_cols=id.x, names_from = id.y, values_from = zscore) %>%
    as.data.frame()
  rownames(res) <- res$id.x
  res$id.x <- NULL
  res <- as.matrix(res)
  rdesc <- gct@rdesc
  cdesc <- gct@cdesc

  if (length(colnames(gct@rdesc)) == 1) {
    rdesc["dummy"] <- "X"
  }

  if (length(colnames(res)) == 1) {
    rdesc["dummy"] <- "X"
  }

  newgct <- new("GCT",
    mat = res,
    rid = rownames(res),
    cid = colnames(res),
    rdesc = rdesc[rownames(res), ],
    cdesc = cdesc[colnames(res), ]
  )

  return(newgct)
}

infer_ordered_factor <- function(column) {
  # Function to infer ordering for a given column vector

  # Extract numeric values from strings
  numeric_values <- as.numeric(gsub("[^0-9\\.]+", "", column))

  # gsub("[^0-9]+", "", column)
  # Find the minimum numeric value
  min_num <- min(numeric_values, na.rm = TRUE)
  if (is.infinite(min_num)) {
    min_num <- 0  # Default to 0 if no numeric values are found
  }

  # Initialize order values with numeric values
  order_values <- numeric_values

  # Indices of non-numeric values
  non_numeric_indices <- which(is.na(order_values))

  if (length(non_numeric_indices) > 0) {
    non_numeric_values <- tolower(column[non_numeric_indices])

    # Assign special order value for "veh" or "vehicle"
    is_vehicle <- grepl("^veh$|^vehicle$", non_numeric_values)
    order_values[non_numeric_indices[is_vehicle]] <- min_num - 2  # Highest priority

    is_ctrl <- grepl("^veh$|^control$", non_numeric_values)
    order_values[non_numeric_indices[is_ctrl]] <- min_num - 2  # Highest priority

    # Assign next priority to other non-numeric values
    order_values[non_numeric_indices[!is_vehicle]] <- min_num - 1
  }

  # browser()

  # Create a data frame for sorting
  df <- data.frame(
    original_value = column,
    order_value = order_values,
    stringsAsFactors = FALSE
  )

  # Remove duplicates while preserving order
  df_unique <- df[!duplicated(df$original_value), ]

  # Sort the data frame by order_value
  df_ordered <- df_unique[order(df_unique$order_value), ]

  # Create an ordered factor with levels in the sorted order
  ordered_factor <- factor(
    column,
    levels = df_ordered$original_value,
    ordered = TRUE
  )

  return(ordered_factor)
}



filter_observations <- function(df, column, threshold_ratio) {
  # Validate inputs
  if (!is.data.frame(df)) {
    stop("Input 'df' must be a data frame.")
  }
  if (!(column %in% colnames(df))) {
    stop(paste("Column", column, "not found in the data frame."))
  }
  if (!is.numeric(threshold_ratio) || threshold_ratio <= 0 || threshold_ratio > 1) {
    stop("threshold_ratio should be a numeric value between 0 and 1.")
  }
  if (!"site_id" %in% colnames(df)) {
    stop("Column 'site_id' not found in the data frame.")
  }
  if (!"log2_val" %in% colnames(df)) {
    stop("Column 'log2_val' not found in the data frame.")
  }

  # Step 1: Group by site_id and the specified column, then count observations
  counts <- df %>%
    filter(!is.na(log2_val)) %>%
    group_by(site_id, !!sym(column)) %>%
    summarize(n = n(), .groups = 'drop')

  # Step 2: Compute total counts per column factor
  # total_counts <- df %>%
  #   group_by(id.y, !!sym(column)) %>%
  #   summarize(total_n = sum(n), .groups = 'drop')
  total_counts <-  df %>% distinct(name, .keep_all = T) %>% group_by(!!sym(column)) %>% summarize(total_n = n())

  # Step 3: Join counts with total_counts and calculate ratio per group
  counts <- counts %>%
    left_join(total_counts, by = column) %>%
    mutate(ratio = n / total_n)

  # Step 4: Inspect counts and ratios
  message("Counts and Ratios per Group:")
  print(counts)

  # Step 5: Filter groups based on the threshold_ratio
  tokeep <- counts %>%
    filter(ratio >= threshold_ratio)

  message(paste("Number of groups meeting the threshold_ratio of", threshold_ratio, ":", nrow(tokeep)))

  # Step 6: Inspect 'tokeep' DataFrame
  if(nrow(tokeep) > 0){
    print(tokeep)
  } else {
    message("No groups meet the threshold ratio. Consider lowering the threshold_ratio.")
  }

  # Step 7: Filter the original dataframe to keep only the desired site_id and column combinations
  df_filtered <- df %>%
    inner_join(tokeep, by = c("site_id", column))

  message(paste("Number of rows after filtering:", nrow(df_filtered)))

  return(df_filtered)
}



condense_groups <- function(.mat, include_annotations = TRUE) {
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

  # Define the annotation columns if needed
  required_non_numeric_exclusions <- if (include_annotations) c("ms_lit", "lt_lit") else character(0)

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
    condensed_data %<>% rename(sitename = aggregated_sitename)
    if ('ms_lit' %in% colnames(condensed_data)){
      condensed_data %<>% rename(ms_lit = aggregated_ms_lit, lt_lit = aggregated_lt_lit)
    }
  }

  return(condensed_data)
}



# Define the select_top_n function
select_top_n <- function(data,
                        sort_column = "t_value",
                        top_n = 10,
                        sort_order = "desc",
                        distinct = TRUE,
                        allow_duplicates = FALSE) {
  # Arguments:
  # data: Data frame to operate on
  # sort_column: Column name to sort/filter on
  # top_n: Number of top entries to select
  # sort_order: "asc" for ascending, "desc" for descending
  # distinct: If TRUE, select distinct values
  # allow_duplicates: If TRUE, allow multiple rows with the same sort_column value

  # Check if sort_column exists
  if (!(sort_column %in% colnames(data))) {
    stop(paste("Column", sort_column, "not found in the data frame."))
  }

  # Check sort_order validity
  if (!(sort_order %in% c("asc", "desc"))) {
    stop("sort_order must be either 'asc' or 'desc'.")
  }

  # Arrange data based on sort_order
  if (sort_order == "desc") {
    arranged_data <- data %>% arrange(desc(.data[[sort_column]]))
  } else {
    arranged_data <- data %>% arrange(.data[[sort_column]])
  }

  # Select distinct values if required
  if (distinct) {
    arranged_data <- arranged_data %>% distinct(.data[[sort_column]], .keep_all = TRUE)
  }

  # Select top_n values
  # If duplicates are not allowed, ensure we keep n distinct values
  if (!allow_duplicates) {
    top_values <- arranged_data %>% distinct(.data[[sort_column]], .keep_all = TRUE) %>% head(n = top_n)
  } else{
    top_values <- arranged_data %>% head(n = top_n) %>% pull(.data[[sort_column]])
    }
  filtered_data <- arranged_data %>% filter(.data[[sort_column]] %in% top_values[[sort_column]])

  return(filtered_data)
}
