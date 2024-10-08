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
    message(paste0('unique values are ', unique_vals))
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
    tidyr::pivot_wider(names_from = id.y, values_from = zscore) %>%
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
  numeric_values <- as.numeric(gsub("[^0-9.-]+", "", column))

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

    # Assign next priority to other non-numeric values
    order_values[non_numeric_indices[!is_vehicle]] <- min_num - 1
  }

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


filter_observations <- function(df, column, threshold_value) {
  # Validate threshold_value
  if (!is.numeric(threshold_value) || threshold_value <= 0) {
    stop("threshold_value should be a positive numeric value")
  }

  # Group by site_id and the specified column, then count observations
  counts <- df %>%
    dplyr::filter(!is.na(log2_val)) %>%
    dplyr::group_by(site_id, .data[[column]]) %>%
    dplyr::summarize(n = n(), .groups = 'drop')

  # Determine if threshold is an absolute count or a fraction
  if (threshold_value < 1) {
    threshold <- threshold_value * max(counts$n)
  } else {
    threshold <- threshold_value
  }

  # Filter groups that meet or exceed the threshold
  tokeep <- counts %>% dplyr::filter(n >= threshold)

  # Filter the original dataframe to keep only the desired site_ids
  df_filtered <- df %>% dplyr::filter(site_id %in% tokeep$site_id)
  return(df_filtered)
}

