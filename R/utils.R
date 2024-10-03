
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
                         "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
                         "#bcbd22", "#17becf")

  # Initialize an empty list to store the result
  color_list <- list()

  # Iterate over each column
  for (col_name in columns) {
    unique_vals <- unique(df[[col_name]])
    n_vals <- length(unique_vals)

    # Assign colors based on the number of unique values, recycling matplotlib colors if needed
    colors_assigned <- rep(matplotlib_colors, length.out = n_vals)

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
  out <- scale(value)

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

