suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(here))

sys.source(file.path(here::here("R", "lazyloader.R")), envir = environment())
util_tools <- get_tool_env("utils")

# Make common test variables visible in this environment when sourced via sys.source
if (!exists("gct", envir = environment()) && exists("gct", envir = .GlobalEnv)) {
  assign("gct", get("gct", envir = .GlobalEnv), envir = environment())
}
if (!exists("config", envir = environment()) && exists("config", envir = .GlobalEnv)) {
  assign("config", get("config", envir = .GlobalEnv), envir = environment())
}


impute_with_draw <- function(values, seed = 123, mean_adjust = -1.8, sd_adjust = 0.8) {
  # Set seed for reproducibility
  set.seed(seed)

  # Ensure the input is numeric and has at least one non-NA value
  if (!is.numeric(values) || all(is.na(values))) {
    stop("Input must be a numeric vector with at least one non-NA value.")
  }

  # Return early if there are no missing values
  if (all(!is.na(values))) {
    return(values)
  }

  # Calculate mean and standard deviation of non-NA values
  original_mean <- mean(values, na.rm = TRUE)
  original_sd <- sd(values, na.rm = TRUE)

  # Adjust mean and sd
  adjusted_mean <- original_mean + mean_adjust
  adjusted_sd <- original_sd * sd_adjust

  # Number of missing values
  n <- sum(is.na(values))

  # Generate 8 random draws and average them
  random_draw_final <- rowMeans(
    matrix(rnorm(n * 8, mean = adjusted_mean, sd = adjusted_sd), nrow = 8)
  )

  # Replace NA values with the random draws
  values[is.na(values)] <- random_draw_final

  # Return the imputed vector
  return(values)
}

impute_with_draw <- util_tools$impute_with_draw


parse_named_contrasts <- function(strings) {
  if (is.null(strings) || length(strings) == 0) {
    return(NULL)
  }

  parsed <- lapply(strings, function(entry) {
    entry <- as.character(entry)
    parts <- strsplit(entry, "=", fixed = TRUE)[[1]]

    if (length(parts) > 1) {
      label <- trimws(parts[1])
      expr <- trimws(paste(parts[-1], collapse = "="))
    } else {
      label <- trimws(entry)
      expr <- trimws(entry)
    }

    if (!nzchar(label)) {
      label <- expr
    }

    list(label = label, expression = expr)
  })

  labels <- vapply(parsed, `[[`, character(1), "label")
  expressions <- vapply(parsed, `[[`, character(1), "expression")

  data.frame(
    label = make.unique(labels),
    expression = expressions,
    stringsAsFactors = FALSE
  )
}


run_limma <- function(gct, config, do_impute = TRUE) {
  .formula <- config$limma$formula %||% as.formula("~0 + group")
  print(paste0("formula is : ", .formula))
  mod <- model.matrix(.formula, gct@cdesc)

  # rownames(mod) <- meta$name
  .mat <- mat(gct)
  colnames(mod) <- make.names(colnames(mod))
  print(mod)
  # or impute here
  if (do_impute) .mat <-  impute_with_draw(.mat)
  .mat[is.na(.mat)] <- 0

  generate_contrasts <- function(meta, group_var) {
    # Get unique levels of the grouping variable
    if (!group_var %in% colnames(meta)) {
      warning(paste0("Grouping variable '", group_var, "' not found in metadata"))
      return()
    }
    group_levels <- unique(meta[[group_var]])

    # Generate all pairwise combinations of group levels
    contrast_combinations <- combn(group_levels, 2, simplify = FALSE)

    # Create contrasts for each combination
    # Heuristic ordering: put the longer label first (falls back to alphabetical)
    contrast_strings <- lapply(contrast_combinations, function(combo) {
      a <- as.character(combo[1]); b <- as.character(combo[2])
      if (nchar(a) < nchar(b) || (nchar(a) == nchar(b) && a < b)) {
        paste(paste0(group_var, b), "-", paste0(group_var, a))
      } else {
        paste(paste0(group_var, a), "-", paste0(group_var, b))
      }
    }) %>% unlist()
    return(contrast_strings)
  }

  # Determine grouping variable for automatic contrasts
  group_var <- config$limma$group_var
  if (is.null(group_var)) {
    # Try to infer from the formula (drop constants like 0/1)
    fvars <- setdiff(all.vars(.formula), c("0", "1"))
    if (length(fvars) >= 1) {
      group_var <- fvars[[length(fvars)]]
    } else {
      group_var <- "group"
    }
  }

  contrast_input <- config$limma$contrasts
  if (is.null(contrast_input) || length(contrast_input) == 0) {
    contrast_input <- generate_contrasts(gct@cdesc, group_var)
  }

  contrast_defs <- parse_named_contrasts(contrast_input)
  if (is.null(contrast_defs) || nrow(contrast_defs) == 0) {
    stop("No LIMMA contrasts could be derived from configuration or metadata.")
  }

  contrast_exprs <- contrast_defs$expression
  contrasts <- makeContrasts(contrasts = contrast_exprs, levels = mod)
  colnames(contrasts) <- contrast_defs$label
  contrast_lookup <- structure(contrast_defs$expression, names = contrast_defs$label)

  fit <- limma::lmFit(.mat, mod)
  fit <- contrasts.fit(fit, contrasts)
  fitres <- fit %>% limma::eBayes(robust = T, trend = T)

  toptables <- colnames(contrasts) %>%
    purrr::map(~ {
      limma::topTable(fitres, coef = .x, number = Inf, sort.by = "none")
    }) %>%
    setNames(nm = colnames(contrasts))

  .data <- gct %>% melt_gct(keep_rdesc = T, keep_cdesc = F)

  toptables2 <- toptables %>% purrr::imap(~ {
    # assertthat::are_equal(nrow(.x), nrow(.mat))  # this doesn't have to be true
    .table <- .x
    .contrast <- .y
    .expr <- contrast_lookup[[.contrast]] %||% .contrast
    .table %<>% rownames_to_column(var = "id")
    .table$contrast <- .contrast
    .table$contrast_expression <- .expr

    # # recalculate adj p val on unique tests - estimated by unique t values.
    # .recal <- .table %>%
    #   distinct(t, .keep_all = T) %>%
    #   mutate(adj.P.Val = p.adjust(adj.P.Val))
    # browser()
    # .table <- left_join(
    #   .table %>% select(-adj.P.Val),
    #   .recal %>% select(t, adj.P.Val),
    # )

    .table2 <- right_join(gct@rdesc, .table, by = "id")
    .table2$contrast_expression <- .expr
    return(.table2)
  })

  attr(toptables2, "contrast_expression") <- contrast_lookup
  ## Build a text summary for logs/files
  dm_cols <- paste(colnames(mod), collapse = ", ")
  grp_levels <- unique(gct@cdesc[[group_var]]) %||% NA
  grp_levels <- paste(grp_levels, collapse = ", ")
  contrast_lines <- paste(paste(names(contrast_lookup), contrast_lookup, sep = ": "), collapse = "\n")
  summary_text <- paste0(
    "LIMMA summary\n",
    "formula: ", deparse(.formula), "\n",
    "group_var: ", group_var, "\n",
    "design columns: ", dm_cols, "\n",
    "levels(", group_var, "): ", grp_levels, "\n",
    "contrasts:\n",
    contrast_lines, "\n"
  )
  attr(toptables2, "limma_summary") <- summary_text
  return(toptables2)
}

plot_p_dist <- function(df, main_title = ""){

  par(mfrow = c(1, 2))  # 1 row, 2 columns

  # Plot raw p-values
  p_values <- df[["P.Value"]]
  p_adjusted <- df[["adj.P.Val"]]
  raw_hist <- hist(p_values, breaks = 30, plot = FALSE)
  adj_hist <- hist(p_adjusted, breaks = 30, plot = FALSE)
  max_count <- max(c(raw_hist$counts, adj_hist$counts), na.rm = TRUE)

  hist(
    p_values,
    breaks = raw_hist$breaks,
    main = "Raw P-values",
    xlab = "P-value",
    ylab = "Count",
    col = "lightblue",
    border = "black",
    ylim = c(0, max_count)
  )

  mtext(main_title, side = 1, line = -1, outer=T, adj = 0, cex = 0.7)
  # put at bottom left

  # Plot adjusted p-values
  hist(
    p_adjusted,
    breaks = adj_hist$breaks,
    main = "Adjusted P-values",
    xlab = "Adjusted P-value",
    ylab = "Count",
    col = "lightcoral",
    border = "black",
    ylim = c(0, max_count)
  )
  #mtext(main_title, outer = TRUE, cex = 1.5, line = -1)

  # Reset plot layout
  par(mfrow = c(1, 1))


}


#
