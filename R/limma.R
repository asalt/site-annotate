suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(here))

source(file.path(here("R", "lazyloader.R")))
util_tools <- get_tool_env("utils")


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
    contrast_strings <- lapply(contrast_combinations, function(combo) {
      paste(paste0(group_var, combo[1]), "-", paste0(group_var, combo[2]))
    }) %>% unlist()
    return(contrast_strings)
  }

  contrast_strings <- config$limma$contrasts %||% generate_contrasts(gct@cdesc, config$limma$group_var %||% "group")
  contrasts <- makeContrasts(contrasts = contrast_strings, levels = mod)

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
    .table %<>% rownames_to_column(var = "id")
    .table$contrast <- .contrast

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
    return(.table2)
  })

  return(toptables2)
}

plot_p_dist <- function(df, main_title = ""){

  par(mfrow = c(1, 2))  # 1 row, 2 columns

  # Plot raw p-values
  p_values <- df[["P.Value"]]
  p_adjusted <- df[["adj.P.Val"]]
  hist(
    p_values,
    breaks = 30,
    main = "Raw P-values",
    xlab = "P-value",
    ylab = "Count",
    col = "lightblue",
    border = "black"
  )

  mtext(main_title, side = 1, line = -1, outer=T, adj = 0, cex = 0.7)
  # put at bottom left

  # Plot adjusted p-values
  hist(
    p_adjusted,
    breaks = 30,
    main = "Adjusted P-values",
    xlab = "Adjusted P-value",
    ylab = "Count",
    col = "lightcoral",
    border = "black"
  )
  #mtext(main_title, outer = TRUE, cex = 1.5, line = -1)

  # Reset plot layout
  par(mfrow = c(1, 1))


}


#
