library(here)

source(file.path(here(), "R", "lazyloader.R"))

# this only works for sitename level
do_pca <- function(gct,  scale = T, removeVar=0.1) {

  library(PCAtools)

  .mat <- mat(gct)
  #.mat[is.na(.mat)] <- min(.mat, na.rm = T)
  .mat[is.na(.mat)] <- mean(.mat, na.rm = T)
  if ("sitename" %in% colnames(gct@rdesc)){
      rownames(.mat) <- gct@rdesc[rownames(.mat), "sitename"] %>% make.unique()
  }

  .mat <- .mat[!duplicated(.mat), ]  # drop duplicates
  
  pca_obj <- PCAtools::pca(
    .mat,
    metadata = gct@cdesc,
    center = T,
    scale = scale,
    removeVar=removeVar
  )

  return(pca_obj)
}

plot_biplot <- function(
    pca_object,
    top_pc = 3,
    showLoadings = T,
    labSize = 1.8,
    pointSize = 3.2,
    sizeLoadingsNames = 2,
    colby = NULL, # or a string like 'group'
    shape = NULL, # or a string like 'group'
    encircle = !is.null(colby),
    ntopLoadings = 5,
    title = "",
    basename = "pca_biplot_",
    fig_width = 6.2,
    fig_height = 5.8,
    ...) {

  if (!is.logical(encircle) || length(encircle) != 1 || is.na(encircle)) {
    stop("`encircle` must be a single logical value (TRUE or FALSE). Received: ", encircle)
  }

  if (!is.null(shape) && !is.null(pca_object$metadata) && shape %in% colnames(pca_object$metadata)) {
      pca_object$metadata[[shape]] <- as.factor(pca_object$metadata[[shape]])
  }



  vec <- paste0("PC", 1:top_pc)
  # vec <- c("PC1", "PC2", "PC3") # "PC4")
  pcs <- combn(vec, 2) %>%
    as.data.frame() %>%
    as.list()
  #

  # log_msg(info=paste0('colby is: ', colby))
  # log_msg(info=paste0('metadata is : ', pca_object$metadata))
  if (!is.null(colby) &&
    !is.na(colby) &&
    !is.null(pca_object$metadata) &&
    !colby %in% colnames(pca_object$metadata)) {
    warning(paste0(colby, " not found in metadata"))
    colby <- NULL
    encircle <- F
  }

  if (!is.null(shape) &&
    !is.null(pca_object$metadata) &&
    !shape %in% colnames(pca_object$metadata)) {
    warning(paste0(shape, " not found in metadata"))
    shape <- NULL
  }

  n_features <- nrow(pca_object$loadings)
  footnote <- paste0("n = ", as.character(n_features))
  if (nzchar(title)) {
    title_final <- paste0(title, " (n=", n_features, ")")
  } else {
    title_final <- paste0("PCA (n=", n_features, ")")
  }


  plts <- pcs %>%
    purrr::map(
      ~ {
        # stopifnot(~COLBY%in%colnames(.metadata))
        .x1 <- .x[[1]]
        .x2 <- .x[[2]]

        if (!.x1 %in% names(pca_object$rotated)) {
          warning("not enough PCs")
          return()
        }

        if (!.x2 %in% names(pca_object$rotated)) {
          warning("not enough PCs")
          return()
        }

        .max_x <- max(
          pca_object$rotated[[.x1]] %>% abs(),
          pca_object$rotated[[.x2]] %>% abs()
        )
        .max_y <- .max_x

        plt <- PCAtools::biplot(
          pca_object,
          x = .x1,
          y = .x2,
          showLoadings = showLoadings,
          labSize = labSize,
          pointSize = pointSize,
          alphaLoadingsArrow = 0.5,
          sizeLoadingsNames = sizeLoadingsNames,
          colby = colby,
          shape = shape,
          # shape="source",
          legendPosition = "right",
          legendLabSize = 9,
          encircle = encircle,
          title = title_final,
          max.overlaps = Inf,
          maxoverlapsConnectors = Inf,
          ntopLoadings = ntopLoadings %||% 5,
        ) +
          xlim(-.max_x, .max_x) +
          ylim(-.max_y, .max_y) +
          coord_fixed(ratio = 1) +
          labs(caption = footnote) +
          geom_hline(yintercept = 0, color = "grey50", show.legend = NA) +
          geom_vline(xintercept = 0, color = "grey50", show.legend = NA) +
          colorspace::scale_color_discrete_qualitative(palette = "Dynamic") +
          colorspace::scale_fill_discrete_qualitative(palette = "Dynamic") +
          scale_shape_manual(values = c(16, 17, 15, 7, 9, 12, 13, 14)) +  # this fails if more than 8
          theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
          )

        filename <- paste0(basename, "_", as.character(n_features), "_", .x1, "_", .x2, ".pdf")
        ggplot2::ggsave(filename, plt, width = fig_width, height = fig_height)
        return(plt)
      } # exit inner
    ) # exit outer
  return(plts)
}



get_top_loadings <- function(pcaobj,
                             components = c("PC1", "PC2", "PC3", "PC4"),
                             rangeRetain = 0.05,
                             limit = Inf
                             ) { # directly lifted from https://github.com/kevinblighe/PCAtools/blob/50f240ba76fe07e2a546998e76dc8aa1c3570693/R/plotloadings.R#L7 and edited
  x <- pcaobj$loadings[, components, drop = FALSE]
  membership <- list()
  retain <- c()
  # could this be rewritten with tidy group_by and arrange and slice top largest diff?
  for (i in seq_along(components)) {
    # for each PC, based on the loadings range, calculate the rangeRetain
    # fraction of the range
    offset <- (max(x[, i]) - min(x[, i])) * rangeRetain

    # to obtain upper and lower cut-offs, offset max and min by the offset
    uppercutoff <- max(x[, i]) - offset
    lowercutoff <- min(x[, i]) + offset

    # build a consensus list of variables that will be included
    retain_vals <- c(
      which(x[, i] >= uppercutoff),
      which(x[, i] <= lowercutoff)
    ) %>% unique()
    if (limit < Inf) {
      retain_vals <- order(abs(x[retain_vals, i]), decreasing=T) %>% head(n=limit)
    }

    retain <- unique(c(
      retain,
      retain_vals
    ))

    membership[[paste0("PC", i, "_member")]] <- retain_vals
  }
  membership_df <- stack(membership)
  membership_df[["var"]] <- rownames(x)[membership_df$values]
  membership_df[["boolvalue"]] <- !is.na(membership_df$var)
  membership_dfw <- membership_df %>% pivot_wider(id_cols = c("values", "var"), names_from = "ind", values_from = "boolvalue")
  membership_dfw[is.na(membership_dfw)] <- FALSE

  message("-- variables retained:")
  message(paste0(rownames(x)[retain], collapse = ", "))
  x_final <- x[retain, , drop = FALSE]
  # create a plot object (2-col df) of PC names and
  # explained variance of each
  x_final <- data.frame(rownames(x_final), x_final[, components, drop = FALSE])
  colnames(x_final)[1] <- "var"

  x_final <- x_final %>% left_join(membership_dfw)
  x_final$values <- NULL
  x_final %<>% as.data.frame
  rownames(x_final) <- x_final$var
  return(x_final)
}



make_heatmap_from_loadings <- function(
    gct_object,
    pca_object,
    components = c("PC1", "PC2", "PC3", "PC4"),
    save_func = NULL,
    ...){



  top_loadings <- get_top_loadings(pca_object, components = components, rangeRetain = 1.05, limit = Inf)
  # top_loadings returns a dataframe with rownames as features and a column called "var" that is equal to rownames(top_loadings)
  # and floats PC1, 2, 3 .. and boolean PC1_member, PC2_member, PC3_member, ... (member being within the top 5%)

  if ("sitename" %in% colnames(gct_object@rdesc)){  # this needs to not be hardcoded or is not usable for other dtypes,  but is the case for heatmap too
    site_ids <- gct_object@rdesc %>% dplyr::filter(sitename %in% rownames(top_loadings))
  }
  submat <- gct %>% filter(pathway %in% rownames(top_loadings))

  ra <- ComplexHeatmap::rowAnnotation(
    PC1 = top_loadings$PC1_member %>% as.character() %>% anno_simple(col = c("TRUE" = "black", "FALSE" = "white")),
    PC2 = top_loadings$PC2_member %>% as.character() %>% anno_simple(col = c("TRUE" = "black", "FALSE" = "white")),
    PC3 = top_loadings$PC3_member %>% as.character() %>% anno_simple(col = c("TRUE" = "black", "FALSE" = "white")),
    PC4 = top_loadings$PC4_member %>% as.character() %>% anno_simple(col = c("TRUE" = "black", "FALSE" = "white"))
  )

  ht <- plot_tools$plot_results_one_collection(
    df = submat,
    metadata = pca_object$metadata, # guaranteed to match the colnames of df which was used originally to make the pca obj
    # pathway_metadata = loadings,
    row_annotation = ra,
    limit = Inf,
    pstat_cutoff = 1, # we want no filtering
    save_func = save_func,
    title = "Top 5% Loadings",
    ...
  )

  top_loadings <- get_top_loadings(pca_object, components = components, rangeRetain = 1, limit = 5)
  submat <- gct_object %>% filter(pathway %in% rownames(top_loadings))
  # submat %<>% mutate(pathway = factor(pathway, levels = rownames(top_loadings), ordered=TRUE)) %>% arrange(pathway)

  ra <- ComplexHeatmap::rowAnnotation(
    PC1 = top_loadings$PC1_member %>% as.character() %>% anno_simple(col = c("TRUE" = "black", "FALSE" = "white")),
    PC2 = top_loadings$PC2_member %>% as.character() %>% anno_simple(col = c("TRUE" = "black", "FALSE" = "white")),
    PC3 = top_loadings$PC3_member %>% as.character() %>% anno_simple(col = c("TRUE" = "black", "FALSE" = "white")),
    PC4 = top_loadings$PC4_member %>% as.character() %>% anno_simple(col = c("TRUE" = "black", "FALSE" = "white"))
  )
  maybe_metadata <- pca_object$metadata
  if ("dummy" %in% colnames(maybe_metadata)) {
    maybe_metadata <- NULL
  }

  ht <- plot_tools$plot_results_one_collection(
    df = submat,
    metadata = maybe_metadata, # guaranteed to match the colnames of df which was used originally to make the pca obj
    # pathway_metadata = loadings,
    row_annotation = ra,
    limit = Inf,
    pstat_cutoff = 1, # we want no filtering
    save_func = save_func,
    title = "Top 5 Loadings Per PC",
    ...
  )
}
