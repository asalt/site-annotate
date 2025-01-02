suppressPackageStartupMessages(library(rlang))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(cmapR))
suppressPackageStartupMessages(library(janitor))

is_absolute_path <- function(path) {
  # On Unix-like systems, an absolute path starts with '/'
  # On Windows, it usually starts with a drive letter, like 'C:/'
  if (.Platform$OS.type == "windows") {
    return(grepl("^[A-Za-z]:[\\/]", path))
  } else {
    return(startsWith(path, "/"))
  }
}
# Function to convert a file path to absolute if it's not already
ensure_absolute_path <- function(file_path, root_dir) {
  if (!is_absolute_path(file_path)) {
    return(normalizePath(file.path(root_dir, file_path), mustWork = FALSE))
  }
  return(file_path)
}


convert_tmt_label <- function(shorthand) {
  normalized_input <- str_remove_all(shorthand, "TMT[_]?")
  res <- case_when(
    normalized_input == "126" ~ "TMT_126",
    normalized_input == "131" ~ "TMT_131_N",
    normalized_input == "134" ~ "TMT_134_N",
    grepl("N$", normalized_input) ~ sub("(\\d+)_?(N)", "TMT_\\1_N", normalized_input),
    grepl("C$", normalized_input) ~ sub("(\\d+)_?(C)", "TMT_\\1_C", normalized_input),
    TRUE ~ paste0("TMT_", normalized_input)
  )
  return(res)
}

set_name <- function(f) {
  path_split <- f %>%
    fs::path_ext_remove() %>%
    fs::path_split() %>%
    unlist()
  path_split_len <- length(path_split)
  if (path_split_len == 1) {
    res <- path_split[1]
  } else {
    res <- paste(path_split[path_split_len - 1] %>% make.names(),
      path_split[path_split_len] %>% make.names(),
      sep = "_"
    )
  }
  res <- res %>% stringr::str_remove(("\\.$"))
  return(res)
}


extract_key_value_pairs <- function(input_string) {
  # Split the string by '|'
  elements <- unlist(strsplit(input_string, "\\|"))

  # Initialize empty lists for keys and values
  keys <- character()
  values <- character()

  # Iterate over the elements to form key-value pairs
  for (i in seq_along(elements)) {
    # When i is odd, it's a key, when even, it's a value
    if (i %% 2 == 1) {
      # Next element should be value unless we are at the end or next is a key
      if (i < length(elements) && !(elements[i + 1] %in% c("ENSP", "ENST", "ENSG", "geneid", "taxon", "symbol"))) {
        keys <- c(keys, elements[i])
        values <- c(values, elements[i + 1])
      } else {
        # Handle cases where a key might not have a following value or next is a key
        next()
        # keys <- c(keys, elements[i])
        # values <- c(values, NA)
      }
    }
  }

  # Combine into a named vector
  named_vector <- setNames(values, keys)

  return(named_vector)
}

is_null_config <- function(val, nullstrings = c("None", "")) {
  if (is.null(val)) return(NULL)
  if (all(val %in% nullstrings)) return(NULL)
  return(val)
}

load_config <- function(file_path, root_dir = NULL) {
  # Read the config file
  if (!is.null(file_path) && file_path == "None") {
    file_path <- NULL
  }

  config_file <-  file_path %||% file.path(here("config"), "base.toml")
  config <- RcppTOML::parseTOML(config_file)

  if (is.null(config$params)) {
    warn("No parameters found in the config file.")
    return()
  }

  config$params$batch <- is_null_config(config$params$batch)
  # if (config$params$batch == "None" || config$params$batch == "") {
  #   config$params$batch <- NULL
  # }

  config$params$limma$formula <- is_null_config(config$params$limma$formula)
  if (!is.null(config$params$limma$formula)) {
    config$params$limma$formula %<>% as.formula
  }

  config$params$limma$contrasts <- is_null_config(config$params$limma$contrasts)

  config$params$heatmap$cut_by <- is_null_config(config$params$heatmap$cut_by)
  config$params$heatmap$z_score <- config$params$heatmap$z_score %||% TRUE
  config$params$heatmap$z_score_by %<>% is_null_config(.) %||% NULL


  config$params$norm$normalize <- config$params$norm$normalize %||% FALSE
  config$params$norm$log_transform <- config$params$norm$log_transform %||% FALSE

  config$params$totalprotein$do <- config$params$totalprotein$do %||% FALSE

  if (is.null(config$params$advanced)) {
    config$params$advanced <- list()
  }
  config$params$advanced$replace <- is_null_config(config$params$advanced$replace) %||% FALSE


  if (is.null(root_dir)) {
    root_dir <- getwd()
  }

  # handle any lists of files that may have relative paths
  if (is.null(config$params$heatmap$selection)){
    config$params$heatmap$selection <- list()
  }

  if (!is.null(config$params$heatmap$selection$file_list)){
    # iterate through each one and make sure the paths are defined as absolute
    # the paths are the values of the list, the keys are the names for use for each file
    # Use purrr::map to iterate through the file_list and convert paths to absolute
    config$params$heatmap$selection$file_list <- map(
      config$params$heatmap$selection$file_list,
      function(file_entry) {
        # Apply ensure_absolute_path function to each path in the file_entry
        map(file_entry, ~ ensure_absolute_path(.x, root_dir))
      }
    )
    config$params$heatmap$selection$file_list <- purrr::list_flatten(config$params$heatmap$selection$file_list)
  }


  config$params$advanced$replace <- config$params$advanced$replace  %||% FALSE

  return(config$params)
}

read_expr <- function(.x){
  df <- readr::read_tsv(.x,
                        #  n_max=4000,
                        show_col_types = FALSE
                        ) %>%
    mutate(expr_file = .x)
  # mutate(modi_id = basename(id) %>% str_extract(pattern = "(sty|k)_[0-9]+_[0-9]+")) %>%

  df %<>% mutate(sitename = if_else(is.na(sitename2), sitename, sitename2))
  if ("geneid" %in% colnames(df)) df %<>% mutate(geneid = as.character(geneid))

  if ("cst_cat#" %in% colnames(df)) df %<>% mutate(`cst_cat#` = as.character(`cst_cat#`))

  if ("all_possible_positions" %in% colnames(df)) df$all_possible_positions <- NULL
  if ("all_possible_positions_psp" %in% colnames(df)) df$all_possible_positions_psp <- NULL

  if ("ms_lit" %in% colnames(df)) df %<>% mutate(`ms_lit` = replace_na(ms_lit, 0))
  if ("lt_lit" %in% colnames(df)) df %<>% mutate(`lt_lit` = replace_na(lt_lit, 0))

  if (!"fifteenmer" %in% colnames(df)) stop('fifteenmer missing')
  if (!"ENSP" %in% colnames(df)) stop('ENSP missing, this is a protein id ')
  if (!"sitename" %in% colnames(df)) stop('sitename missing')
  if (!"taxon" %in% colnames(df)) df$taxon <- "<NULL>"
  df %<>% mutate(site_id = str_c(sitename, ENSP, fifteenmer, sep='_'))

  return(df)
}

maybe_match_arg <- function(arg, choices, several.ok = TRUE) {
    tryCatch({
      match.arg(arg, choices, several.ok = several.ok)
    },
    error = function(e) {
      return(NULL)
  })
}

get_file_reader <- function(the_file){
  # first get the extension, then choose appropriat reader
  ext <- tools::file_ext(the_file)
  if (ext == "tsv" || ext == "txt"){
    reader <- readr::read_tsv
  }
  if (ext == "csv"){
    reader <- readr::read_csv
  }
  if (ext == "xlsx"){
    reader <- readxl::read_xlsx
  }
  else {
    reader <- readr::read_tsv
  }
  return(reader)
}

read_genelist_file <- function(file_path,
                               id_col = "id",
                                possible_id_cols = c("id", "GeneID", "gene"),
                                possible_contrast_cols = c("contrast", "group"),
                                 ...){
  # Read the file
  reader <- get_file_reader(file_path)
  df <- reader(file_path, show_col_types = FALSE) %>% janitor::clean_names()

  if (!id_col %in% colnames(df)) {
    for (possible_id_col in possible_id_cols) {
      if (possible_id_col %in% colnames(df)) {
        id_col <- possible_id_col
        df[["id"]] <- df[[id_col]]
        break
      }
    }
  }

  # now sort by another column if it exists
  possible_col_for_sort <-maybe_match_arg(colnames(df), c("pvalue", "p_value"), several.ok = TRUE)
  if (!is.null(possible_col_for_sort)) {
    if (length(possible_col_for_sort) > 1) possible_col_for_sort <- possible_col_for_sort[1]
    df <- df %>% arrange(!!sym(possible_col_for_sort))
  }

  possible_col_for_facet <- maybe_match_arg(colnames(df), possible_contrast_cols, several.ok = TRUE)
  if (!is.null(possible_col_for_facet)) {
    if (length(possible_col_for_facet) > 1) possible_col_for_facet <- possible_col_for_facet[1]
    vals <- unique(df[[possible_col_for_facet]])
    df <- purrr::map(vals,
      ~ df %>% filter(!!sym(possible_col_for_facet) == .x)
    ) %>% set_names(vals)
  }

  return(df)


}


create_gct_from_datal <- function(datal){

    # 2.1. Create the Expression Matrix
    # Select necessary columns: 'site_id' and expression values
    # Here, 'name' represents the sample identifier

    expression_matrix <- datal %>%
      select(site_id, name, log2_val) %>%
      pivot_wider(names_from = name,
                  values_from = log2_val,
                  values_fn = function(x) median(x, na.rm = T) # sometimes the merging of different tables creates duplicate records
                  # e.g. acquired data merged with phosphositeplus .

      ) %>%
      column_to_rownames("site_id") %>%
      as.matrix()

    row_metadata <- datal %>%
      select(site_id,
             sitename,
             lt_lit, ms_lit, ms_cst,
             gene, acc_id, hu_chr_loc, mod_rsd, organism, mw_kd, domain, fifteenmer, ENSP, ENST, ENSG, geneid, taxon, symbol,
             AA,
             protein_start,
             protein_end,
             protein_start_psp
             ) %>%
      distinct(site_id, .keep_all = TRUE) %>%
      column_to_rownames("site_id")
    row_metadata$id <- rownames(row_metadata)


    # Example:
    column_metadata <- meta %>%
      # select(name, any_of(c("cell", "time", "treatment", "group", "replicate"))) %>%
      distinct(name, .keep_all = TRUE) %>%
      column_to_rownames("name")
    column_metadata$id <- rownames(column_metadata)

    # 3.1. Create the GCT Object
    gct_object <- new("GCT",
      mat = expression_matrix,
      rid = rownames(expression_matrix),
      cid = colnames(expression_matrix),
      rdesc = row_metadata,
      cdesc = column_metadata
    )
    return(gct_object)
}


match_comparisons <- function(dir1, dir2){
  dir1_files <- fs::dir_ls(dir1, glob = "*.tsv")
  dir2_files <- fs::dir_ls(dir2, glob = "*.tsv")

  if (length(dir1_files) == 0) {
    warning("No files found in dir1. skipping")
    return()
  }

  if (length(dir2_files) == 0) {
    warning("No files found in dir2. skipping")
    return()
  }

  # example file names

}
