suppressPackageStartupMessages(library(rlang))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))



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

get_config <- function(file_path) {
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

  if (is.null(config$params$advanced)) {
    config$params$advanced <- list()
  }
  config$params$advanced$replace <- is_null_config(config$params$advanced$replace) %||% FALSE


  return(config$params)
}
