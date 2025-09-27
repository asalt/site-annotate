# setup_environment.R


# load_from_cache <- function() {
#   source("cache.R")
#   load_environment_from_db(db_path)
# }

# setup_environment <- function() {

# load_from_cache()
# source("cache.R")

# load_environment_from_db(db_path)

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
# #
# source("lazyloader.R")
# #
# return()
# }



setup_environment <- function(
    db_path = "cache.sqlite",
    packages_to_load = c(
      "rlang", "tidyverse", "magrittr", "purrr",
      "ComplexHeatmap", "circlize", # "reactable",
      "tibble", "grid",
      "RColorBrewer", "here"
    )) {
  log_tools <- tryCatch(get_tool_env("logging"), error = function(e) NULL)
  if (is.null(log_tools)) {
    log_info <- function(...) invisible(NULL)
    log_warn <- function(...) invisible(NULL)
    log_error <- function(...) invisible(NULL)
  } else {
    log_info <- log_tools$log_info
    log_warn <- log_tools$log_warn
    log_error <- log_tools$log_error
  }

  log_info("Initializing R environment", context = list(db_path = db_path, packages = length(packages_to_load)))
  source("cache.R")
  con <- initialize_cache_db(db_path, close = FALSE)
  # con <- initialize_cache_db(db_path, close = T)
  # con <- NULL

  for (package in packages_to_load) {
    log_info("Loading package", context = list(package = package))
    tryCatch(
      {
        handle_package_cache(package, db_path, con = con)
      },
      error = function(e) {
        log_warn("Failed to cache package", context = list(package = package, error = e$message))
      }
    )
  }

  # now check if packages_to_load are loaded via .packages() call
  loaded_pkgs <- .packages()
  for (pkg in packages_to_load) {
    if (!(pkg %in% loaded_pkgs)) {
      log_error("Package failed to load", context = list(package = pkg))
      stop()
    }
  }

  # browser()

  # close the connection
  if (!is.null(con)) {
    DBI::dbDisconnect(con)
    log_info("Closed cache connection", context = list(db_path = db_path))
  }

  log_info("Environment setup complete", context = list(packages_loaded = length(packages_to_load)))
}
