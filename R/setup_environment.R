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
    packages_to_load = c("rlang", "tidyverse", "magrittr", "purrr",
                         "ComplexHeatmap", "circlize", #"reactable",
                         "tibble",
                         "RColorBrewer", "here")
) {

  source("cache.R")
  con <- initialize_cache_db(db_path, close = FALSE)
  # con <- initialize_cache_db(db_path, close = T)
  # con <- NULL

  for (package in packages_to_load) {
    message(paste("Processing package:", package))
    tryCatch({
      handle_package_cache(package, db_path, con=con)
    }, error = function(e) {
      message(paste("Failed to process package:", package, "-", e$message))
    })
  }

  # now check if packages_to_load are loaded via .packages() call
  loaded_pkgs <- .packages()
  for (pkg in packages_to_load) {
    if (!(pkg %in% loaded_pkgs)) {
      message(paste("Package", pkg, "not loaded."))
      stop()
    }
  }

# browser()

  # close the connection
  if (!is.null(con)) {
    DBI::dbDisconnect(con)
    message("Connection closed.")
  }

  message("Environment setup complete.")
}