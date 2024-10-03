suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(lightparser))

do_summary <- function(file, ...){
    if (!file.exists(file)) stop(paste0("file ", file, " does not exist!"))
    tbl <- lightparser::split_to_tbl(file)
}

