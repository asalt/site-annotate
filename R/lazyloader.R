# lazyloader.R
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(fs))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(magrittr))

# Only define get_tool_env if it doesn't already exist
if (!exists("get_tool_env", envir = .GlobalEnv)) {
  get_tool_env <- local({
    # Resolve R tooling directory with a robust fallback
    root_R_dir <- tryCatch(here::here("R"), error = function(e) "R")
    # Build tool list without relying on magrittr pipes or here() at call sites
    files <- fs::dir_ls(path = root_R_dir, glob = "*.R")
    all_tools <- basename(files)
    all_tools <- gsub("\\.R$", "", all_tools)

    print(getwd())
    print(all_tools)

    # all_tools <- c("heatmap", "io", "utils")
    #heatmap.R  io.R  lazyloader.R  reduce.R  summarize_rmd.R  utils.R

    # Initialize cache with empty environments
    tools_cache <- vector("list", length(all_tools))
    names(tools_cache) <- all_tools
    for (tool_name in all_tools) {
      tools_cache[[tool_name]] <- new.env()
    }

    # Keep track of tools currently being loaded
    tools_loading <- list()

    function(tool_name) {
      if (!tool_name %in% names(tools_cache)) {
        stop("Unknown tool: ", tool_name)
      }

      # Return environment if already loaded
      if (exists(".__loaded__", envir = tools_cache[[tool_name]])) {
        return(tools_cache[[tool_name]])
      }

      # Return environment if already loading (circular dependency)
      if (tool_name %in% tools_loading) {
        return(tools_cache[[tool_name]])
      }

      # Mark as currently loading
      tools_loading <<- c(tools_loading, tool_name)

      # Source the tool file into its environment
      src_dir <- root_R_dir
      source_file <- paste0(tool_name, ".R")
      full_path <- file.path(src_dir, source_file)
      if (!file.exists(full_path)) {
        tools_loading <<- setdiff(tools_loading, tool_name)
        stop("Source file does not exist for tool: ", tool_name)
      }

      # Provide access to get_tool_env
      env <- tools_cache[[tool_name]]
      env$get_tool_env <- get_tool_env
      env$tools_cache <- tools_cache

      # Source the file
      source(full_path, local = env)

      # Mark as loaded
      assign(".__loaded__", TRUE, envir = env)

      # Remove from loading state
      tools_loading <<- setdiff(tools_loading, tool_name)

      return(env)
    }
  })

  # Assign get_tool_env to the global environment
  assign("get_tool_env", get_tool_env, envir = .GlobalEnv)
}

# Also ensure get_tool_env is bound in the current sourcing environment
if (!exists("get_tool_env", envir = environment())) {
  if (exists("get_tool_env", envir = .GlobalEnv)) {
    assign("get_tool_env", get("get_tool_env", envir = .GlobalEnv), envir = environment())
  }
}
