# logging.R
# Lightweight structured logging utilities for site-annotate R workflows.

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) {
    y
  } else {
    x
  }
}

log_level_rank <- c(
  TRACE = 10,
  DEBUG = 20,
  INFO = 30,
  SKIP = 30,
  WARN = 40,
  ERROR = 50,
  FATAL = 60
)

#' Determine whether a log level should be emitted under the current threshold.
log_is_enabled <- function(level) {
  level <- toupper(level %||% "INFO")
  lvl_rank <- log_level_rank[[level]] %||% log_level_rank[["INFO"]]
  threshold_level <- toupper(getOption("site_annotate.log_level", "INFO"))
  threshold_rank <- log_level_rank[[threshold_level]] %||% log_level_rank[["INFO"]]
  lvl_rank >= threshold_rank
}

#' Format arbitrary R objects for log context.
format_log_value <- function(value) {
  if (is.null(value) || length(value) == 0) {
    return("<NA>")
  }
  if (is.list(value)) {
    value <- unlist(value, recursive = TRUE, use.names = FALSE)
  }
  if (length(value) > 1) {
    value <- paste(value, collapse = ",")
  }
  if (is.logical(value)) {
    return(ifelse(value, "TRUE", "FALSE"))
  }
  if (is.numeric(value)) {
    return(format(value, digits = 6, trim = TRUE))
  }
  value <- as.character(value)
  value <- encodeString(value, quote = "")
  value
}

#' Turn a context list into a string of key=value pairs.
format_log_context <- function(context) {
  if (is.null(context) || length(context) == 0) {
    return("")
  }
  if (is.null(names(context))) {
    names(context) <- rep("", length(context))
  }
  pieces <- Map(function(name, value) {
    formatted <- format_log_value(value)
    if (!nzchar(name)) {
      return(formatted)
    }
    sprintf("%s=%s", name, formatted)
  }, names(context), context)
  paste(unlist(pieces), collapse = " ")
}

#' Core logging function. Use the convenience wrappers below instead.
log_emit <- function(level = "INFO", msg, ..., context = NULL) {
  if (!log_is_enabled(level)) {
    return(invisible(NULL))
  }
  if (length(list(...)) > 0) {
    msg <- sprintf(msg, ...)
  }
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  ctx <- format_log_context(context)
  prefix <- sprintf("[%s] [%s]", timestamp, toupper(level))
  line <- trimws(paste(prefix, msg, ctx))

  message(line)

  log_file <- getOption("site_annotate.log_file")
  if (!is.null(log_file) && nzchar(log_file)) {
    tryCatch(
      {
        dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
        cat(line, "\n", file = log_file, append = TRUE)
      },
      error = function(e) {
        message(sprintf("[site-annotate][WARN] Failed to write log file %s: %s", log_file, e$message))
      }
    )
  }
  invisible(line)
}

log_trace <- function(msg, ..., context = NULL) {
  log_emit("TRACE", msg, ..., context = context)
}

log_debug <- function(msg, ..., context = NULL) {
  log_emit("DEBUG", msg, ..., context = context)
}

log_info <- function(msg, ..., context = NULL) {
  log_emit("INFO", msg, ..., context = context)
}

log_skip <- function(msg, ..., context = NULL) {
  log_emit("SKIP", msg, ..., context = context)
}

log_warn <- function(msg, ..., context = NULL) {
  log_emit("WARN", msg, ..., context = context)
}

log_error <- function(msg, ..., context = NULL) {
  log_emit("ERROR", msg, ..., context = context)
}

set_log_level <- function(level = "INFO") {
  options(site_annotate.log_level = toupper(level))
  invisible(level)
}

set_log_file <- function(path) {
  options(site_annotate.log_file = path)
  invisible(path)
}

# Export helpers for other modules
assign("log_is_enabled", log_is_enabled, envir = environment())
assign("log_emit", log_emit, envir = environment())
assign("log_trace", log_trace, envir = environment())
assign("log_debug", log_debug, envir = environment())
assign("log_info", log_info, envir = environment())
assign("log_skip", log_skip, envir = environment())
assign("log_warn", log_warn, envir = environment())
assign("log_error", log_error, envir = environment())
assign("set_log_level", set_log_level, envir = environment())
assign("set_log_file", set_log_file, envir = environment())
