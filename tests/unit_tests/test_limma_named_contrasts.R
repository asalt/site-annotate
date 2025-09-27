suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(here))

load_limma_env <- function() {
  env <- new.env(parent = baseenv())
  sys.source(here("R", "limma.R"), envir = env)
  env
}

load_run_env <- function() {
  env <- new.env(parent = baseenv())
  sys.source(here("R", "run.R"), envir = env)
  env
}

build_test_gct <- function() {
  skip_if_not_installed("cmapR")
  mat <- matrix(
    c(
      2, 4, 6, 8,
      1, 3, 5, 7
    ),
    nrow = 2,
    byrow = TRUE
  )
  rid <- paste0("gene", seq_len(nrow(mat)))
  cid <- paste0("sample", seq_len(ncol(mat)))
  group <- rep(c("A", "B"), each = 2)
  cdesc <- data.frame(group = group, row.names = cid)
  rdesc <- data.frame(id = rid, row.names = rid)
  cmapR::GCT(mat = mat, rid = rid, cid = cid, cdesc = cdesc, rdesc = rdesc)
}


with_limma_env <- function(code) {
  env <- load_limma_env()
  eval(substitute(code), envir = env)
}


test_that("parse_named_contrasts handles named and unnamed entries", {
  result <- with_limma_env({
    parse_named_contrasts(c(
      "friendly=groupA - groupB",
      "groupC - groupD",
      "duplicate=groupC - groupD"
    ))
  })

  expect_equal(
    result$label,
    c("friendly", "groupC - groupD", "duplicate")
  )
  expect_equal(
    result$expression,
    c("groupA - groupB", "groupC - groupD", "groupC - groupD")
  )
})


create_limma_config <- function() {
  list(
    limma = list(
      formula = as.formula("~0 + group"),
      contrasts = c("friendly=groupA - groupB"),
      group_var = "group",
      top_n = c(2)
    )
  )
}


test_that("run_limma preserves label and expression metadata", {
  skip_if_not_installed("limma")
  gct <- build_test_gct()
  config <- create_limma_config()

  output <- with_limma_env({
    run_limma(gct, config, do_impute = FALSE)
  })

  expect_named(output, "friendly")
  expect_true(all(output[["friendly"]]$contrast == "friendly"))
  expect_true(all(output[["friendly"]]$contrast_expression == "groupA - groupB"))

  contrast_attr <- attr(output, "contrast_expression")
  expect_equal(contrast_attr[["friendly"]], "groupA - groupB")
})


test_that("sanitize_path_component removes unsafe characters", {
  run_env <- load_run_env()
  sanitize <- run_env$sanitize_path_component

  expect_equal(sanitize("A vs B"), "A_vs_B")
  expect_equal(sanitize("F|oo"), "F_oo")
  expect_equal(sanitize("  leading"), "leading")
  expect_equal(sanitize(""), "unnamed")
})
