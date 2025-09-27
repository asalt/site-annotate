suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(withr))
suppressPackageStartupMessages(library(here))

source_cache <- function(tmp_dir) {
  cache_env <- new.env(parent = parent.frame())
  withr::local_dir(tmp_dir)
  source(here("R/cache.R"), local = cache_env)
  cache_env
}

setup_counter <- function(cache_env) {
  cache_env$counter <- 0
  cache_env$inc <- function() {
    cache_env$counter <- cache_env$counter + 1
    cache_env$counter
  }
}


test_that("cached_by stores small objects in configured DB", {
  tmp <- withr::local_tempdir()
  db_path <- file.path(tmp, "objects.sqlite")
  rds_dir <- file.path(tmp, "rds")

  withr::local_options(
    cached_by_db_path = db_path,
    cached_by_rds_dir = rds_dir
  )

  cache_env <- source_cache(tmp)
  setup_counter(cache_env)

  res1 <- evalq(inc() %cached_by% "db-test", envir = cache_env)
  expect_equal(res1, 1)
  expect_true(file.exists(db_path))
  expect_true(dir.exists(rds_dir))

  res2 <- evalq(inc() %cached_by% "db-test", envir = cache_env)
  expect_equal(res2, 1)
  expect_equal(cache_env$counter, 1)
})


test_that("cached_by writes to RDS when objects exceed threshold", {
  tmp <- withr::local_tempdir()
  db_path <- file.path(tmp, "objects.sqlite")
  rds_dir <- file.path(tmp, "rds")

  withr::local_options(
    cached_by_db_path = db_path,
    cached_by_rds_dir = rds_dir,
    cached_by_db_threshold_bytes = 1L  # ensure everything goes to RDS
  )

  cache_env <- source_cache(tmp)
  setup_counter(cache_env)

  res1 <- evalq({ inc(); rep(1, 10) } %cached_by% "rds-test", envir = cache_env)
  expect_equal(cache_env$counter, 1)
  files <- list.files(rds_dir, pattern = "\\.rds$", full.names = TRUE)
  expect_length(files, 1)
  expect_true(file.exists(files[[1]]))

  res2 <- evalq({ inc(); rep(1, 10) } %cached_by% "rds-test", envir = cache_env)
  expect_equal(res2, res1)
  expect_equal(cache_env$counter, 1)
})


test_that("cached_by falls back to RDS when DB write fails", {
  tmp <- withr::local_tempdir()
  db_path <- file.path(tmp, "objects.sqlite")
  rds_dir <- file.path(tmp, "rds")

  withr::local_options(
    cached_by_db_path = db_path,
    cached_by_rds_dir = rds_dir
  )

  cache_env <- source_cache(tmp)
  setup_counter(cache_env)

  cache_env$save_object_to_cache <- function(...) stop("boom")

  expect_message(
    res1 <- evalq({ inc(); letters[1:3] } %cached_by% "fallback-test", envir = cache_env),
    "falling back to RDS",
    fixed = TRUE
  )
  expect_equal(cache_env$counter, 1)
  files <- list.files(rds_dir, pattern = "\\.rds$", full.names = TRUE)
  expect_length(files, 1)

  res2 <- evalq({ inc(); letters[1:3] } %cached_by% "fallback-test", envir = cache_env)
  expect_equal(res2, res1)
  expect_equal(cache_env$counter, 1)
})
