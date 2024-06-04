suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(withr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tibble))
# suppressPackageStartupMessages(library(cmapR))
suppressPackageStartupMessages(library(here))

src_dir <- file.path(here("R"))
io_tools <- new.env()
source(file.path(src_dir, "./io.R"), local = io_tools)


test_that("convert_tmt_label_test_equal_12xN|C", {
  io_tools$convert_tmt_label("TMT_127_N") %>% expect_equal("TMT_127_N")
  io_tools$convert_tmt_label("TMT_127_C") %>% expect_equal("TMT_127_C")
})

test_that("convert_tmt_label_test_equal_126", {
  io_tools$convert_tmt_label("TMT_126") %>% expect_equal("TMT_126")
})

test_that("convert_tmt_label_test_equal_131", {
  io_tools$convert_tmt_label("TMT_131") %>% expect_equal("TMT_131_N")
})

test_that("convert_tmt_label_test_equal_134", {
  io_tools$convert_tmt_label("TMT_134") %>% expect_equal("TMT_134_N")
})


test_that("convert_tmt_label_test_short", {
  io_tools$convert_tmt_label("127_N") %>% expect_equal("TMT_127_N")
  io_tools$convert_tmt_label("127_C") %>% expect_equal("TMT_127_C")
  io_tools$convert_tmt_label("126") %>% expect_equal("TMT_126")
  io_tools$convert_tmt_label("131") %>% expect_equal("TMT_131_N")
  io_tools$convert_tmt_label("131_C") %>% expect_equal("TMT_131_C")
})

test_that("convert_tmt_label_test_var", {
  io_tools$convert_tmt_label("127N") %>% expect_equal("TMT_127_N")
  io_tools$convert_tmt_label("127C") %>% expect_equal("TMT_127_C")
  io_tools$convert_tmt_label("126") %>% expect_equal("TMT_126")
  io_tools$convert_tmt_label("131C") %>% expect_equal("TMT_131_C")
})


#
test_that("set_name_test", {
  io_tools$set_name("test.txt") %>% expect_equal("test")
})

#
test_that("set_name_withparent", {
  io_tools$set_name("parent/test.txt") %>% expect_equal("parent_test")
})

#
test_that("set_name_withparents", {
  io_tools$set_name("grandparent/parent/test.txt") %>% expect_equal("parent_test")
})

#
test_that("set_name_withparents", {
  io_tools$set_name("parent/test.") %>% expect_equal("parent_test")
})

# == extract key value pairs ==

test_that("extract_key_value_pairs_test", {
  input_string <- "ENSP|ENSMUSP00000001963|ENST|ENSMUST00000001963|ENSG|ENSMUSG00000020681|geneid|11421|taxon|10090|symbol|Ace|Ace"
  expected_names <- c("ENSP", "ENST", "ENSG", "geneid", "taxon", "symbol")
  expected_values <- c("ENSMUSP00000001963", "ENSMUST00000001963", "ENSMUSG00000020681", "11421", "10090", "Ace")
  res <- io_tools$extract_key_value_pairs(input_string)
  for (name in expected_names) {
    rlang::has_name(res, name) %>% expect_true()
  }
  for (value in expected_values) {
    expect_true(value %in% res)
  }
})
