# test_utils.R
library(testthat)
library(here)


print(here())
#source(here("R/lazyloader.R"))

testthat::test_that("test_utils", {
  source(here("R/utils.R"))
  expect_true(TRUE)
})


test_that("impute_with_draw handles no NA input", {
  expect_equal(impute_with_draw(c(1, 2, 3)), c(1, 2, 3))
})

test_that("impute_with_draw handles all NA input", {
  expect_error(impute_with_draw(c(NA, NA, NA)))
})


test_that("impute_with_draw handles some NA input", {
  expect_no_error(impute_with_draw(c(1, NA, 3)))
})

test_that("impute_with_draw handles a large matrix", {
  # Set seed for reproducibility
  set.seed(123)

  # Create a 10x10 matrix with some NA values
  test_matrix <- matrix(rnorm(100, mean = 5, sd = 2), nrow = 10)
  test_matrix[sample(1:100, size = 20)] <- NA  # Introduce 20 random NAs

  # Ensure there are NA values
  expect_true(any(is.na(test_matrix)))

  # Impute the matrix
  # imputed_matrix <- apply(test_matrix, 1, function(row) impute_with_draw(row))
  imputed_matrix <- impute_with_draw(test_matrix)
  imputed_matrix <- t(imputed_matrix)  # Transpose back to original dimensions

  # Check that the imputed matrix has no NAs
  expect_false(any(is.na(imputed_matrix)))

  # Check that imputed values are not uniform across the matrix
  unique_values <- unique(imputed_matrix)
  expect_true(length(unique_values) > 1)

  # Check that the imputed values differ from the original NAs
  for (i in 1:nrow(test_matrix)) {
    for (j in 1:ncol(test_matrix)) {
      if (is.na(test_matrix[i, j])) {
        # Ensure the imputed value is not NA
        expect_false(is.na(imputed_matrix[i, j]))
        # Ensure the imputed value differs from a uniform value (e.g., all the same)
        expect_false(all(imputed_matrix[i, ] == imputed_matrix[i, j]))
      }
    }
  }

  # Check if imputed values maintain the general scale of the data
  expect_true(all(imputed_matrix >= min(test_matrix, na.rm = TRUE) - 3))
  expect_true(all(imputed_matrix <= max(test_matrix, na.rm = TRUE) + 3))

  # # Optional visualization (useful during development)
  # print("Original Matrix with NAs:")
  # print(test_matrix)
  # print("Imputed Matrix:")
  # print(imputed_matrix)

})