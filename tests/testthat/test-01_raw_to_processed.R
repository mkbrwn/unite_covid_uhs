# Test script for UNITE COVID data analysis

library(testthat)

# Load the processed data
source("src/01_raw_to_processed.R")

# Test that the processed data is not empty
test_that("Processed data is not empty", {
  expect_true(nrow(UNITE_2020_corrected) > 0)
  expect_true(nrow(UNITE_2021_corrected) > 0)
})

# Test that Exposure(Corticosteroid and outcome (60 day mortality) data is not NA
test_that("Corticosteroid and outcome (60 day mortality) var is not NA", {
  expect_true(all(!is.na(UNITE_2020_corrected$ICU_CORTICO_YN)))
  expect_true(all(!is.na(UNITE_2020_corrected$OUTCOME_LD)))
})

# Test that Corticosteroids has 2 levels
test_that("Corticosteroids has 2 levels", {
  expect_equal(length(unique(UNITE_2020_corrected$ICU_CORTICO_YN)), 2)
  })
