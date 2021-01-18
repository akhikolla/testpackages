library(testthat)
library(AnaCoDa)

context("Utility")

test_that("general utility functions", {
  expect_equal(testUtility(), 0)
})
