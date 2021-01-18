library(testthat)
library(AnaCoDa)

context("Parameter")

test_that("general parameter functions", {
  expect_equal(testParameter("UnitTestingData/testMCMCROCFiles"), 0)
})

# TODO: Implement the following
#p <-new(Parameter)
