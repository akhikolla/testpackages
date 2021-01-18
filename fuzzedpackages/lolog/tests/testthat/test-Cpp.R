
context("C++ tests")

#library(lolog)
library(testthat)



test_that("C++", {
  runLologCppTests()
  expect_true(TRUE)
})
