

library(testthat)
library(imputeMulti)

context("multinomial_stats")

test_that("possible obs returns correctly", {
  data(tract2221)
  x_p <- multinomial_stats(tract2221[,1:5], "possible.obs")
  
  expect_equal(class(x_p), "data.frame")
  expect_equal(names(x_p), names(tract2221)[1:5])
  expect_true(all(unlist(lapply(x_p, is.factor))))
  expect_true(all(levels(x_p[[1]]) %in% levels(tract2221[[1]])))
  expect_true(all(levels(x_p[[2]]) %in% levels(tract2221[[2]])))
  expect_true(all(levels(x_p[[3]]) %in% levels(tract2221[[3]])))
  expect_true(all(levels(x_p[[4]]) %in% levels(tract2221[[4]])))
  expect_true(all(levels(x_p[[5]]) %in% levels(tract2221[[5]])))
})

test_that("possible obs returns correctly", {
  data(tract2221)
  x_y <- multinomial_stats(tract2221[,1:5], "x_y")
  
  expect_equal(class(x_y), "data.frame")
  expect_equal(names(x_y), c(names(tract2221)[1:5], "counts"))
  expect_true(all(unlist(lapply(x_y[,-6], is.factor))))
  expect_true(all(levels(x_y[[1]]) %in% levels(tract2221[[1]])))
  expect_true(all(levels(x_y[[2]]) %in% levels(tract2221[[2]])))
  expect_true(all(levels(x_y[[3]]) %in% levels(tract2221[[3]])))
  expect_true(all(levels(x_y[[4]]) %in% levels(tract2221[[4]])))
  expect_true(all(levels(x_y[[5]]) %in% levels(tract2221[[5]])))
  
  expect_equal(sum(x_y$counts), sum(complete.cases(tract2221[,1:5])))
})

test_that("possible obs returns correctly", {
  data(tract2221)
  x_y <- multinomial_stats(tract2221[,1:5], "z_Os_y")
  
  expect_equal(class(x_y), "data.frame")
  expect_equal(names(x_y), c(names(tract2221)[1:5], "counts"))
  expect_true(all(unlist(lapply(x_y[,-6], is.factor))))
  expect_true(all(levels(x_y[[1]]) %in% levels(tract2221[[1]])))
  expect_true(all(levels(x_y[[2]]) %in% levels(tract2221[[2]])))
  expect_true(all(levels(x_y[[3]]) %in% levels(tract2221[[3]])))
  expect_true(all(levels(x_y[[4]]) %in% levels(tract2221[[4]])))
  expect_true(all(levels(x_y[[5]]) %in% levels(tract2221[[5]])))
  
  expect_equal(sum(x_y$counts), sum(!complete.cases(tract2221[,1:5])))
})