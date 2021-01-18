library(testthat)
library(combiter)
library(foreach)

context("use with foreach")

test_that("foreach and as.list", {
  expect_equal(foreach(i = iperm(5)) %do% i, as.list(iperm(5)))

  expect_equal(foreach(i = icomb(7,3)) %do% sum(i),
               lapply(as.list(icomb(7,3)), sum))

  expect_equal(foreach(i = isubset(6)) %do% length(i),
               lapply(as.list(isubset(6)), length))
})


