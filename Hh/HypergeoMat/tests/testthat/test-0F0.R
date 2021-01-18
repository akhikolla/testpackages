context("0F0")

test_that("0F0 is the exponential of the trace", {
  X <- toeplitz(c(3,2,1))/10
  obtained <- hypergeomPFQ(m = 10, a = NULL, b = NULL, x = X)
  expected <- exp(sum(diag(X)))
  expect_equal(obtained, expected)
  #
  X <- toeplitz(c(3i,2,1))/10
  obtained <- hypergeomPFQ(m = 10, a = NULL, b = NULL, x = X)
  expected <- exp(sum(diag(X)))
  expect_equal(obtained, expected)
})

