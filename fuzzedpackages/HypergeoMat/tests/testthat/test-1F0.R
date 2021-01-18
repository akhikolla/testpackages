context("1F0")

test_that("1F0 is det(I-X)^(-a)", {
  X <- toeplitz(c(3,2,1))/100
  obtained <- hypergeomPFQ(m = 15, a = 3, b = NULL, x = X)
  expected <- det(diag(3)-X)^(-3)
  expect_equal(obtained, expected)
  #
  X <- toeplitz(c(2,2,1))/100
  obtained <- hypergeomPFQ(m = 15, a = 4i, b = NULL, x = X)
  expected <- det(diag(3)-X)^(-4i)
  expect_equal(obtained, expected)
  #
  X <- toeplitz(c(3i,2,1))/100
  obtained <- hypergeomPFQ(m = 15, a = 3, b = NULL, x = X)
  expected <- complexplus::Det(diag(3)-X)^(-3)
  expect_equal(obtained, expected)
  #
  X <- toeplitz(c(2i,2,1))/100
  obtained <- hypergeomPFQ(m = 15, a = 4i, b = NULL, x = X)
  expected <- complexplus::Det(diag(3)-X)^(-4i)
  expect_equal(obtained, expected)
})


