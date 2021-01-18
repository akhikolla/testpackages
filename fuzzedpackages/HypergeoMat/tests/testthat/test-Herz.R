context("Herz's relations")

test_that("Herz's relation for 1F1", {
  X <- toeplitz(c(3,2,1))/100
  o1 <- hypergeomPFQ(m=15, a = 2, b = 3, x = X)
  o2 <- exp(sum(diag(X))) * hypergeomPFQ(m=15, a = 3-2, b = 3, x = -X)
  expect_equal(o1, o2)
  #
  X <- toeplitz(c(3i,2,1))/100
  o1 <- hypergeomPFQ(m=15, a = 2i, b = 3, x = X)
  o2 <- exp(sum(diag(X))) * hypergeomPFQ(m=15, a = 3-2i, b = 3, x = -X)
  expect_equal(o1, o2)
})

test_that("Herz's relation for 2F1", {
  X <- toeplitz(c(3,2,1))/100
  o1 <- hypergeomPFQ(m=15, a = c(1,2i), b = 3, x = X)
  o2 <- det(diag(3)-X)^(-2i) *
    hypergeomPFQ(15, a = c(3-1,2i), b = 3, -X%*%solve(diag(3)-X))
  expect_equal(o1, o2)
  #
  X <- toeplitz(c(3,2,1))/100
  o1 <- hypergeomPFQ(m=15, a = c(1,2), b = 3, x = X)
  o2 <- det(diag(3)-X)^(-2) *
    hypergeomPFQ(15, a = c(3-1,2), b = 3, -X%*%solve(diag(3)-X))
  expect_equal(o1, o2)
  #
  X <- toeplitz(c(3i,2,1))/100
  o1 <- hypergeomPFQ(m=15, a = c(1,2i), b = 3, x = X)
  o2 <- complexplus::Det(diag(3)-X)^(-2i) *
    hypergeomPFQ(15, a = c(3-1,2i), b = 3, -X%*%solve(diag(3)-X))
  expect_equal(o1, o2)
  #
  X <- toeplitz(c(3i,2,1))/100
  o1 <- hypergeomPFQ(m=15, a = c(1,2), b = 3, x = X)
  o2 <- complexplus::Det(diag(3)-X)^(-2) *
    hypergeomPFQ(15, a = c(3-1,2), b = 3, -X%*%solve(diag(3)-X))
  expect_equal(o1, o2)
})
