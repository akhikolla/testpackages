context("Errors")

test_that("errors on m", {
  expect_error(hypergeomPFQ(m = -1, a = c(1,2), b = c(3,4), x = c(5,6)))
  expect_error(hypergeomPFQ(m = c(1,2), a = c(1,2), b = c(3,4), x = c(5,6)))
})

test_that("errors on alpha", {
  expect_error(hypergeomPFQ(m = 5, a = c(1,2), b = c(3,4), x = c(5,6),
                            alpha = "a"))
  expect_error(hypergeomPFQ(m = 5, a = c(1,2), b = c(3,4), x = c(5,6),
                            alpha = -1))
  expect_error(hypergeomPFQ(m = 5, a = c(1,2), b = c(3,4), x = c(5,6),
                            alpha = c(1,2)))
})

test_that("errors on x", {
  expect_error(hypergeomPFQ(m = 5, a = c(1,2), b = c(3,4),
                            x = cbind(c(1,1),c(2,3),c(4,5))))
  expect_error(hypergeomPFQ(m = 5, a = c(1,2), b = c(3,4), x = "a"))
  expect_error(hypergeomPFQ(m = 5, a = c(1,2), b = c(3,4), x = list(5,6)))
  expect_error(hypergeomPFQ(m = 5, a = c(1,2), b = c(3,4), x = NULL))
  expect_error(hypergeomPFQ(m = 5, a = c(1,2), b = c(3,4), x = numeric(0L)))
})

test_that("errors lmvgamma", {
  expect_error(lmvgamma(-1+2i, 3))
  expect_error(lmvgamma(1+2i, c(2,3)))
  expect_error(lmvgamma(c(3,4), 3))
})

test_that("errors mvgamma", {
  expect_error(mvgamma(c(-1+2i,3), 3))
  expect_error(mvgamma(1+2i, c(2,3)))
  expect_error(mvgamma(-2, 3))
})

test_that("errors mvbeta", {
  expect_error(mvbeta(a = c(1,2), b = 3, p = 4))
  expect_error(mvbeta(a = 2, b = 3, p = c(3,4)))
})

test_that("errors IncBeta", {
  expect_error(IncBeta(m = 5, a = 0.5, b = 2, x = c(1,2,3)))
  expect_error(IncBeta(m = 5, a = c(1,5), b = 2, x = c(1,2,3)))
  expect_error(IncBeta(m = 5, a = 5, b = 2, x = diag(c(0.5,-0.5))))
  expect_error(IncBeta(m = 5, a = 5, b = 2, x = diag(c(1.1,0.5))))
})

test_that("errors IncGamma", {
  expect_error(IncGamma(m = 5, a = 1, x = toeplitz(3:1)))
})

test_that("errors BesselA", {
  expect_error(BesselA(m = 5, x = diag(c(1,2+2i)), nu = -1+1i))
})
