context("Multivariate Gamma function")

test_that("Comparison with gamma", {
  expect_equal(mvgamma(-1.5,1), gamma(-1.5))
  expect_equal(mvgamma(-2.5,1), gamma(-2.5))
})

test_that("Complex value for p=1", {
  y <- 2
  o1 <- Mod(mvgamma(y*1i,1))
  o2 <- sqrt(pi/y/sinh(pi*y))
  expect_equal(o1, o2)
})
