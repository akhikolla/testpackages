context("mf trapezoidal sup")

test_that("mf trapezoidal sup degrees", {
  mf <- new(mf_trapezoidal_sup, 0, 1)

  expect_equal(mf$degree(0), 0)
  expect_equal(mf$degree(0.5), 0.5)
  expect_equal(mf$degree(1), 1)
})
