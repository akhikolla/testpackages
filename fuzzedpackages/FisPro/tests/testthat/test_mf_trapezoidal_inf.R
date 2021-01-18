context("mf trapezoidal inf")

test_that("mf trapezoidal inf degrees", {
  mf <- new(mf_trapezoidal_inf, 0, 1)

  expect_equal(mf$degree(0), 1)
  expect_equal(mf$degree(0.5), 0.5)
  expect_equal(mf$degree(1), 0)
})
