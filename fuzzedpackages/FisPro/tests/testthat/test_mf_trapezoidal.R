context("mf trapezoidal")

test_that("mf trapezoidal degrees", {
  mf <- new(mf_trapezoidal, 0, 1, 2, 3)

  expect_equal(mf$degree(0), 0)
  expect_equal(mf$degree(0.5), 0.5)
  expect_equal(mf$degree(1), 1)
  expect_equal(mf$degree(2), 1)
  expect_equal(mf$degree(2.5), 0.5)
  expect_equal(mf$degree(3), 0)
})
