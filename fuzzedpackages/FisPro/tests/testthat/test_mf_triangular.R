context("mf triangular")

test_that("mf triangular degrees", {
  mf <- new(mf_triangular, 0, 1, 2)

  expect_equal(mf$degree(0), 0)
  expect_equal(mf$degree(0.5), 0.5)
  expect_equal(mf$degree(1), 1)
  expect_equal(mf$degree(1.5), 0.5)
  expect_equal(mf$degree(2), 0)
})
