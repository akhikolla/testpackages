context("test-logit")

test_that("logit produce a 'NaN' when misused", {
  expect_equal(logit(NaN), NaN)
  expect_equal(logit(Inf), NaN)
  expect_equal(logit(-Inf), NaN)

  # Trick to remove the warnings from the test report
  expect_warning(res1 <- logit(2), "NaN")
  expect_warning(res2 <- logit(-1), "NaN")

  expect_equal(logit(res1), NaN)
  expect_equal(logit(res2), NaN)
})

test_that("logit works", {
  expect_equal(logit(0), -Inf)
  expect_equal(logit(1), Inf)
  expect_equal(logit(0.5), 0)
  expect_equal(logit(0.1), -2.197225, tolerance = 1e-6)
})
