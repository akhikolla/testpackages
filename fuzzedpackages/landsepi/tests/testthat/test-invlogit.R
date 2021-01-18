context("test-invlogit")

test_that("invlogit produce a 'NaN' when misused", {
  expect_equal(invlogit(NaN), NaN)
  expect_equal(invlogit(Inf), NaN)
  expect_equal(invlogit(1000), NaN)
})

test_that("invlogit works", {
  expect_equal(invlogit(0), 0.5)
  expect_equal(invlogit(100), 1)
  expect_equal(invlogit(0.5), 0.6224593, tolerance = 1e-7)
  expect_equal(invlogit(-Inf), 0)
})
