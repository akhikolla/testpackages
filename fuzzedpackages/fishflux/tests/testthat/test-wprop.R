
test_that("Simple corner cases", {
  expect_gt(min(fishflux::wprop(family="Scaridae")), 0)
  expect_length(fishflux::wprop(family="Scaridae"), 2)
  expect_equal(nrow(fishflux::wprop(family="Scaridae")), 1)
  expect_s3_class(fishflux::wprop(family="Scaridae"), "data.frame")
})
