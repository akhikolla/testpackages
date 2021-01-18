# example

test_that("Simple corner cases", {
  needs_api()
  check_api()
  l <- fishflux::trophic_level("Zebrasoma scopas")
  expect_gt(min(l$trophic_level), 0)
  expect_length(l, 3)
  expect_equal(nrow(l), 1)
  expect_true(is.numeric(l$trophic_level))
  expect_true(!is.numeric(l$species))
  expect_true(!is.numeric(l$level))
})
