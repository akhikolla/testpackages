# example
m <- suppressMessages(suppressWarnings(fishflux::metabolism(family = "Pomacentridae", temp = 27, troph_m = 2)))

test_that("Simple corner cases", {
  expect_gt(min(m), 0)
  expect_length(m, 6)
  expect_equal(nrow(m), 1)
  expect_true(is.numeric(m$f0_m))
  expect_true(is.numeric(m$f0_sd))
  expect_true(is.numeric(m$alpha_m))
  expect_true(is.numeric(m$alpha_sd))
  expect_true(is.numeric(m$b0_m))
  expect_true(is.numeric(m$b0_sd))
})
