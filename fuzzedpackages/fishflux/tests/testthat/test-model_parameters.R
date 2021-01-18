# example

test_that("Simple corner cases", {
  needs_api()
  check_api()
  mp <- suppressMessages(suppressWarnings(
    fishflux::model_parameters(
      sp = "Scarus psittacus", family = "Scaridae", temp = 27)))
  expect_length(mp, 15)
  expect_equal(nrow(mp), 1)
  expect_true(!is.numeric(mp$species))
  expect_true(is.numeric(mp$t0))
  expect_true(is.numeric(mp$Linf))
  expect_true(is.numeric(mp$k))
  expect_true(is.numeric(mp$asp))
  expect_true(is.numeric(mp$troph))
  expect_true(is.numeric(mp$lwa_m))
  expect_true(is.numeric(mp$lwa_sd))
  expect_true(is.numeric(mp$lwb_m))
  expect_true(is.numeric(mp$lwb_sd))
  expect_true(is.numeric(mp$mdw_m))
  expect_true(is.numeric(mp$f0_m))
  expect_true(is.numeric(mp$f0_sd))
  expect_true(is.numeric(mp$alpha_m))
  expect_true(is.numeric(mp$alpha_sd))
})
