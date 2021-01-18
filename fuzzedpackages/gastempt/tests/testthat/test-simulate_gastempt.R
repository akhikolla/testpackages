context("Test simulation of data")

library(testthat)
library(assertthat)
suppressWarnings(RNGversion("3.5.0"))

test_that("Default call of simulate_gastempt must return plausible values",{
  d = simulate_gastempt(seed = 4711)
  expect_equal(names(d), c("record", "data", "stan_data"))
  expect_match(comment(d$data), "linexp")
  data = d$data
  record = d$record
  expect_equal(nrow(record), 10)
  expect_gte(nrow(data), 100)
  expect_true(all(data$vol >= 0))
  expect_lt(median(data$vol[data$minute > max(data$minute)*0.7]), 200)
})

test_that("When max_minute is explicitly given, it must be used",{
  d = simulate_gastempt(max_minute = 88, seed = 4711)
  expect_equal(max(d$data$minute), 60)
})


test_that("Noise = 0 must issue a warning",{
  expect_warning(d <- simulate_gastempt(noise = 0, seed = 4711),
                 "might fail")
  expect_equal(names(d), c("record", "data", "stan_data"))
  data = d$data
  record = d$record
  expect_equal(nrow(record), 10)
  expect_gte(nrow(data), 100)
  # Without noise, first
  expect_equal(data$vol[data$minute == 0], record$v0)
})

test_that("student_t_df < 0 must issue a warning",{
  set.seed(4711)
  expect_warning(d <- simulate_gastempt(student_t_df = 1), "freedom")
})

test_that("When data are missing, records must be shortened",{
  d = simulate_gastempt(missing = 0.5, seed = 4711)
  data = d$data
  record = d$record
  expect_equal(nrow(record), 10)
  expect_gte(nrow(data), 50)
})

test_that("Warning must be issued for invalid values of missing",{
  expect_warning(simulate_gastempt(missing = 1, seed = 4711), "0.5")
  expect_warning(simulate_gastempt(missing = -1, seed = 4711), "0")
})

test_that("Should stop when linexp input parameters are invalid",{
  expect_error(simulate_gastempt(kappa_mean = -1),"greater")
  expect_error(simulate_gastempt(kappa_mean = 1, kappa_std = 1),"greater")

  expect_error(simulate_gastempt(tempt_mean = -1),"greater")
  expect_error(simulate_gastempt(tempt_mean = 100, tempt_std = 100),"greater")

  expect_error(simulate_gastempt(v0_mean = -1),"greater")
  expect_error(simulate_gastempt(v0_mean = 400, v0_std = 210),"greater")
})


# powexp
test_that(
  "Default call of simulate_gastempt for powexp must return plausible values",{
  set.seed(4711)
  d = simulate_gastempt(model = powexp, seed = 4711)
  expect_equal(names(d), c("record", "data", "stan_data"))
  data = d$data
  record = d$record
  expect_equal(nrow(record), 10)
  expect_gte(nrow(data), 80)
  expect_true(all(data$vol > 0))
})

test_that("Noise = 0 must issue a warning and powexp should not overshoot",{
  expect_warning(d <- simulate_gastempt(
    n_records = 4, v0_mean = 400, v0_std = 0,
    beta_mean = 4,
    noise = 0, model = powexp, seed = 4711), "might fail")
  data = d$data
  record = d$record
  expect_equal(names(d), c("record", "data", "stan_data"))
  expect_true(all(data$vol >= 0))
  expect_true(all(data$vol <= 400))
  # Without noise, first
  expect_equal(data$vol[data$minute == 0], record$v0)
})


test_that("Should stop when powexp input parameters are invalid",{
  set.seed(4711)
  expect_error(simulate_gastempt(beta_mean = -1, model = powexp),"greater")
  expect_error(simulate_gastempt(beta_mean = 1, beta_std = 1,
                                 model = powexp),"greater")

  expect_error(simulate_gastempt(tempt_mean = -1, model = powexp), "greater")
  expect_error(simulate_gastempt(tempt_mean = 100, tempt_std = 100,
                                 model = powexp),"greater")

  expect_error(simulate_gastempt(v0_mean = -1, model = powexp),"greater")
  expect_error(simulate_gastempt(v0_mean = 400, model = powexp, v0_std = 210),"greater")
})
