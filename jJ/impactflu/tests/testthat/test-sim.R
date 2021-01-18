# Tests of simulation functions
# Arseniy Khvorov
# Created 2019/12/24
# Last edit 2019/12/24

library(testthat)

test_that("count generation works", {
  set.seed(1)
  test_counts <- generate_counts(
    init_pop_size = 1e6, n_timepoints = 304,
    overall_prop = 0.55, mean = 100, sd = 50
  )
  expect_equal(length(test_counts), 304)
  expect_equal(sum(test_counts) / 1e6, 0.55, tol = 0.01)
})

test_that("date generation works", {
  timepoints <- 1:5
  start <- lubridate::ymd("2000-01-01")
  expect_equal(
    generate_dates(timepoints, start, "day"),
    c(
      lubridate::ymd("2000-01-01"), lubridate::ymd("2000-01-02"),
      lubridate::ymd("2000-01-03"), lubridate::ymd("2000-01-04"),
      lubridate::ymd("2000-01-05")
    )
  )
  expect_equal(
    generate_dates(timepoints, start, "month"),
    c(
      lubridate::ymd("2000-01-01"), lubridate::ymd("2000-02-01"),
      lubridate::ymd("2000-03-01"), lubridate::ymd("2000-04-01"),
      lubridate::ymd("2000-05-01")
    )
  )
  expect_equal(
    generate_dates(timepoints, start, "year"),
    c(
      lubridate::ymd("2000-01-01"), lubridate::ymd("2001-01-01"),
      lubridate::ymd("2002-01-01"), lubridate::ymd("2003-01-01"),
      lubridate::ymd("2004-01-01")
    )
  )
  expect_error(
    generate_dates(timepoints, start, "wrong"),
    "unrecognised unit 'wrong'"
  )
})

test_that("simulation works with no lag", {
  pop <- sim_reference(
    init_pop_size = 1e6L,
    vaccinations = generate_counts(1e6L, 304L, 0.55, 100, 50),
    cases_novac = generate_counts(1e6L, 304L, 0.12, 190, 35),
    ve = 0.48,
    lag = 0L,
    seed = 1L,
    deterministic = TRUE
  )
  expect_named(
    pop, c(
      "timepoint", "vaccinations", "cases_novac", "ve", "pflu",
      "popn", "pvac", "b", "A", "B", "C", "D", "E", "F", "cases", "avert"
  ))
  expect_equal(attr(pop, "seed"), 1L)
  expect_equal(attr(pop, "init_pop_size"), 1e6L)
  expect_equal(attr(pop, "lag"), 0L)
  with(pop, {
    expect_equal(timepoint, 1L:304L)
    expect_equal(pflu, cases_novac / dplyr::lag(popn, default = 1e6L))
    expect_equal(
      cases,
      as.integer(round(pflu * dplyr::lag(A, default = 1e6L), 0)) +
      as.integer(round(pflu * dplyr::lag(C, default = 0L), 0))
    )
    expect_equal(popn, dplyr::lag(popn, default = 1e6L) - cases_novac)
    expect_equal(avert, cases_novac - cases)
    expect_equal(
      pvac, vaccinations /
        (dplyr::lag(A, default = 1e6L) + dplyr::lag(E, default = 0L))
    )
    expect_equal(b, as.integer(round(pvac * dplyr::lag(A, default = 1e6L), 0)))
    expect_equal(
      A, dplyr::lag(A, default = 1e6L) -
        as.integer(round(dplyr::lag(A, default = 1e6L) * pflu), 0) - b
    )
    expect_equal(B, rep(0L, 304))
    expect_equal(
      C, dplyr::lag(C, default = 0L) -
        as.integer(round(dplyr::lag(C, default = 0L) * pflu, 0)) +
        as.integer(round(b * (1 - ve), 0))
    )
    expect_equal(
      D, dplyr::lag(D, default = 0L) + b - as.integer(round(b * (1 - ve), 0))
    )
    expect_equal(
      E, dplyr::lag(E, default = 0L) +
        as.integer(round(dplyr::lag(A, default = 1e6L) * pflu, 0)) -
        as.integer(round(dplyr::lag(E, default = 0L) * pvac, 0))
    )
    expect_equal(
      `F`, dplyr::lag(`F`, default = 0L) +
        as.integer(round(dplyr::lag(C, default = 0L) * pflu, 0)) +
        as.integer(round(dplyr::lag(E, default = 0L) * pvac, 0))
    )
  })
})

test_that("simulation works with lag", {
  pop <- sim_reference(
    init_pop_size = 1e6L,
    vaccinations = generate_counts(1e6L, 304L, 0.55, 100, 50),
    cases_novac = rep(0L, 304),
    ve = 0.48,
    lag = 1L,
    seed = 1L,
    deterministic = TRUE
  )
  expect_equal(pop$B, pop$b)
  pop <- sim_reference(
    init_pop_size = 1e6L,
    vaccinations = generate_counts(1e6L, 304L, 0.55, 100, 50),
    cases_novac = rep(0L, 304),
    ve = 0.48,
    lag = 2L,
    seed = 1L,
    deterministic = TRUE
  )
  expect_equal(pop$B, pop$b + dplyr::lag(pop$b, default = 0L))
})

test_that("random simulation works", {
  pop <- sim_reference(
    init_pop_size = 1e6L,
    vaccinations = generate_counts(1e6L, 304L, 0.55, 100, 50),
    cases_novac = generate_counts(1e6L, 304L, 0.12, 190, 35),
    ve = 0.48,
    lag = 0L,
    seed = 1L,
    deterministic = FALSE
  )
  sum1 <- sum(pop$avert)
  pop <- sim_reference(
    init_pop_size = 1e6L,
    vaccinations = generate_counts(1e6L, 304L, 0.55, 100, 50),
    cases_novac = generate_counts(1e6L, 304L, 0.12, 190, 35),
    ve = 0.48,
    lag = 0L,
    seed = 1L,
    deterministic = FALSE
  )
  sum2 <- sum(pop$avert)
  expect_equal(sum1, sum2)
  pop <- sim_reference(
    init_pop_size = 1e6L,
    vaccinations = generate_counts(1e6L, 304L, 0.55, 100, 50),
    cases_novac = generate_counts(1e6L, 304L, 0.12, 190, 35),
    ve = 0.48,
    lag = 0L,
    seed = 2L,
    deterministic = FALSE
  )
  sum3 <- sum(pop$avert)
  expect_true(sum1 != sum3)
})
