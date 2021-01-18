# Tests of methods
# Arseniy Khvorov
# Created 2020/01/03
# Last edit 2020/01/03

library(dplyr)

pop <- sim_reference(
  init_pop_size = 1e6L,
  vaccinations = generate_counts(1e6L, 304L, 0.55, 100, 50),
  cases_novac = generate_counts(1e6L, 304L, 0.12, 190, 35),
  ve = 0.48,
  lag = 0L,
  seed = 1L,
  deterministic = TRUE
)
pop_monthly <- pop %>%
  mutate(
    dates = generate_dates(timepoint, lubridate::ymd("2017/08/01"), "day"),
    month = lubridate::month(dates),
    year = lubridate::year(dates)
  ) %>%
  group_by(year, month) %>%
  summarise(
    vaccinations = sum(vaccinations),
    cases = sum(cases),
    ve = mean(ve)
  ) %>%
  ungroup()

test_that("method1 works", {
  m1 <- method1(
    attr(pop, "init_pop_size"),
    pop_monthly$vaccinations, pop_monthly$cases, pop_monthly$ve
  )
  with(m1, {
    expect_equal(vaccinations, pop_monthly$vaccinations)
    expect_equal(cases, pop_monthly$cases)
    expect_equal(ve, pop_monthly$ve)
    expect_equal(vc_lag, (pvac + lag(pvac, default = 0)) / 2)
    expect_equal(
      pops, as.integer(round(
        (lag(pops, default = attr(pop, "init_pop_size")) -
          lag(cases, default = 0L)) * (1 - vc_lag * ve)
      ), 0)
    )
    expect_equal(pflu, cases / pops)
    expect_equal(
      popn, lag(popn, default = attr(pop, "init_pop_size")) -
        lag(cases_novac, default = 0L)
    )
    expect_equal(cases_novac, as.integer(round(popn * pflu), 0))
    expect_equal(avert, cases_novac - cases)
  })
})

test_that("method3 works", {
  m3 <- method3(
    attr(pop, "init_pop_size"),
    pop_monthly$vaccinations, pop_monthly$cases, pop_monthly$ve
  )
  with(m3, {
    expect_equal(vaccinations, pop_monthly$vaccinations)
    expect_equal(cases, pop_monthly$cases)
    expect_equal(ve, pop_monthly$ve)
    expect_equal(
      pflu, cases / (
        lag(A, default = attr(pop, "init_pop_size")) + lag(C, default = 0L)
    ))
    expect_equal(
      pvac, vaccinations / (
        lag(A, default = attr(pop, "init_pop_size")) + lag(E, default = 0L)
      ))
    expect_equal(
      b, as.integer(round(
        pvac * lag(A, default = attr(pop, "init_pop_size"))
      ), 0)
    )
    expect_equal(
      A, as.integer(round(
        (1 - pflu) * lag(A, default = attr(pop, "init_pop_size")) - b
      ), 0)
    )
    expect_equal(
      C, as.integer(round(
        (1 - pflu) * lag(C, default = 0L) + b * (1 - ve)
      ), 0)
    )
    expect_equal(
      D, as.integer(round(
        lag(D, default = 0L) + b * ve
      ), 0)
    )
    expect_equal(
      E, as.integer(round(
        (1 - pvac) * lag(E, default = 0L) +
          lag(A, default = attr(pop, "init_pop_size")) * pflu
      ), 0)
    )
    expect_equal(
      `F`, as.integer(round(
        lag(`F`, default = 0L) +
          lag(C, default = 0L) * pflu + lag(E, default = 0L) * pvac
      ), 0)
    )
    expect_equal(
      cases_novac, as.integer(round(
        lag(popn, default = attr(pop, "init_pop_size")) * pflu
      ), 0)
    )
    expect_equal(
      popn, lag(popn, default = attr(pop, "init_pop_size")) - cases_novac
    )
    expect_equal(avert, cases_novac - cases)
  })
})
