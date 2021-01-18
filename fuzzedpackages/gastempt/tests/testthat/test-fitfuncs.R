context("Test fit functions")
tolerance = 5.e-5

test_that("Linexp Functions can be called with par or explicit parameters",{
  v0 = 400
  tempt = 60
  kappa = 2
  t = c(0, 10, 300)

  # linexp
  pars = c(v0 = v0, tempt = tempt, kappa = kappa)
  r = linexp(t, v0, tempt, kappa)
  rp = linexp(t, pars = pars)
  rpd = linexp(t, pars = data.frame(t(pars)))
  expect_equal(r, rp)
  expect_equal(r, rpd)
  expect_equal(length(r), 3)
  expect_error(linexp(t, v0, tempt, kappa, pars), "Either")
  expect_error(linexp(t, tempt = tempt, kappa = kappa, pars = pars), "Either")

  # linexp_slope
  pars = c(v0 = v0, tempt = tempt, kappa = kappa)
  r = linexp_slope(t, v0, tempt, kappa)
  rp = linexp_slope(t, pars = pars)
  rpd = linexp_slope(t, pars = data.frame(t(pars)))
  expect_equal(r, rp)
  expect_equal(r, rpd)
  expect_equal(length(r), 3)
  expect_error(linexp_slope(t, v0, tempt, kappa, pars), "Either")
  expect_error(linexp_slope(t, tempt = tempt, kappa = kappa, pars = pars), "Either")

  # linexp_auc
  pars = c(v0 = v0, tempt = tempt, kappa = kappa)
  r = linexp_auc(v0, tempt, kappa)
  rp = linexp_auc(pars = pars)
  rpd = linexp_auc(pars = data.frame(t(pars)))
  expect_equal(r, rp)
  expect_equal(r, rpd)
  expect_error(linexp_auc(v0, tempt, kappa, pars), "Either")
  expect_error(linexp_auc(tempt = tempt, kappa = kappa, pars = pars), "Either")

  # linexp_log
  logtempt = log(tempt)
  logkappa = log(kappa)
  pars = c(v0 = v0, logtempt = logtempt, logkappa = logkappa)
  r = linexp_log(t, v0, logtempt, logkappa)
  rp = linexp_log(t, pars = pars)
  rpd = linexp_log(t, pars = data.frame(t(pars)))
  expect_equal(r, rp)
  expect_equal(r, rpd)
  expect_equal(length(r), 3)
  expect_error(linexp_log(t, v0, logtempt, logkappa, pars), "Either")
  expect_error(linexp_log(t, logtempt = logtempt, logkappa = logkappa,
                          pars = pars), "Either")
})

test_that("Powexp Functions can be called with pars or explicit parameters",{
  v0 = 400
  tempt = 60
  beta = 2
  t = c(0, 10, 300)

  # powexp
  pars = c(v0 = v0, tempt = tempt, beta = beta)
  r = powexp(t, v0, tempt, beta)
  rp = powexp(t, pars = pars)
  rpd = powexp(t, pars = data.frame(t(pars)))
  expect_equal(r, rp)
  expect_equal(r, rpd)
  expect_equal(length(r), 3)
  expect_error(powexp(t, v0, tempt, beta, pars), "Either")
  expect_error(powexp(t, tempt = tempt, beta = beta, pars = pars), "Either")

  # powexp_slope
  pars = c(v0 = v0, tempt = tempt, beta = beta)
  r = powexp_slope(t, v0, tempt, beta)
  rp = powexp_slope(t, pars = pars)
  expect_equal(r, rp)
  expect_equal(length(r), 3)
  expect_error(powexp_slope(t, v0, tempt, beta, pars), "Either")
  expect_error(powexp_slope(t, tempt = tempt, beta = beta, pars = pars), "Either")


  # There is no powexp_auc

  # powexp_log
  logtempt = log(tempt)
  logbeta = log(beta)
  pars = c(v0 = v0, logtempt = logtempt, logbeta = logbeta)
  r = powexp_log(t, v0, logtempt, logbeta)
  rp = powexp_log(t, pars = pars)
  rpd = powexp_log(t, pars = data.frame(t(pars)))
  expect_equal(r, rp)
  expect_equal(r, rpd)
  expect_equal(length(r), 3)
  expect_error(powexp_log(t, v0, logtempt, logbeta, pars), "Either")
  expect_error(powexp_log(t, logtempt = logtempt, logbeta = logbeta,
                          pars = pars), "Either")
})


test_that("Functions at t=0 must return initial volume",{
  v0 = 400
  tempt = 60
  kappa = 2
  t = c(0, 10, 300)
  # linexp
  r = linexp(t, v0, tempt, kappa)
  expect_equal(r[1], v0)
  expect_gt(r[2], r[1])
  expect_lt(r[3], r[1])

  rlog = linexp_log(t, v0, log(tempt), log(kappa))
  expect_equal(rlog[1], v0)
  expect_gt(rlog[2], rlog[1])
  expect_lt(rlog[3], rlog[1])

  expect_equal(r, rlog)

  # Power exponential
  beta = 2
  r = powexp(t, v0, tempt, beta)
  expect_equal(r[1], v0)
  expect_lt(r[2], r[1])
  expect_lt(r[3], r[1])

  rlog = powexp_log(t, v0, log(tempt), log(beta))
  expect_equal(rlog[1], v0)
  expect_lt(rlog[2], r[1])
  expect_lt(rlog[3], r[1])

  expect_equal(r, rlog)

})

test_that("Limiting cases of slopes are correct",{
  tempt = 60
  kappa = 1
  # linexp
  t = c(0, 10, 300)
  r = linexp_slope(t, tempt = tempt, kappa = 1)
  expect_equal(r[1], 0)
  expect_lt(r[2], 0)
  expect_lt(r[3], 0)
  r = linexp_slope(0, tempt = tempt, kappa = 0)*tempt
  expect_equal(r, -1)

  # special values
  expect_equal(powexp_slope(0, 100, 100,  1), -1)
  expect_equal(powexp_slope(0, 1, 100,  1), -0.01)
  expect_equal(powexp_slope(0, 100, 100,  2), 0)
  expect_equal(powexp_slope(0, 100, 100,  3), 0)

    r = powexp_slope(1e-5, tempt = tempt, beta = 1)*tempt
  expect_equal(r, -1, tolerance = tolerance)

})


