context("Test parameter extraction function")

tolerance = 1.e-5
suppressWarnings(RNGversion("3.5.0"))

checklinexp = function(tempt, kappa){
  # Non-transformed
  # Without v0
  tt0 = t50(c(tempt = tempt, kappa = kappa))
  expect_equal(attr(tt0,"fun"), linexp)
  v50 = as.numeric(attr(tt0,"fun")(tt0, tempt = tempt, kappa = kappa))
  expect_equal(v50, 0.5,  tolerance = tolerance)
  expect_is(attr(tt0,"slope"), "numeric")
  expect_is(attr(tt0,"auc"), "numeric")

  # With v0
  v0 = 400
  tt = t50(c(v0 = v0, tempt = tempt, kappa = kappa))
  expect_equal(attr(tt,"fun"), linexp)
  v50 = as.numeric(attr(tt,"fun")(tt, v0 = v0, tempt = tempt,
                                  kappa = kappa))
  expect_equal(v50, v0/2,  tolerance = tolerance)
  expect_equal(tt0, tt, check.attributes = FALSE, tolerance = tolerance )

  expect_is(attr(tt,"slope"), "numeric")
  slope = attr(tt,"slope")
  # Compute slope directly
  slope_0 = 10*diff(linexp(c(tt + 0.05, tt - 0.05), v0 = v0,
                           tempt = tempt, kappa = kappa))
  # This is very approximate only
  expect_equal(slope_0, slope, tolerance = 1e-3)

  auc = attr(tt,"auc")
  auc_0 = integrate(linexp, 0, Inf, v0, tempt, kappa)$value
  expect_equal(auc, auc_0, tolerance = 1e-3)

if (FALSE) {
  x = 0:400
  plot(x, linexp(x, v0, tempt, kappa), type ="l", ylim=c(0, 400))
  a = as.numeric(v50 + slope_0*tt)
  abline(a = a, b = -slope_0)
  abline(a = a, b = slope)
  abline(a = 200, b = 0)
}

  # Transformed
  ttt = t50(c(v0 = v0, logtempt = log(tempt), logkappa = log(kappa)))
  expect_equal(attr(ttt,"fun"), linexp_log)
  v50 = as.numeric(attr(ttt,"fun")(ttt, logtempt = log(tempt),
                                  logkappa = log(kappa)))
  expect_equal(v50, 0.5, tolerance = tolerance)
  expect_is(attr(ttt,"slope"), "numeric")

  expect_equal(tt0, ttt, check.attributes = FALSE, tolerance = tolerance )
}

test_that("linexp t50 returns correct value",{
  checklinexp(120, 2)
  checklinexp(120, 5)
  checklinexp(20, .1)
  checklinexp(.5, .5)
})




checkpowexp = function(tempt, beta){
  # Non-transformed
  # Without v0
  tt0 = t50(c(tempt = tempt, beta = beta))
  expect_equal(attr(tt0,"fun"), powexp)
  v50 = as.numeric(attr(tt0,"fun")(tt0, tempt = tempt, beta = beta))
  expect_equal(v50, 0.5,  tolerance = tolerance)
  expect_is(attr(tt0,"slope"), "numeric")

  # With v0
  v0 = 400
  tt = t50(c(v0 = v0, tempt = tempt, beta = beta))
  expect_equal(attr(tt,"fun"), powexp)
  v50 = as.numeric(attr(tt,"fun")(tt, v0 = v0, tempt = tempt, beta = beta))
  expect_equal(v50, v0/2,  tolerance = tolerance)
  expect_is(attr(tt,"slope"), "numeric")

  expect_equal(tt0, tt, check.attributes = FALSE, tolerance = tolerance )

  # Transformed
  ttt = t50(c(logtempt = log(tempt), logbeta = log(beta)))
  expect_equal(attr(ttt,"fun"), powexp_log)
  v50 = as.numeric(attr(ttt,"fun")(ttt, logtempt = log(tempt),
                                  logbeta = log(beta)))
  expect_equal(v50, 0.5, tolerance = tolerance)
  expect_is(attr(ttt,"slope"), "numeric")

  expect_equal(tt0, ttt, check.attributes = FALSE, tolerance = tolerance )
}

test_that("powexp t50 returns correct value",{
  checkpowexp(120, 2)
  checkpowexp(120, 5)
  checkpowexp(20, .1)
  checkpowexp(.5, .5)
})

test_that("t50 for data frame returns a data frame with column t50",{
  n = 3
  set.seed(4711)
  x = data.frame(
    record = letters[1:n],
    v0 = rnorm(n, 400, 10),
    tempt = rnorm(n, 100, 20),
    kappa = rnorm(n, 0.7, 0.1)
  )
  ret = t50(x)
  expect_is(ret, "data.frame")
  expect_equal(nrow(ret), 3)
  expect_false(any(sapply(ret, function(x) any(is.na(x)))))
  expect_equal(names(ret), c(names(x), c("t50", "slope_t50", "auc")))
})


