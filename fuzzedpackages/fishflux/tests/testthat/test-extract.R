library(fishflux)
model <- suppressWarnings(cnp_model_mcmc(TL = 10, param = list(Qc_m = 40, Qn_m = 10, Qp_m = 4, theta_m = 3)))

test_that("simple corner cases", {
  expect_is(extract(model, "Ip"), "data.frame")
  expect_gt(min(extract(model, c("Ip", "Gc", "Fn"))), 0)
  expect_error(extract(model, "name"))
  expect_equal(length(extract(model, "Ip")), 8)
})
