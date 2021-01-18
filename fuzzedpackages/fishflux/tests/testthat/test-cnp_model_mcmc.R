context("cnp_model_mcmc")

library(fishflux)
model <- suppressWarnings(cnp_model_mcmc(TL = 10, param = list(Qc_m = 40, Qn_m = 10, Qp_m = 4, theta_m = 3)))

test_that("Simple corner cases", {
    expect_is(model, "list")
    expect_warning(cnp_model_mcmc(TL = 10, param = list(Qc_m = 40, Qn_m = 10, Qp_m = 4, theta_m = 3)))
    expect_error(cnp_model_mcmc(TL = 10))
    expect_error(cnp_model_mcmc(TL = 10, param = list(Qc_m = -40, Qn_m = 10, Qp_m = 4, theta_m = 3)))
    expect_equal(length(model), 2)
    expect_is(model[[1]], "stanfit")
    expect_is(model[[2]], "data.frame")
})
