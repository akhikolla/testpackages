#example
mod <- suppressWarnings(fishflux::cnp_model_mcmc(TL = 5:10,
	                                             param = list(Qc_m = 40, Qn_m = 10, Qp_m = 4,
                                                              Dc_sd = 0.1, Dn_sd = 0.05, Dp_sd = 0.05)))
lim <- fishflux::limitation(mod)

test_that("Simple corner cases", {
  expect_length(lim, 3)
  expect_equal(nrow(lim), 18)
  expect_gte(min(lim$prop_lim), 0)
  expect_gte(min(lim$tl), 0)
  expect_true(is.character(lim$nutrient))
})
