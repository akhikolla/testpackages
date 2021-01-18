# example
mod <- suppressMessages(suppressWarnings(fishflux::cnp_model_mcmc(TL = 5:15, param = list(
  Qc_m = 40, Qn_m = 10, Qp_m = 4, Dn_sd = 0.05))))
plt <- fishflux::plot_cnp(mod = mod, y = c("Fp", "Gp", "Wp", "Ip"),
                   x = "tl", probs = c(0.5, 0.8))

test_that("Simple corner cases", {
  expect_true(ggplot2::is.ggplot(plt))
})
