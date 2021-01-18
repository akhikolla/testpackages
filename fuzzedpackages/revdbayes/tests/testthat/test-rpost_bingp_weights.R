context("binGP rpost: weights vs no weights")

# We check that the values simulated using rpost() are identical when weights
# is not supplied and when weights is a vector of ones.

# Set a tolerance for the comparison of the simulated values
my_tol <- 1e-5

u <- stats::quantile(gom, probs = 0.65)
fp <- set_prior(prior = "flat", model = "gp", min_xi = -1)
bp <- set_bin_prior(prior = "jeffreys")
set.seed(14092019)
res1 <- rpost(n = 1000, model = "bingp", prior = fp, thresh = u, data = gom,
              bin_prior = bp)
set.seed(14092019)
res2 <- rpost(n = 1000, model = "bingp", prior = fp, thresh = u, data = gom,
              bin_prior = bp, weights = rep(1, length(gom)))

testthat::test_that("binGP rpost: weights vs no weights", {
  testthat::expect_equal(res1$sim_vals, res2$sim_vals, tolerance = my_tol)
})
