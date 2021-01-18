context("GP: upper")

# We check that if we specify a finite upper endpoint for the GP distribution
# then the values posterior samples is consistent with this

# The primary motivation for this comes from the threshr package.  If we
# Box-Cox(lambda) transform the data prior to a GP analysis then we know that
# the upper end point of the transformed distribution is < -1/lambda,
# and should therefore impose this constraint in the posterior sampling.
# We use the "flatflat" prior distribution that is used in threshr for this
# purpose.

# Two potential problems:
#
# 1. The mode of the posterior is likely to be negative, perhaps strongly.
#    If it is strongly negative then it is advisable to use trans = "BC".
#
# 2. It is possible that the (unconstrained) mode of the posterior lies outside
#    of the support of the posterior once the the constraint has been imposed.
#    In this event we will get a warning (at best) and if results are produced
#    they are not to be trusted.  There may be similar problems if the
#    unconstrained mode lies inside the constrained support but close to its
#    boundary.

# Set a prior: flat for GP parameters, Haldane for P(exceedance)
prior_args <- list(prior = "flatflat", bin_prior = "haldane",
                   h_prior = list(min_xi = -Inf))

set.seed(5092019)

lambda_vec <- c(-0.5, -0.1)
n <- 1000
# Sample from a half-normal distribution
x <- abs(rnorm(10000))
# Threshold
thresh <- stats::quantile(x, probs = 0.7)

for (i in 1:length(lambda_vec)) {
  # Box-Cox transform: data and threshold
  lambda <- lambda_vec[i]
  y <- (x ^ lambda - 1) / lambda
  u <- (thresh ^ lambda - 1) / lambda
#  y <- threshr:::bc(x, lambda)
#  u <- threshr:::bc(thresh, lambda)
  # Set prior
  fp <- set_prior(prior = "flatflat", model = "gp", upper = - 1 / lambda - u)
  res <- rpost(n = 1000, model = "gp", prior = fp, thresh = u, data = y,
               trans = "BC")
  sigma <- res$sim_vals[, "sigma[u]"]
  xi <- res$sim_vals[, "xi"]
  xf <- u - sigma / xi

  testthat::test_that(paste0("Half-normal, Box-Cox lambda = ", lambda), {
    testthat::expect_true(all(xf <= -1 / lambda))
  })
}
