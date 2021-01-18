context("kgaps_imt")

# Check that calling kgaps_imt() with vector arguments u and k gives
# the same results as calling kgaps_imt() repeatedly with scalar arguments

u <- stats::quantile(newlyn, probs = c(0.85, 0.90, 0.95))
k_vals <- 1:4
all_res <- kgaps_imt(newlyn, u, k_vals)
all_IMT <- all_res$imt
all_p <- all_res$p
all_theta <- all_res$theta

ind_IMT <- ind_p <- ind_theta <- all_IMT
for (i in 1:length(u)) {
  for (j in 1:length(k_vals)) {
    temp <- kgaps_imt(newlyn, u = u[i], k = k_vals[j])
    ind_IMT[i, j] <- temp$imt
    ind_p[i, j] <- temp$p
    ind_theta[i, j] <- temp$theta
  }
}

my_tol <- 1e-5

test_that("IMT values agree", {
  testthat::expect_equal(all_IMT, ind_IMT, tolerance = my_tol)
})
test_that("p-values agree", {
  testthat::expect_equal(all_p, ind_p, tolerance = my_tol)
})
test_that("MLEs of theta values agree", {
  testthat::expect_equal(all_theta, ind_theta, tolerance = my_tol)
})
