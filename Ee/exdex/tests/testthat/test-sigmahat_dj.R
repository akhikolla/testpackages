context("sigmahat_dj")

# Check that the various ways of calculating a summation involved in the
# estimation of sigma^2_dj on page 2319 of Berghaus and Bucher (2018)
# give the same result.

# Calculate the contributions to the sum from each (disjoint) block
# Each row corresponds to a different block
# Each column corresponds to a different way to calculate the summation
# The 1st column is used in estimation of the sampling variance of thetahat

# The permitted range of b for these data is 15 - 196
my_tol <- 1e-5
b_vec <- c(15, 100, 196)
for (i in 1:length(b_vec)) {
  sigma_mat <- spm_sigmahat_dj(data = newlyn, b = b_vec[i], check = TRUE)
  test_that("Usum2 = Usum", {
    testthat::expect_equal(sigma_mat[, 2], sigma_mat[, 1], tolerance = my_tol)
  })
  test_that("Usum3 = Usum", {
    testthat::expect_equal(sigma_mat[, 3], sigma_mat[, 1], tolerance = my_tol)
  })
  test_that("Usum4 = Usum", {
    testthat::expect_equal(sigma_mat[, 4], sigma_mat[, 1], tolerance = my_tol)
  })
}

# Check that the R function ests_sigmahat_dj() gives the same results as the
# Rcpp function cpp_sigma2hat_dj

