context("confint.spm")

my_tol <- 1e-5

# =================================== spm ====================================

res <- spm(newlyn, 100)

ci1 <- confint(res, type = "cholesky")
ci2 <- confint(res, type = "spectral")
test_that("sliding: cholesky and spectral are identical", {
  testthat::expect_identical(ci1$cis, ci2$cis)
})

# Check that the estimates of theta in res and returned from
# chandwich::adjust_loglik()

test_that("estimates of theta agree, sliding", {
  testthat::expect_equal(res$theta_sl, ci1$theta, tolerance = my_tol)
})

ci1 <- confint(res, maxima = "disjoint", interval_type = "both",
               type = "cholesky")
ci2 <- confint(res, maxima = "disjoint", interval_type = "both",
               type = "spectral")
test_that("disjoint: cholesky and spectral are identical", {
  testthat::expect_identical(ci1$cis, ci2$cis)
})

# Check estimates of theta

test_that("estimates of theta agree, disjoint", {
  testthat::expect_equal(res$theta_dj, ci1$theta, tolerance = my_tol)
})

ci3 <- confint(res, maxima = "disjoint", type = "cholesky",
               interval_type = "both", conf_scale = "log")

which_rows <- c("N2015lik", "BB2018lik")
test_that("spm lik intervals don't depend on conf_scale", {
  testthat::expect_identical(ci1$cis[which_rows, ], ci3$cis[which_rows, ])
})

# ============================= plot.confint.spm =============================

context("plot.confint_spm")

# Check that plot.confint_spm works

# No legend position via legend_pos

cis <- confint(res, interval_type = "both")
ciplot <- plot(cis)
test_that("plot.confint_spm works, sliding", {
  testthat::expect_identical(ciplot, NULL)
})

ciplot <- plot(cis, estimator = "BB2018", title = "BB2018 only")
test_that("plot.confint_spm works, sliding, BB2018 only, add leg title", {
  testthat::expect_identical(ciplot, NULL)
})

ciplot <- plot(cis, estimator = "BB2018b", legend = c("cool", "neat"))
test_that("plot.confint_spm works, sliding, BB2018b only, add legend", {
  testthat::expect_identical(ciplot, NULL)
})

ciplot <- plot(cis, estimator = c("N2015", "BB2018"),
               main = "N2015 and BB2018", legend = c("cool", "neat"),
               title = "2 ests")
test_that("plot.confint_spm works, sliding, 2 ests, user legend & title", {
  testthat::expect_identical(ciplot, NULL)
})

# Change legend position via legend_pos

cis <- confint(res, interval_type = "both")
ciplot <- plot(cis, legend_pos = "bottomright")
test_that("plot.confint_spm works, sliding", {
  testthat::expect_identical(ciplot, NULL)
})

ciplot <- plot(cis, estimator = "BB2018", title = "BB2018 only",
               legend_pos = "bottomright")
test_that("plot.confint_spm works, sliding, BB2018 only, add leg title & pos", {
  testthat::expect_identical(ciplot, NULL)
})

ciplot <- plot(cis, estimator = "BB2018b", legend = c("cool", "neat"),
               legend_pos = "bottomright")
test_that("plot.confint_spm works, sliding, BB2018b only, add leg & pos", {
  testthat::expect_identical(ciplot, NULL)
})

ciplot <- plot(cis, estimator = c("BB2018", "BB2018b"),
               main = "BB2018 and BB2018b", legend = c("cool", "neat"),
               title = "2 ests", legend_pos = "bottomright")
test_that("plot.confint_spm works, sliding, 2 ests, user leg & title & pos", {
  testthat::expect_identical(ciplot, NULL)
})

cis <- confint(res, interval_type = "both", maxima = "disjoint")
ciplot <- plot(cis, xlab = "my xlab", lwd = 2, col = "blue")
test_that("plot.confint_spm works, user plot args, disjoint", {
  testthat::expect_identical(ciplot, NULL)
})

# Extreme example where b is so small that the SEs for the sliding blocks
# version of the estimator cannot be calculated

# b = 7 makes BB2018 SE missing, and the plot should work
res_small_b <- spm(newlyn, 7)
cis <- confint(res_small_b, interval_type = "lik")
ciplot <- plot(cis, xlab = "my xlab", lwd = 2, col = "blue")
test_that("plot.confint_spm works, user plot args, disjoint", {
  testthat::expect_identical(ciplot, NULL)
})

# b = 4 makes both N2015 and BB2018 SE missing, and plot.confint_spm()
# should stop
res_small_b <- spm(newlyn, 4)
cis <- confint(res_small_b, interval_type = "lik")
ciplot <- try(plot(cis, xlab = "my xlab", lwd = 2, col = "blue"),
              silent = TRUE)
test_that("plot.confint_spm works, user plot args, disjoint", {
  testthat::expect_identical(class(ciplot), "try-error")
})

# ================================== kgaps ===================================

context("confint.kgaps")

u <- quantile(newlyn, probs = 0.90)

res <- kgaps(newlyn, u)
res1 <- confint(res)
res2 <- confint(res, conf_scale = "log")
test_that("kgaps lik intervals don't depend on conf_scale", {
  testthat::expect_identical(res1["lik", ], res2["lik", ])
})

# Repeat for inc_cens = TRUE

res <- kgaps(newlyn, u, inc_cens = TRUE)
res1 <- confint(res)
res2 <- confint(res, conf_scale = "log")
test_that("kgaps lik intervals don't depend on conf_scale", {
  testthat::expect_identical(res1["lik", ], res2["lik", ])
})

