context("spm_R_quick vs spm_check vs spm")

# Use only a subset of the newlyn data, for speed.
# With b = 180 there are only two sets of disjoint block maxima

test_data <- newlyn[1:2881]

# Check that spm_R_quick(), faster but not very transparent, gives the same
# results as spm_check(), slower but more transparent
# Check that spm(), which uses Rcpp, gives the same results as spm_R_quick()

# 8 cases:
# bias_adjust in c("BB3", "BB1", "N", "none")
bias_adjust_vec <- c("BB3", "BB1", "N", "none")
# which_dj in c("last", "first")
which_dj_vec <- c("last", "first")
# block size: pick a big one so that the tests aren't slow
# The permitted range of b for these data is 15 - 196
# We must respect this here because spm_check() doesn't check b
# It also doesn't check the format of the data
b <- 180
# Tolerance
my_tol <- 1e-5

for (i in 1:4){
  for (j in 1:2) {
    res <- spm_R_quick(test_data, b = b,
                       bias_adjust = bias_adjust_vec[i],
                       which_dj = which_dj_vec[j])
    res_sl <- spm_check(test_data, b = b, sliding = TRUE,
                        bias_adjust = bias_adjust_vec[i],
                        which_dj = which_dj_vec[j])
    res_dj <- spm_check(test_data, b = b, sliding = FALSE,
                        bias_adjust = bias_adjust_vec[i],
                        which_dj = which_dj_vec[j])
    # spm_R_quick vs spm_check
    my_text <- paste("spm_R_quick vs spm_check", bias_adjust_vec[i],
                     which_dj_vec[j])
    test_that(paste(my_text, "sliding, theta"), {
      testthat::expect_equal(res$theta_sl, res_sl$theta, tolerance = my_tol)
    })
    test_that(paste(my_text, "disjoint, theta"), {
      testthat::expect_equal(res$theta_dj, res_dj$theta, tolerance = my_tol)
    })
    test_that(paste(my_text, "sliding, se"), {
      testthat::expect_equal(res$se_sl, res_sl$se, tolerance = my_tol)
    })
    test_that(paste(my_text, "disjoint, se"), {
      testthat::expect_equal(res$se_dj, res_dj$se, tolerance = my_tol)
    })
    test_that(paste(my_text, "sliding, bias"), {
      testthat::expect_equal(res$bias_sl, res_sl$bias_val, tolerance = my_tol)
    })
    test_that(paste(my_text, "disjoint, bias"), {
      testthat::expect_equal(res$bias_dj, res_dj$bias_val, tolerance = my_tol)
    })
    test_that(paste(my_text, "sliding, data_sl"), {
      testthat::expect_equal(summary(res$data_sl),
                             summary(cbind(N2015 = res_sl$N2015_data,
                                           BB2018 = res_sl$BB2018_data)),
                             tolerance = my_tol)
    })
    test_that(paste(my_text, "sliding, data_dj"), {
      testthat::expect_equal(summary(res$data_dj),
                             summary(cbind(N2015 = res_dj$N2015_data,
                                           BB2018 = res_dj$BB2018_data)),
                             tolerance = my_tol)
    })
    # spm_R_quick vs spm
    # Note: only spm() returns estimates of BB2018b (BB2018 - 1 / b),
    #       hence the use of [1:2] below
    my_text <- paste("spm vs spm_R_quick", bias_adjust_vec[i],
                     which_dj_vec[j])
    res_c <- spm(test_data, b = b,
                 bias_adjust = bias_adjust_vec[i],
                 which_dj = which_dj_vec[j])
    test_that(paste(my_text, "sliding, theta"), {
      testthat::expect_equal(res$theta_sl, res_c$theta_sl[1:2],
                             tolerance = my_tol)
    })
    test_that(paste(my_text, "disjoint, theta"), {
      testthat::expect_equal(res$theta_dj, res_c$theta_dj[1:2],
                             tolerance = my_tol)
    })
    test_that(paste(my_text, "sliding, se"), {
      testthat::expect_equal(res$se_sl, res_c$se_sl[1:2],
                             tolerance = my_tol)
    })
    test_that(paste(my_text, "disjoint, se"), {
      testthat::expect_equal(res$se_dj, res_c$se_dj[1:2],
                             tolerance = my_tol)
    })
    test_that(paste(my_text, "sliding, bias"), {
      testthat::expect_equal(res$bias_sl, res_c$bias_sl[1:2],
                             tolerance = my_tol)
    })
    test_that(paste(my_text, "disjoint, bias"), {
      testthat::expect_equal(res$bias_dj, res_c$bias_dj[1:2],
                             tolerance = my_tol)
    })
    test_that(paste(my_text, "sliding, data_sl"), {
      testthat::expect_equal(summary(res$data_sl),
                             summary(res_c$data_sl),
                             tolerance = my_tol)
    })
    test_that(paste(my_text, "sliding, data_dj"), {
      testthat::expect_equal(summary(res$data_dj),
                             summary(res_c$data_dj),
                             tolerance = my_tol)
    })
  }
}

###############################################################################

context("spm: equivalence of BB2018 when bias_adjust = ''BB1'' and ''N''")

b <- 180
resBB1 <- spm_R_quick(test_data, b = b, bias_adjust = "BB1")
resN <- spm_R_quick(test_data, b = b, bias_adjust = "N")

test_that(paste("BB1 vs N, b is OK"), {
  testthat::expect_equal(resBB1$theta_dj["BB2018"],
                         resN$theta_dj["BB2018"], tolerance = my_tol)
})

# ============================ summary.spm ====================================

context("summary.spm")

theta <- spm(newlyn, 20)
res <- summary(theta)
test_that(paste("No warning when b is large enough"), {
  testthat::expect_identical(res$warning, NULL)
})

theta <- spm(newlyn, 7)
res <- summary(theta)
test_that(paste("Warning when b is small enough"), {
  testthat::expect_identical(class(res$warning), "character")
})

