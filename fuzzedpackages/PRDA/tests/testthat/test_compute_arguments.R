#####################################
####    Test compute arguments   ####
#####################################

library(PRDA)

#----    input checks    ----

context("Compute arguments")

x <- c(-32, -49, 157, 116, 114, 170, 162, 220, -55, 97, 151, 160, 145, 78, 116)
y <- c(-30, 102, 184, 28, -89, -40, -43, -129, 53, 119, 44, 35, 189, 137, -68)
diff <- x-y
mx <- mean(x)
my <- mean(y)
nx <- length(x)
ny <- length(y)
ny2 <- ny+2
sig_level <- .10
mu <- .5
mu2 <- .5

# For case ratio_sd
var1 <- 2^2
var2 <- 1


#----    compute_df    ----

test_that("evaluate the correct df", {
  # Redefine function to avoid specify arguments each the times
  test_compute_df <- function(effect_type = "cohen_d", sample_n1 = nx,
                              sample_n2 = NULL, ratio_sd = 1, test_method){
    compute_df(effect_type, sample_n1, sample_n2, ratio_sd, test_method)
  }

  expect_equal(test_compute_df(sample_n2 = ny2, test_method = "one_sample"),
               nx-1L)
  expect_equal(test_compute_df(sample_n2 = ny, test_method = "paired"),
               nx-1L)
  expect_equal(test_compute_df(sample_n2 = ny2, test_method = "two_sample"),
               nx+ny2-2L)
  # case ratio_sd = 1
  expect_equal(test_compute_df(sample_n2 = ny2, test_method = "welch"),
               ((nx+ny2)^2*(nx-1)*(ny2-1))/(nx^2*(nx-1)+ ny2^2*(ny2-1)))
  # case ratio_sd = 2
  expect_equal(test_compute_df(sample_n2 = ny2, test_method = "welch", ratio_sd = 2),
               (var1/nx + var2/ny2)^2/((var1/nx)^2/(nx-1) + (var2/ny2)^2/(ny2-1)))

  expect_equal(test_compute_df(effect_type = "correlation", sample_n1 = nx),
               nx-2L)
  })

#----    compute_critical_effect    ----

test_that("evaluate the correct critical effect", {
  # Redefine function to avoid specify arguments each the times
  test_compute_critical_effect <- function(effect_type = "cohen_d",
              sample_n1 = nx, sample_n2 = NULL, test_method, sig_level = .10,
              alternative, ratio_sd = 1, mu = 0){
    compute_critical_effect(effect_type, sample_n1, sample_n2, test_method,
                            sig_level, alternative, ratio_sd, mu)
  }

  expect_equal(test_compute_critical_effect(
    sample_n2 = ny2, test_method = "one_sample", alternative = "two_sided"),
    list(df = nx-1, critical_effect = qt(1-sig_level/2, df = nx-1)/sqrt(nx)))

  expect_equal(test_compute_critical_effect(
    sample_n2 = ny, test_method = "paired", alternative = "less"),
    list(df = nx-1, critical_effect = qt(sig_level, df = nx-1)/sqrt(nx)))
  expect_equal(test_compute_critical_effect(
    sample_n2 = ny, test_method = "paired", alternative = "greater", mu = mu),
    list(df = nx-1, critical_effect = qt(1-sig_level, df = nx-1)/sqrt(nx)+mu))

  expect_equal(test_compute_critical_effect(
    sample_n2 = ny2, test_method = "two_sample", alternative = "two_sided"),
    list(df = nx+ny2-2, critical_effect = qt(1-sig_level/2, df = nx+ny2-2)*sqrt((nx+ny2)/(nx*ny2))))
  expect_equal(test_compute_critical_effect(
    sample_n2 = ny2, test_method = "two_sample", alternative = "less", mu = mu2),
    list(df = nx+ny2-2, critical_effect = qt(sig_level, df = nx+ny2-2)*sqrt((nx+ny2)/(nx*ny2)) + mu2))

  # case ratio_sd = 1
  expect_equal(test_compute_critical_effect(
    sample_n2 = ny2, test_method = "welch", alternative = "greater"),
    list(df = ((nx+ny2)^2 * (nx-1) * (ny2-1))/(nx^2*(nx-1) + ny2^2*(ny2-1)),
         critical_effect = qt(1-sig_level, df = ((nx+ny2)^2 * (nx-1) * (ny2-1))/
                                (nx^2*(nx-1) + ny2^2*(ny2-1)))*sqrt((nx+ny2)/(nx*ny2))))
  expect_equal(test_compute_critical_effect(
    sample_n2 = ny2, test_method = "welch", alternative = "two_sided", mu = mu2),
    list(df = ((nx+ny2)^2 * (nx-1) * (ny2-1))/(nx^2*(nx-1) + ny2^2*(ny2-1)),
         critical_effect = qt(1-sig_level/2, df = ((nx+ny2)^2 * (nx-1) * (ny2-1))/
                                (nx^2*(nx-1) + ny2^2*(ny2-1)))*sqrt((nx+ny2)/(nx*ny2)) + mu2))

  # case ratio_sd = 2
  expect_equal(test_compute_critical_effect(
    sample_n2 = ny2, test_method = "welch", alternative = "greater", ratio_sd = 2),
    list(df = (var1/nx + var2/ny2)^2/((var1/nx)^2/(nx-1) + (var2/ny2)^2/(ny2-1)),
         critical_effect = qt(1-sig_level, df = (var1/nx + var2/ny2)^2/
                                ((var1/nx)^2/(nx-1) + (var2/ny2)^2/(ny2-1)))*sqrt(2/(nx*ny2) * (var1*ny2 +var2*nx)/(var1 + var2))))

  expect_equal(test_compute_critical_effect(
    effect_type = "correlation", alternative = "two_sided"),
    list(df = nx-2, critical_effect = qt(1-sig_level/2, df = nx-2)/sqrt(nx-2+qt(1-sig_level/2, df = nx-2)^2)))
  expect_equal(test_compute_critical_effect(
    effect_type = "correlation", alternative = "greater"),
    list(df = nx-2, critical_effect = qt(1-sig_level, df = nx-2)/sqrt(nx-2+qt(1-sig_level, df = nx-2)^2)))
  expect_equal(test_compute_critical_effect(
    effect_type = "correlation", alternative = "less"),
    list(df = nx-2, critical_effect = qt(sig_level, df = nx-2)/sqrt(nx-2+qt(sig_level, df = nx-2)^2)))

  })

#----
