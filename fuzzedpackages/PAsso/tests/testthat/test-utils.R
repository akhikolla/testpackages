context("sure: Utility functions")


test_that("utility functions work for \"clm\" objects", {

  # Skips
  skip_on_cran()
  skip_if_not_installed("ordinal")

  # Load data
  data(df1)

  # Fit cumulative link model
  fit_logit <- ordinal::clm(y ~ x + I(x ^ 2), data = df1, link = "logit")
  fit_probit <- ordinal::clm(y ~ x + I(x ^ 2), data = df1, link = "probit")
  fit_loglog <- ordinal::clm(y ~ x + I(x ^ 2), data = df1, link = "loglog")
  fit_cloglog <- ordinal::clm(y ~ x + I(x ^ 2), data = df1, link = "cloglog")
  fit_cauchit <- ordinal::clm(y ~ x + I(x ^ 2), data = df1, link = "cauchit")

  # Expectations
  expect_equal(length(getBounds(fit_logit)), 5)
  expect_identical(getResponseValues(fit_logit), as.integer(df1$y))
  expect_equal(ncat(fit_logit), 4)
  expect_equal(getDistributionFunction(fit_logit), plogis)
  expect_equal(getDistributionFunction(fit_probit), pnorm)
  expect_equal(getDistributionFunction(fit_loglog), pgumbel)
  expect_equal(getDistributionFunction(fit_cloglog), pGumbel)
  expect_equal(getDistributionFunction(fit_cauchit), pcauchy)
  expect_equal(getDistributionName(fit_logit), "logis")
  expect_equal(getDistributionName(fit_probit), "norm")
  expect_equal(getDistributionName(fit_loglog), "gumbel")
  expect_equal(getDistributionName(fit_cloglog), "Gumbel")
  expect_equal(getDistributionName(fit_cauchit), "cauchy")
  expect_equal(getQuantileFunction(fit_logit), qlogis)
  expect_equal(getQuantileFunction(fit_probit), qnorm)
  expect_equal(getQuantileFunction(fit_loglog), qgumbel)
  expect_equal(getQuantileFunction(fit_cloglog), qGumbel)
  expect_equal(getQuantileFunction(fit_cauchit), qcauchy)

})


test_that("utility functions work for \"glm\" objects", {

  # Skips
  skip_on_cran()

  # Load data
  data(df1)

  # Fit binary probit model
  # BUG? only works for binary? But the response has 4 levels
  suppressWarnings(
    fit <- stats::glm(y ~ x + I(x ^ 2), data = df1,
                      family = binomial(link = "probit"))
  )

  # Expectations
  # expect_null(getBounds(fit))
  expect_equal(getDistributionFunction(fit), pnorm)
  expect_equal(getDistributionName(fit), "norm")
  expect_equal(getQuantileFunction(fit), qnorm)
  expect_identical(getResponseValues(fit), as.integer(df1$y))
  # expect_equal(ncat(fit), 4)

})


test_that("utility functions work for \"lrm\" objects", {

  # Skips
  skip_on_cran()
  skip_if_not_installed("rms")

  # Load data
  data(df1)

  # Fit cumulative link model
  fit <- rms::lrm(y ~ x, data = df1)

  # Expectations
  expect_equal(length(getBounds(fit)), 5)
  expect_equal(getDistributionFunction(fit), plogis)
  expect_equal(getDistributionName(fit), "logis")
  expect_equal(getQuantileFunction(fit), qlogis)
  expect_identical(getResponseValues(fit), as.integer(df1$y))
  expect_equal(ncat(fit), 4)

})


test_that("utility functions work for \"orm\" objects", {

  # Skips
  skip_on_cran()
  skip_if_not_installed("rms")

  # Load data
  data(df1)

  # Fit cumulative link model
  fit <- rms::orm(y ~ x, data = df1, family = probit)

  # Fit cumulative link models
  fit_logit <- rms::orm(y ~ x, data = df1, family = logistic)
  fit_probit <- rms::orm(y ~ x, data = df1, family = probit)
  fit_loglog <- rms::orm(y ~ x, data = df1, family = loglog)
  fit_cloglog <- rms::orm(y ~ x, data = df1, family = cloglog)
  fit_cauchit <- rms::orm(y ~ x, data = df1, family = cauchit)

  # Expectations
  expect_equal(length(getBounds(fit_logit)), 5)
  expect_identical(getResponseValues(fit_logit), as.integer(df1$y))
  expect_equal(ncat(fit_logit), 4)
  expect_equal(getDistributionFunction(fit_logit), plogis)
  expect_equal(getDistributionFunction(fit_probit), pnorm)
  expect_equal(getDistributionFunction(fit_loglog), pgumbel)
  expect_equal(getDistributionFunction(fit_cloglog), pGumbel)
  expect_equal(getDistributionFunction(fit_cauchit), pcauchy)
  expect_equal(getDistributionName(fit_logit), "logis")
  expect_equal(getDistributionName(fit_probit), "norm")
  expect_equal(getDistributionName(fit_loglog), "gumbel")
  expect_equal(getDistributionName(fit_cloglog), "Gumbel")
  expect_equal(getDistributionName(fit_cauchit), "cauchy")
  expect_equal(getQuantileFunction(fit_logit), qlogis)
  expect_equal(getQuantileFunction(fit_probit), qnorm)
  expect_equal(getQuantileFunction(fit_loglog), qgumbel)
  expect_equal(getQuantileFunction(fit_cloglog), qGumbel)
  expect_equal(getQuantileFunction(fit_cauchit), qcauchy)

})


test_that("utility functions work for \"polr\" objects", {

  # Skips
  skip_on_cran()
  skip_if_not_installed("MASS")

  # Load data
  data(df1)

  # Fit cumulative link models
  fit_logit <- MASS::polr(y ~ x + I(x ^ 2), data = df1, method = "logistic")
  fit_probit <- MASS::polr(y ~ x + I(x ^ 2), data = df1, method = "probit")
  fit_loglog <- MASS::polr(y ~ x + I(x ^ 2), data = df1, method = "loglog")
  fit_cloglog <- MASS::polr(y ~ x + I(x ^ 2), data = df1, method = "cloglog")
  fit_cauchit <- MASS::polr(y ~ x + I(x ^ 2), data = df1, method = "cauchit")

  # Expectations
  expect_equal(length(getBounds(fit_logit)), 5)
  expect_identical(getResponseValues(fit_logit), as.integer(df1$y))
  expect_equal(ncat(fit_logit), 4)
  expect_equal(getDistributionFunction(fit_logit), plogis)
  expect_equal(getDistributionFunction(fit_probit), pnorm)
  expect_equal(getDistributionFunction(fit_loglog), pgumbel)
  expect_equal(getDistributionFunction(fit_cloglog), pGumbel)
  expect_equal(getDistributionFunction(fit_cauchit), pcauchy)
  expect_equal(getDistributionName(fit_logit), "logis")
  expect_equal(getDistributionName(fit_probit), "norm")
  expect_equal(getDistributionName(fit_loglog), "gumbel")
  expect_equal(getDistributionName(fit_cloglog), "Gumbel")
  expect_equal(getDistributionName(fit_cauchit), "cauchy")
  expect_equal(getQuantileFunction(fit_logit), qlogis)
  expect_equal(getQuantileFunction(fit_probit), qnorm)
  expect_equal(getQuantileFunction(fit_loglog), qgumbel)
  expect_equal(getQuantileFunction(fit_cloglog), qGumbel)
  expect_equal(getQuantileFunction(fit_cauchit), qcauchy)

})


test_that("utility functions work for \"vglm\" objects", {

  # Skips
  skip_on_cran()
  skip_if_not_installed("VGAM")

  # Load data
  data(df1)

  # Fit cumulative link models
  suppressWarnings(
    fit_logit <- VGAM::vglm(y ~ x + I(x ^ 2), data = df1,
                            family = VGAM::cumulative(link = "logit",
                                                      parallel = TRUE))
  )
  suppressWarnings(
    fit_probit <- VGAM::vglm(y ~ x + I(x ^ 2), data = df1,
                             family = VGAM::cumulative(link = "probit",
                                                       parallel = TRUE))
  )
  # fit_loglog <- VGAM::vglm(y ~ x + I(x ^ 2), data = df1,
  #                          family = VGAM::cumulative(link = "loglog",
  #                                                    parallel = TRUE))
  suppressWarnings(
    fit_cloglog <- VGAM::vglm(y ~ x + I(x ^ 2), data = df1,
                              family = VGAM::cumulative(link = "cloglog",
                                                        parallel = TRUE))
  )
  suppressWarnings(
    fit_cauchit <- VGAM::vglm(y ~ x + I(x ^ 2), data = df1,
                              family = VGAM::cumulative(link = "cauchit",
                                                        parallel = TRUE))
  )

  # Fit an adjacent categories regression model
  suppressWarnings(
    fit_acat <- VGAM::vglm(y ~ x + I(x ^ 2), data = df1,
                            family = VGAM::acat(reverse=TRUE, parallel=TRUE))
  )

  # Expectations
  expect_equal(length(getBounds(fit_logit)), 5)
  expect_identical(getResponseValues(fit_logit), as.integer(df1$y))
  expect_equal(ncat(fit_logit), 4)
  expect_equal(getDistributionFunction(fit_logit), plogis)
  expect_equal(getDistributionFunction(fit_probit), pnorm)
  # expect_equal(getDistributionFunction(fit_loglog), pgumbel)
  expect_equal(getDistributionFunction(fit_cloglog), pGumbel)
  expect_equal(getDistributionFunction(fit_cauchit), pcauchy)
  expect_equal(getDistributionName(fit_logit), "logis")
  expect_equal(getDistributionName(fit_probit), "norm")
  # expect_equal(getDistributionName(fit_loglog), "gumbel")
  expect_equal(getDistributionName(fit_cloglog), "Gumbel")
  expect_equal(getDistributionName(fit_cauchit), "cauchy")
  expect_equal(getQuantileFunction(fit_logit), qlogis)
  expect_equal(getQuantileFunction(fit_probit), qnorm)
  # expect_equal(getQuantileFunction(fit_loglog), qgumbel)
  expect_equal(getQuantileFunction(fit_cloglog), qGumbel)
  expect_equal(getQuantileFunction(fit_cauchit), qcauchy)

})


test_that("generate_residuals function works for different methods", {

  # Skips
  skip_on_cran()

  # Load data
  data(df1)

  # Fit cumulative link models
  suppressWarnings(
    fit_clm <- ordinal::clm(y ~ x, data = df1, link = "logit")
  )
  suppressWarnings(
    resids_clm <- generate_residuals(fit_clm, method = "sign", draws_id = NULL)
  )
  suppressWarnings(
    resids_clm_boot <-
      generate_residuals(fit_clm, method = "sign",
                         draws_id = sample(nobs(fit_clm), replace = TRUE))
  )


  fit_polr <- MASS::polr(y ~ x, data = df1, method = "logistic")
  resids_polr <- generate_residuals(fit_polr, method = "sign", draws_id = NULL)
  resids_polr_boot <-
    generate_residuals(fit_polr, method = "sign",
                       draws_id = sample(nobs(fit_polr), replace = TRUE))

  fit_lrm <- rms::lrm(y ~ x, data = df1)
  resids_lrm <- generate_residuals(fit_lrm, method = "sign", draws_id = NULL)
  resids_lrm_boot <-
    generate_residuals(fit_lrm, method = "sign",
                       draws_id = sample(nobs(fit_lrm), replace = TRUE))

  fit_orm <- rms::orm(y ~ x, data = df1, family = logistic)
  resids_orm <- generate_residuals(fit_orm, method = "sign", draws_id = NULL)
  resids_orm_boot <-
    generate_residuals(fit_orm, method = "sign",
                       draws_id = sample(nobs(fit_orm), replace = TRUE))

  fit_vglm <- VGAM::vglm(y ~ x, data = df1,
                         family = VGAM::cumulative(link = "logit",
                                                   parallel = TRUE))
  resids_vglm <- generate_residuals(fit_vglm, method = "sign", draws_id = NULL)
  resids_vglm_boot <-
    generate_residuals(fit_vglm, method = "sign",
                       draws_id = sample(nobs(fit_vglm), replace = TRUE))


  # Expectations
  expect_equal(length(resids_clm), nrow(df1))
  expect_equal(length(resids_clm_boot), nrow(df1))
  expect_equal(length(resids_polr), nrow(df1))
  expect_equal(length(resids_polr_boot), nrow(df1))
  expect_equal(length(resids_lrm), nrow(df1))
  expect_equal(length(resids_lrm_boot), nrow(df1))
  expect_equal(length(resids_orm), nrow(df1))
  expect_equal(length(resids_orm_boot), nrow(df1))
  expect_equal(length(resids_vglm), nrow(df1))
  expect_equal(length(resids_vglm_boot), nrow(df1))

  expect_null(attr(resids_clm, "draws"))
  expect_null(attr(resids_clm, "draws_id"))
  expect_null(attr(resids_vglm_boot, "draws"))
  expect_null(attr(resids_vglm_boot, "draws_id"))

})


test_that("getMeanResponse works", {

  # Skips
  skip_on_cran()
  skip_if_not_installed("ordinal")
  skip_if_not_installed("rms")
  skip_if_not_installed("MASS")
  skip_if_not_installed("VGAM")

  # Load data
  data(df1)

  # Fit cumulative link models
  fit_clm <- ordinal::clm(y ~ x, data = df1, link = "logit")
  fit_polr <- MASS::polr(y ~ x, data = df1, method = "logistic")
  fit_lrm <- rms::lrm(y ~ x, data = df1)
  fit_orm <- rms::orm(y ~ x, data = df1, family = logistic)
  fit_vglm <- VGAM::vglm(y ~ x, data = df1,
                         family = VGAM::cumulative(link = "logit",
                                                   parallel = TRUE))
  # fit_vglm@misc$reverse
  # unname(stats::coef(fit_clm))
  # unname(stats::coef(fit_polr))
  # unname(stats::coef(fit_lrm))
  # unname(stats::coef(fit_orm))
  # unname(stats::coef(fit_vglm))


  # Mean response
  mr <- cbind(
    "clm" = getMeanResponse(fit_clm),
    "polr" = getMeanResponse(fit_polr),
    "lrm" = getMeanResponse(fit_lrm),
    "orm" = getMeanResponse(fit_orm),
    "vglm" = getMeanResponse(fit_vglm)
  )

  # Compute maximum pairwise difference per row
  max_diff <- apply(mr, MARGIN = 1, FUN = function(x) max(as.numeric(dist(x))))

  # Expectations
  expect_true(max(max_diff) < 1e-05)

})
