context("PAsso: diagnostic.plot")

test_that("diagnostic.plot works for \"PAsso\" objects", {

  # Skips
  skip_on_cran()

  # Load data
  data("ANES2016")

  PAsso_3v <- PAsso(responses = c("PreVote.num", "PID", "selfLR"),
                   adjustments = c("income.num", "age", "edu.year"),
                   data = ANES2016, uni.model = "probit",
                   method = c("kendall"),
                   resids.type = "surrogate", jitter = "latent")

  diag_p1 <- diagnostic.plot(object = PAsso_3v, output = "qq")
  diag_p2 <- diagnostic.plot(object = PAsso_3v, output = "fitted")
  suppressWarnings(
    diag_p3 <- diagnostic.plot(object = PAsso_3v, output = "covariate")
  )
  # Expectations
  expect_s3_class(diag_p1, "gtable")
  expect_s3_class(diag_p1, "gtable")
  expect_s3_class(diag_p1, "gtable")
  expect_s3_class(diag_p2, "gtable")
  expect_s3_class(diag_p3, "gtable")

  # expect_s3_class(diag_p2, "ggmatrix")
  # expect_s3_class(diag_p3, "ggmatrix")

})

test_that("diagnostic.plot works for \"PAsso\" objects with 5 responses", {

  # Skips
  skip_on_cran()

  # multivariate analysis (5 variables) --------------------------------------------------------------------
  PAsso_5v <- PAsso(responses = c("PreVote.num", "PID", "selfLR", "ClinLR", "TrumpLR"),
                    adjustments = c("income.num", "age", "edu.year"),
                    data = ANES2016, uni.model = "logit",
                    method = c("kendall"),
                    resids.type = "surrogate", jitter = "latent")

  diag_p1 <- diagnostic.plot(object = PAsso_5v, output = "qq")

  # Expectations
  expect_s3_class(diag_p1, "gtable")

})

test_that("diagnostic.plot work for \"clm\" objects", {

  # Skips
  skip_on_cran()
  skip_if_not_installed("ordinal")

  # Load data
  data(df1)

  # Fit cumulative link model
  fit <- ordinal::clm(y ~ x + I(x ^ 2), data = df1, link = "logit")

  # diagnostic.plot
  plot_qq <- diagnostic.plot(object = fit, output = "qq")
  plot_fit <- diagnostic.plot(object = fit, output = "fitted")
  plot_cov <- diagnostic.plot(object = fit, output = "covariate")

  # Expectations
  expect_is(plot_qq, "ggplot")
  expect_is(plot_fit, "ggplot")
  expect_is(plot_cov, "ggplot")

})

test_that("diagnostic.plot work for \"glm\" objects", {

  # Skips
  skip_on_cran()

  # Load data
  data(df1)

  # Fit cumulative link model
  fit <- glm(y ~ x + I(x ^ 2), data = df1, family = binomial)

  # diagnostic.plot
  plot_qq <- diagnostic.plot(object = fit, output = "qq")
  plot_fit <- diagnostic.plot(object = fit, output = "fitted")
  plot_cov <- diagnostic.plot(object = fit, output = "covariate")

  # Expectations
  expect_is(plot_qq, "ggplot")
  expect_is(plot_fit, "ggplot")
  expect_is(plot_cov, "ggplot")


  data("ANES2016")

  fit1 <- ordinal::clm(formula = as.factor(PreVote.num) ~ income.num + age + edu.year,
                       data = ANES2016, link = "logit")

  # diagnostic.plot
  plot_qq_1 <- diagnostic.plot(object = fit1, output = "qq")
  plot_fit_1 <- diagnostic.plot(object = fit1, output = "fitted")
  plot_cov_1 <- diagnostic.plot(object = fit1, output = "covariate")
  # Expectations
  expect_is(plot_qq_1, "ggplot")
  expect_is(plot_qq_1, "ggplot")
  expect_is(plot_qq_1, "ggplot")

})

test_that("diagnostic.plot work for \"lrm\" objects", {

  # Skips
  skip_on_cran()
  skip_if_not_installed("rms")

  # Load data
  data(df1)

  # Fit cumulative link model
  fit <- rms::lrm(y ~ x, data = df1)

  # diagnostic.plot
  plot_qq <- diagnostic.plot(object = fit, output = "qq")
  plot_fit <- diagnostic.plot(object = fit, output = "fitted")
  plot_cov <- diagnostic.plot(object = fit, output = "covariate")

  # Expectations
  expect_is(plot_qq, "ggplot")
  expect_is(plot_fit, "ggplot")
  expect_is(plot_cov, "ggplot")
})

test_that("diagnostic.plot work for \"orm\" objects", {

  # Skips
  skip_on_cran()
  skip_if_not_installed("rms")

  # Load data
  data(df1)

  # Fit cumulative link model
  fit <- rms::orm(y ~ x, data = df1, family = logistic)

  # diagnostic.plot
  plot_qq <- diagnostic.plot(object = fit, output = "qq")
  plot_fit <- diagnostic.plot(object = fit, output = "fitted")
  plot_cov <- diagnostic.plot(object = fit, output = "covariate")

  # Expectations
  expect_is(plot_qq, "ggplot")
  expect_is(plot_fit, "ggplot")
  expect_is(plot_cov, "ggplot")
})

test_that("diagnostic.plot work for \"polr\" objects", {

  # Skips
  skip_on_cran()
  skip_if_not_installed("MASS")

  # Load data
  data(df1)

  # Fit cumulative link model
  fit <- MASS::polr(y ~ x + I(x ^ 2), data = df1, method = "logistic")

  # diagnostic.plot
  plot_qq <- diagnostic.plot(object = fit, output = "qq")
  plot_fit <- diagnostic.plot(object = fit, output = "fitted")
  plot_cov <- diagnostic.plot(object = fit, output = "covariate")

  # Expectations
  expect_is(plot_qq, "ggplot")
  expect_is(plot_fit, "ggplot")
  expect_is(plot_cov, "ggplot")
})

test_that("diagnostic.plot work for \"vglm\" objects", {

  # Skips
  skip_on_cran()
  skip_if_not_installed("VGAM")

  # Load data
  data(df1)

  # Fit cumulative link model
  suppressWarnings(
    fit <- VGAM::vglm(y ~ x + I(x ^ 2), data = df1,
                      family = VGAM::cumulative(link = "logit",
                                                parallel = TRUE))
  )

  # Expectations
  expect_warning(diagnostic.plot(object = fit, output = "qq"), "does not know")
  expect_warning(diagnostic.plot(object = fit, output = "fitted"), "does not know")
  expect_warning(diagnostic.plot(object = fit, output = "covariate"), "does not know")

})

test_that("diagnostic.plot work for \"clm\" objects with different link functions", {

  # Skips
  skip_on_cran()
  skip_if_not_installed("ordinal")

  # Load data
  data(df1)

  # Fit cumulative link models
  fit1 <- ordinal::clm(y ~ x + I(x ^ 2), data = df1, link = "logit")
  fit2 <- ordinal::clm(y ~ x + I(x ^ 2), data = df1, link = "probit")
  fit3 <- ordinal::clm(y ~ x + I(x ^ 2), data = df1, link = "loglog")
  fit4 <- ordinal::clm(y ~ x + I(x ^ 2), data = df1, link = "cloglog")
  fit5 <- ordinal::clm(y ~ x + I(x ^ 2), data = df1, link = "cauchit")

  # diagnostic.plot
  plot_qq_1 <- diagnostic.plot(object = fit1, output = "qq")
  plot_fit_1 <- diagnostic.plot(object = fit1, output = "fitted")
  plot_cov_1 <- diagnostic.plot(object = fit1, output = "covariate")
  plot_qq_2 <- diagnostic.plot(object = fit2, output = "qq")
  plot_fit_2 <- diagnostic.plot(object = fit2, output = "fitted")
  plot_cov_2 <- diagnostic.plot(object = fit2, output = "covariate")
  plot_qq_3 <- diagnostic.plot(object = fit3, output = "qq")
  plot_fit_3 <- diagnostic.plot(object = fit3, output = "fitted")
  plot_cov_3 <- diagnostic.plot(object = fit3, output = "covariate")
  plot_qq_4 <- diagnostic.plot(object = fit4, output = "qq")
  plot_fit_4 <- diagnostic.plot(object = fit4, output = "fitted")
  plot_cov_4 <- diagnostic.plot(object = fit4, output = "covariate")
  plot_qq_5 <- diagnostic.plot(object = fit5, output = "qq")
  plot_fit_5 <- diagnostic.plot(object = fit5, output = "fitted")
  plot_cov_5 <- diagnostic.plot(object = fit5, output = "covariate")

  # Expectations
  expect_is(plot_qq_1, "ggplot")
  expect_is(plot_fit_1, "ggplot")
  expect_is(plot_cov_1, "ggplot")
  expect_is(plot_qq_2, "ggplot")
  expect_is(plot_fit_2, "ggplot")
  expect_is(plot_cov_2, "ggplot")
  expect_is(plot_qq_3, "ggplot")
  expect_is(plot_fit_3, "ggplot")
  expect_is(plot_cov_3, "ggplot")
  expect_is(plot_qq_4, "ggplot")
  expect_is(plot_fit_4, "ggplot")
  expect_is(plot_cov_4, "ggplot")
  expect_is(plot_qq_5, "ggplot")
  expect_is(plot_fit_5, "ggplot")
  expect_is(plot_cov_5, "ggplot")
})
