context("PAsso: S3 method 'residuals()' for Surrogate residuals")

test_that("residuals work for \"PAsso\" objects", {

  # Skips
  skip_on_cran()
  skip_if_not_installed("ordinal")

  # Load data
  data("ANES2016")

  # Load data
  PAsso_1 <- PAsso(responses = c("PreVote.num", "PID"),
                   adjustments = c("income.num", "age", "edu.year"),
                   data = ANES2016
                   # association = c("partial"),
                   # models = c("probit", "probit"),
                   # method = c("kendall"),
                   # resids.method = "latent", fitted.models = NULL,
                   # rep_num = 20
  )

  # Compute residuals
  res1 <- residuals(PAsso_1)

  # Expectations
  expect_equal(dim(res1)[1], nrow(ANES2016))
  expect_is(attr(res1, "draws"), "array")
  expect_null(attr(res1, "draws_id"))
  expect_equal(dim(attr(res1, "draws")), c(nrow(ANES2016), 20, 2))

  expect_equal(attr(res1, "arguments"), c("partial", "kendall", "surrogate", "latent", "probability"))

  # USE adjacent_cate of Shaobo;
  # USE draw_id, in residuals.PAsso

})

test_that("residuals work for \"clm\" objects", {

  # Skips
  skip_on_cran()
  skip_if_not_installed("ordinal")

  # Load data
  data(df1)

  # Fit cumulative link model
  fit <- ordinal::clm(y ~ x + I(x ^ 2), data = df1, link = "logit")

  # Compute residuals
  res1 <- residuals(fit)
  res2 <- residuals(fit, nsim = 10)
  res_PR1 <- residuals(fit, type = "sign")
  res_PR2 <- residuals(fit, type = "sign", nsim = 10)
  res_GR1 <- residuals(fit, type = "general")
  res_GR2 <- residuals(fit, type = "general", nsim = 10)
  res_DR1 <- residuals(fit, type = "deviance")
  res_DR2 <- residuals(fit, type = "deviance", nsim = 10)


  # Expectations
  expect_equal(length(res1), nrow(df1))
  expect_equal(length(res2), nrow(df1))
  expect_null(attr(res1, "draws"))
  expect_null(attr(res1, "draws_id"))
  expect_is(attr(res2, "draws"), "matrix")
  expect_is(attr(res2, "draws_id"), "matrix")
  expect_equal(dim(attr(res2, "draws")), c(nrow(df1), 10))
  expect_equal(dim(attr(res2, "draws_id")), c(nrow(df1), 10))

  expect_equal(attr(res1, "arguments")[1], "surrogate")
  expect_equal(attr(res_PR1, "arguments")[1], "sign")
  expect_equal(attr(res_GR1, "arguments")[1], "general")
  expect_equal(attr(res_DR1, "arguments")[1], "deviance")

})

test_that("stats::residuals.glm work for \"glm\" objects with ordinal response, no surrogate approach", {

  # Skips
  skip_on_cran()

  # Load data
  data(df1)

  # Fit cumulative link model
  fit <- glm(y ~ x + I(x ^ 2), data = df1, family = binomial)

  # Compute residuals
  res1 <- residuals(fit)
  res2 <- residuals(fit, nsim = 10)

  # Expectations
  expect_error(res3 <- residuals(fit, type = "surrogate", jitter = "uniform",
                                 jitter.unifrom.scale = "probability"))
  expect_error(res4 <- residuals(fit, type = "surrogate", jitter = "uniform",
                                 jitter.unifrom.scale = "probability", nsim = 10))
  expect_equal(length(res1), nrow(df1))
  expect_equal(length(res2), nrow(df1))
  expect_null(attr(res1, "draws")) # Since use residuals.glm, return a vector of numerical values
  expect_null(attr(res1, "draws_id"))
  expect_null(attr(res2, "draws"), "matrix")
  expect_null(attr(res2, "draws_id"), "matrix")

})

test_that("Surrogate approach: residuals.ord work for \"glm\" objects with ordinal response", {

  # Skips
  skip_on_cran()

  # Load data
  data(df1)

  # Fit cumulative link model
  fit <- glm(y ~ x + I(x ^ 2), data = df1, family = binomial)

  # Try to add "ord" class to the fitted model
  class(fit) <- c("ord", class(fit))

  # Compute residuals
  res1 <- residuals(fit)
  res2 <- residuals(fit, nsim = 10)
  res3 <- residuals(fit, type = "surrogate", jitter = "uniform",
                    jitter.unifrom.scale = "probability")
  res4 <- residuals(fit, type = "surrogate", jitter = "uniform",
                    jitter.unifrom.scale = "probability", nsim = 10)

  # Expectations
  expect_equal(length(res1), nrow(df1))
  expect_equal(length(res2), nrow(df1))
  expect_equal(length(res3), nrow(df1))
  expect_equal(length(res4), nrow(df1))
  expect_null(attr(res1, "draws"))
  expect_null(attr(res1, "draws_id"))
  expect_null(attr(res3, "draws"))
  expect_null(attr(res3, "draws_id"))
  expect_is(attr(res2, "draws"), "matrix")
  expect_is(attr(res2, "draws_id"), "matrix")
  expect_is(attr(res4, "draws"), "matrix")
  expect_is(attr(res4, "draws_id"), "matrix")
  expect_equal(dim(attr(res2, "draws")), c(nrow(df1), 10))
  expect_equal(dim(attr(res2, "draws_id")), c(nrow(df1), 10))
  expect_equal(dim(attr(res4, "draws")), c(nrow(df1), 10))
  expect_equal(dim(attr(res4, "draws_id")), c(nrow(df1), 10))

})

test_that("residuals work for \"lrm\" objects", {

  # Skips
  skip_on_cran()
  skip_if_not_installed("rms")

  # Load data
  data(df1)

  # Fit cumulative link model
  fit <- rms::lrm(y ~ x, data = df1)

  # Compute residuals
  res1 <- residuals(fit)
  res2 <- residuals(fit, nsim = 10)

  # Expectations
  expect_equal(length(res1), nrow(df1))
  expect_equal(length(res2), nrow(df1))
  expect_null(attr(res1, "draws"))
  expect_null(attr(res1, "draws_id"))
  expect_is(attr(res2, "draws"), "matrix")
  expect_is(attr(res2, "draws_id"), "matrix")
  expect_equal(dim(attr(res2, "draws")), c(nrow(df1), 10))
  expect_equal(dim(attr(res2, "draws_id")), c(nrow(df1), 10))

})

test_that("residuals work for \"orm\" objects", {

  # Skips
  skip_on_cran()
  skip_if_not_installed("rms")

  # Load data
  data(df1)

  # Fit cumulative link model
  fit <- rms::orm(y ~ x, data = df1, family = logistic)

  # Compute residuals
  res1 <- residuals(fit)
  res2 <- residuals(fit, nsim = 10)

  # Expectations
  expect_equal(length(res1), nrow(df1))
  expect_equal(length(res2), nrow(df1))
  expect_null(attr(res1, "draws"))
  expect_null(attr(res1, "draws_id"))
  expect_is(attr(res2, "draws"), "matrix")
  expect_is(attr(res2, "draws_id"), "matrix")
  expect_equal(dim(attr(res2, "draws")), c(nrow(df1), 10))
  expect_equal(dim(attr(res2, "draws_id")), c(nrow(df1), 10))

})

test_that("residuals work for \"polr\" objects", {

  # Skips
  skip_on_cran()
  skip_if_not_installed("MASS")

  # Load data
  data(df1)

  # Fit cumulative link model
  fit <- MASS::polr(y ~ x + I(x ^ 2), data = df1, method = "logistic")

  # Compute residuals
  res1 <- residuals(fit)
  res2 <- residuals(fit, nsim = 10)

  # Expectations
  expect_equal(length(res1), nrow(df1))
  expect_equal(length(res2), nrow(df1))
  expect_null(attr(res1, "draws"))
  expect_null(attr(res1, "draws_id"))
  expect_is(attr(res2, "draws"), "matrix")
  expect_is(attr(res2, "draws_id"), "matrix")
  expect_equal(dim(attr(res2, "draws")), c(nrow(df1), 10))
  expect_equal(dim(attr(res2, "draws_id")), c(nrow(df1), 10))

})

test_that("residualsAcat work for \"vglm\" objects", {

  # Skips
  skip_on_cran()
  skip_if_not_installed("rms")
  # library(PAsso)
  # Load data
  data(df1)

  # Fit cumulative link model
  suppressWarnings(
    fit <- VGAM::vglm(y ~ x + I(x ^ 2), data = df1,
                      family = VGAM::cumulative(link = "logit",
                                                parallel = TRUE))
  )

  resids_1 <- residualsAcat(object = fit,
                               type = "surrogate", jitter = "latent",
                               jitter.uniform.scale = "response",
                               nsim = 1)

  # head(temp_resids)

  suppressWarnings(
    fit2 <- VGAM::vglm(y ~ x + I(x ^ 2), data = df1,
                       family=acat(reverse=TRUE, parallel=TRUE))
  )
  resids_2 <- residualsAcat(object = fit2,
                           type = "surrogate", jitter = "latent",
                           jitter.uniform.scale = "response",
                           nsim = 30)

  fitted_temp <- do.call("vglm",
                         list(formula = PID ~ income.num + age + edu.year,
                              data = quote(ANES2016),
                              family = VGAM::acat(reverse = TRUE, parallel = TRUE)))

  temp_resids <- residualsAcat(object = fitted_temp,
                           type = "surrogate", jitter = "latent",
                           jitter.uniform.scale = "response",
                           nsim = 30)
  # dim(attr(temp_resids, "draws"))

  # Expectations
  expect_equal(length(resids_1), nrow(df1))
  expect_equal(length(resids_2), nrow(df1))
  expect_equal(length(temp_resids), nrow(ANES2016))

  expect_null(attr(resids_1, "draws"))
  expect_null(attr(resids_1, "draws_id"))

  expect_is(attr(resids_2, "draws"), "matrix")
  expect_is(attr(resids_2, "draws_id"), "matrix")

  expect_equal(dim(attr(temp_resids, "draws")), c(nrow(ANES2016), 30))
  expect_equal(dim(attr(temp_resids, "draws_id")), c(nrow(ANES2016), 30))
})

test_that("residuals work for \"vglm\" objects", {
  data(df1)

  # Fit cumulative link model
  suppressWarnings(
    fit <- VGAM::vglm(y ~ x + I(x ^ 2), data = df1,
                      family = VGAM::cumulative(link = "logit",
                                                parallel = TRUE))
  )
  # Compute residuals
  res1_1 <- residuals(fit)
  res1_2 <- residuals(fit, type = "surrogate", jitter = "latent",
                      jitter.uniform.scale = "probability")

  suppressWarnings(
    fit2 <- VGAM::vglm(y ~ x + I(x ^ 2), data = df1,
                       family=acat(reverse=TRUE, parallel=TRUE))
  )
  res2_1 <- residuals(fit2)
  res2_2 <- residuals(fit2, nsim=10)
  # summary(res2_1)

  # Expectations
  expect_equal(length(res1_1), nrow(df1))
  expect_equal(length(res1_1), nrow(df1))
  expect_equal(length(res2_1), nrow(df1))
  expect_equal(length(res2_2), nrow(df1))
  expect_null(attr(res1_1, "draws"))
  expect_null(attr(res1_1, "draws_id"))
  expect_null(attr(res1_2, "draws"))
  expect_null(attr(res1_2, "draws_id"))
  expect_null(attr(res2_1, "draws"))
  expect_is(attr(res2_2, "draws"), "matrix")
  expect_equal(dim(attr(res2_2, "draws")), c(nrow(df1), 10))
  expect_equal(dim(attr(res2_2, "draws_id")), c(nrow(df1), 10))

})

test_that("residuals work for \"clm\" objects with different link functions", {

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

  # Compute residuals
  res1 <- residuals(fit1)
  res2 <- residuals(fit2)
  res3 <- residuals(fit3)
  res4 <- residuals(fit4)
  res5 <- residuals(fit5)

  # Expectations
  expect_equal(length(res1), nrow(df1))
  expect_equal(length(res2), nrow(df1))
  expect_equal(length(res3), nrow(df1))
  expect_equal(length(res4), nrow(df1))
  expect_equal(length(res5), nrow(df1))

})
