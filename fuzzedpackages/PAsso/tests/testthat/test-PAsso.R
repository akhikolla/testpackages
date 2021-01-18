context("PAsso: Partial association")


test_that("Advanced PAsso for two responses", {

  # Skips
  skip_on_cran()
  skip_if_not_installed("MASS")
  skip_if_not_installed("tidyverse")

  # Load data
  data("ANES2016")

  # "PAsso" advanced using of the function: The First way (Advanced), input a few models directly ------------------------------
  # Test "PAsso" function: Partial Association by surrogate residuals regression models
  fit.vote<- glm(PreVote.num ~ income.num+ age + edu.year, data = ANES2016,
                 family = binomial(link = "probit"))
  fit.PID<- polr(as.factor(PID) ~ income.num+age+edu.year, data = ANES2016,
                 method="probit", Hess = TRUE)

  PAsso_adv1 <- PAsso(fitted.models=list(fit.vote, fit.PID),
                                  method = c("kendall"),
                                  resids.type = "surrogate")

  # Test jittering
  PAsso_adv1_jit <- PAsso(fitted.models=list(fit.vote, fit.PID),
                                      association = c("partial"),
                                      method = c("kendall"),
                          jitter = "uniform",
                          jitter.uniform.scale = "response")

  # Expectations
  expect_s3_class(PAsso_adv1, "PAsso")
  expect_equal(dim(PAsso_adv1$corr), c(2,2))
  expect_equal(length(attr(PAsso_adv1, "arguments")), 5)


})

test_that("Simple PAsso() for two responses",{
  skip_if_not_installed("MASS")
  skip_if_not_installed("tidyverse")

  PAsso_1 <- PAsso(responses = c("PreVote.num", "PID"),
                   adjustments = c("income.num", "age", "edu.year"),
                   data = ANES2016
                   # association = c("partial"),
                   # models = c("probit", "probit"),
                   # method = c("kendall"),
                   # resids.type = "surrogate", fitted.models = NULL,
                   # rep_num = 30
  )

  p1 <- plot(PAsso_1)

  # Expectations
  expect_s3_class(PAsso_1, "PAsso")
  expect_equal(dim(PAsso_1$corr), c(2,2))
  expect_equal(length(attr(PAsso_1, "arguments")), 5)
  expect_is(p1, "ggmatrix")
  expect_equal(dim(PAsso_1$corr)[1], 2)

  # Expectations
  expect_error(PAsso(responses = c("PreVote.num", "PID", "selfLR"),
                     adjustments = c("income.num", "age", "edu.year"),
                     data = ANES2016, uni.model = "probit",
                     method = c("kendall"),
                     resids.type = "surrogate", jitter = "uniform"),
               "only supported for logit-type models")

  expect_error(PAsso(responses = c("PreVote.num", "PID", "selfLR"),
                     adjustments = c("income.num", "age", "edu.year"),
                     data = ANES2016, uni.model = "probit",
                     method = c("kendall"),
                     resids.type = "surrogate", jitter = "uniform",
                     jitter.uniform.scale = "response"),
               NA)

  expect_error(PAsso(responses = c("PreVote.num", "PID", "selfLR"),
                     adjustments = c("income.num", "age", "edu.year"),
                     data = ANES2016, uni.model = "logit",
                     method = c("kendall"),
                     resids.type = "surrogate", jitter = "uniform",
                     jitter.uniform.scale = "prob"),
               NA)

})

test_that("Simple PAsso() for two responses(Choose one binary response to test residuals.ord)",{
  skip_if_not_installed("MASS")
  skip_if_not_installed("tidyverse")

  PAsso_1 <- PAsso(responses = c("PreVote.num", "PID"),
                   adjustments = c("income.num", "age", "edu.year"),
                   data = ANES2016
                   # association = c("partial"),
                   # models = c("probit", "probit"),
                   # method = c("kendall"),
                   # resids.type = "surrogate", fitted.models = NULL,
                   # rep_num = 30
  )

  p1 <- plot(PAsso_1)

  # Expectations: PAsso object


  # Expectations: PAsso object
  expect_s3_class(PAsso_1, "PAsso")
  expect_s3_class(PAsso_1$fitted.models[[1]], "ord")
  expect_s3_class(PAsso_1$fitted.models[[2]], "polr")
  expect_equal(dim(PAsso_1$corr), c(2,2))
  expect_equal(length(attr(PAsso_1, "arguments")), 5)
  expect_equal(dim(PAsso_1$corr)[1], 2)

  # Expectations: plot
  expect_is(p1, "ggmatrix")


  # Expectations: test compatibility
  expect_error(PAsso(responses = c("PreVote.num", "PID", "selfLR"),
                     adjustments = c("income.num", "age", "edu.year"),
                     data = ANES2016, uni.model = "probit",
                     method = c("kendall"),
                     resids.type = "surrogate", jitter = "uniform"),
               "only supported for logit-type models")

  expect_error(PAsso(responses = c("PreVote.num", "PID", "selfLR"),
                     adjustments = c("income.num", "age", "edu.year"),
                     data = ANES2016, uni.model = "probit",
                     method = c("kendall"),
                     resids.type = "surrogate", jitter = "uniform",
                     jitter.uniform.scale = "response"),
               NA)

  expect_error(PAsso(responses = c("PreVote.num", "PID", "selfLR"),
                     adjustments = c("income.num", "age", "edu.year"),
                     data = ANES2016, uni.model = "logit",
                     method = c("kendall"),
                     resids.type = "surrogate", jitter = "uniform",
                     jitter.uniform.scale = "prob"),
               NA)

})

test_that("Check links: Simple PAsso() for two responses",{
  skip_if_not_installed("MASS")
  skip_if_not_installed("tidyverse")

  # logit models!
  PAsso_1 <- PAsso(responses = c("PreVote.num", "PID"),
                   adjustments = c("income.num", "age", "edu.year"),
                   data = ANES2016, uni.model = "logit")
  links <- attr(PAsso_1, "models")

  # Expectations
  expect_equal(links, rep("logit", 2))
  expect_equal(PAsso_1$fitted.models[[1]]$family$link, "logit")
  expect_equal(PAsso_1$fitted.models[[2]]$method, "logistic")

  # probit models!
  PAsso_1_probit <- PAsso(responses = c("PreVote.num", "PID"),
                   adjustments = c("income.num", "age", "edu.year"),
                   data = ANES2016, uni.model = "probit")
  links <- attr(PAsso_1_probit, "models")

  # Expectations
  expect_equal(links, rep("probit", 2))
  expect_equal(PAsso_1_probit$fitted.models[[1]]$family$link, "probit")
  expect_equal(PAsso_1_probit$fitted.models[[2]]$method, "probit")

  # Complicated combination: logit + probit
  PAsso_1_com <- PAsso(responses = c("PreVote.num", "PID"),
                       adjustments = c("income.num", "age", "edu.year"),
                       data = ANES2016, uni.model = "probit",
                       models = c("logit", "probit"))

  links <- attr(PAsso_1_com, "models")

  # Expectations
  expect_equal(links, c("logit", "probit"))
  expect_equal(PAsso_1_com$fitted.models[[1]]$family$link, "logit")
  expect_equal(PAsso_1_com$fitted.models[[2]]$method, "probit")

})

test_that("Check acat: Simple PAsso() for two responses",{
  skip_if_not_installed("MASS")
  skip_if_not_installed("tidyverse")

  # Load data
  data("ANES2016")
  library(PAsso)
  library(VGAM)
  # acat models!

  PAsso_1 <- PAsso(responses = c("PreVote.num", "PID"),
                   adjustments = c("income.num", "age", "edu.year"),
                   data = ANES2016,
                   models = c("probit", "acat"))
  links <- attr(PAsso_1, "models")
  # links
  # Expectations
  expect_equal(links, c("probit", "acat"))
  expect_equal(PAsso_1$fitted.models[[1]]$family$link, "probit")
  expect_equal(PAsso_1$fitted.models[[2]]@family@vfamily[1], "acat")

  # probit models!
  PAsso_1_probit <- PAsso(responses = c("PreVote.num", "PID"),
                          adjustments = c("income.num", "age", "edu.year"),
                          data = ANES2016, uni.model = "probit")
  links <- attr(PAsso_1_probit, "models")

  # Expectations
  expect_equal(links, rep("probit", 2))
  expect_equal(PAsso_1_probit$fitted.models[[1]]$family$link, "probit")
  expect_equal(PAsso_1_probit$fitted.models[[2]]$method, "probit")

  # Complicated combination: logit + probit
  PAsso_1_com <- PAsso(responses = c("PreVote.num", "PID"),
                       adjustments = c("income.num", "age", "edu.year"),
                       data = ANES2016, uni.model = "probit",
                       models = c("logit", "probit"))

  links <- attr(PAsso_1_com, "models")

  # Expectations
  expect_equal(links, c("logit", "probit"))
  expect_equal(PAsso_1_com$fitted.models[[1]]$family$link, "logit")
  expect_equal(PAsso_1_com$fitted.models[[2]]$method, "probit")

  # Complicated combination: acat + acat + logit
  PAsso_2_com <- PAsso(responses = c("PID", "selfLR", "PreVote.num"),
                       adjustments = c("income.num", "age", "edu.year"),
                       data = ANES2016, uni.model = "probit",
                       models = c("acat", "acat", "logit"))

  links <- attr(PAsso_2_com, "models")

  # Expectations
  expect_equal(links, c("acat", "acat", "logit"))
  expect_equal(PAsso_2_com$fitted.models[[1]]@family@vfamily[1], "acat")
  expect_equal(PAsso_2_com$fitted.models[[2]]@family@vfamily[1], "acat")
  expect_equal(PAsso_2_com$fitted.models[[3]]$family$link, "logit")

  # Complicated combination: acat + acat + acat
  PAsso_3_com <- PAsso(responses = c("PID", "selfLR", "PreVote.num"),
                       adjustments = c("income.num", "age", "edu.year"),
                       data = ANES2016, uni.model = "acat")
  # summary(PAsso_3_com)

  links <- attr(PAsso_3_com, "models")

  # Expectations
  expect_equal(links, c("acat", "acat", "acat"))
  expect_equal(PAsso_3_com$fitted.models[[1]]@family@vfamily[1], "acat")
  expect_equal(PAsso_3_com$fitted.models[[2]]@family@vfamily[1], "acat")
  expect_equal(PAsso_3_com$fitted.models[[3]]@family@vfamily[1], "acat")
})

test_that("Check continuous respones: Simple PAsso() for two responses",{
  skip_if_not_installed("MASS")
  skip_if_not_installed("tidyverse")

  # Load data
  data("ANES2016")
  library(PAsso)
  library(VGAM)
  # acat models!

  PAsso_1 <- PAsso(responses = c("PreVote.num", "PID"),
                   adjustments = c("income.num", "age", "edu.year"),
                   data = ANES2016,
                   models = c("probit", "acat"))
  links <- attr(PAsso_1, "models")
  # links
  # Expectations
  expect_equal(links, c("probit", "acat"))
  expect_equal(PAsso_1$fitted.models[[1]]$family$link, "probit")
  expect_equal(PAsso_1$fitted.models[[2]]@family@vfamily[1], "acat")

})

