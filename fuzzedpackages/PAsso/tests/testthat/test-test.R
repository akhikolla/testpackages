context("PAsso: 'test()' for PAsso object")


test_that("test() for 'PAsso' object without parallel", {

  # Skips
  skip_on_cran()
  skip_if_not_installed("MASS")
  skip_if_not_installed("tidyverse")
  skip_if_not_installed("progress")

  # Load data
  data("ANES2016")

  PAsso_2 <- PAsso(responses = c("PreVote.num", "PID", "selfLR"),
                   adjustments = c("income.num", "age", "edu.year"),
                   data = ANES2016, uni.model = "probit",
                   method = c("kendall"),
                   resids.type = "surrogate", jitter = "latent")

  # test.PAsso function: Conduct inference fo object of "PAsso" class ----------------------------
  system.time(PAsso_2_test <- test(object = PAsso_2, bootstrap_rep=100, H0=0, parallel=F))

  # Expectations
  expect_s3_class(PAsso_2_test, "PAsso.test")
  expect_equal(length(attr(PAsso_2, "arguments")), 5)

})
