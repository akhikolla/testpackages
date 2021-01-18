context("PAsso: plot")

test_that("plot works for \"PAsso\" objects", {

  # Skips
  skip_on_cran()

  # Load data
  data("ANES2016")

  # multivariate analysis (5 variables) --------------------------------------------------------------------
  PAsso_5v <- PAsso(responses = c("PreVote.num", "PID", "selfLR", "ClinLR", "TrumpLR"),
                    adjustments = c("income.num", "age", "edu.year"),
                    data = ANES2016, uni.model = "logit",
                    method = c("kendall"),
                    resids.type = "surrogate", jitter = "latent")

  p2 <- plot(PAsso_5v)
  # Expectations
  expect_s3_class(p2, "ggmatrix")

})

