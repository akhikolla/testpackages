context("PAsso: plot3D")

test_that("plot3D works for \"PAsso\" objects", {

  # Skips
  skip_on_cran()
  skip_if_not_installed("copula")
  skip_if_not_installed("plotly")

  library(copula)
  library(plotly)
  # Load data
  data("ANES2016")

  # multivariate analysis (2 variables) --------------------------------------------------------------------
  skip("Take times and no display")
  PAsso_2v <- PAsso(responses = c("PreVote.num", "PID"),
                    adjustments = c("income.num", "age", "edu.year"),
                    data = ANES2016, uni.model = "logit",
                    method = c("kendall"),
                    resids.type = "surrogate", jitter = "latent")

  testPlots <- plot3D(PAsso_2v, y1="PreVote.num", y2="PID")

  # Expectations
  expect_s3_class(testPlots, "plotly")

})

