test_that("basic run", {
  testthat::skip_on_cran()
  local_edition(3)

  set.seed(48)
  ## complex test after wrong calculation results
  ## RLumCarlo
  TL110 <- expect_s3_class(RLumCarlo::run_MC_TL_DELOC(
    s = 5e12, E = 0.97, R = 5e-10, times = seq(20,400,2),
    N_e = 78, method = "seq"),
    "RLumCarlo_Model_Output")


  TL230 <- expect_s3_class(RLumCarlo::run_MC_TL_DELOC(
    s = 5e14, E = 1.55, R = 5e-10, times = seq(20,400,2),
    N_e = 52, method = "seq"),
    "RLumCarlo_Model_Output")

  TL325 <- expect_s3_class(RLumCarlo::run_MC_TL_DELOC(
    s = 5e13, E = 1.7, R = 5e-10, times = seq(20,400,2),
    N_e = 1909, method = "seq"),
    "RLumCarlo_Model_Output")

  ## combine
  object <- expect_s3_class(c(TL110, TL230, TL325), "RLumCarlo_Model_Output")

  ##test result
  expect_equal(round(sum(object$signal), digits = 0),10141)

})
