context("plot_cusum")
# 
# cs <- cusum(failure_probability = 0.1,
#             patient_outcomes = rbinom(100,1,0.1),
#             limit = 2,
#             odds_multiplier = 2,
#             reset = FALSE)
# cs <- cs[,-1]
# 
# test_that("plot asks for t",
#           expect_error(plot(cs))
# )
# 
# 
# 
# test_that("plot asks for ct", {
#   cs <- cusum(failure_probability = 0.1,
#               patient_outcomes = rbinom(100,1,0.1),
#               limit = 2,
#               odds_multiplier = 2,
#               reset = FALSE)
#   cs <- cs[,-3]
#   expect_error(plot(cs))
#   }
#   )

set.seed(2313)
cs <- cusum(failure_probability = 0.1,
            patient_outcomes = rbinom(100,1,0.1),
            limit = 2,
            odds_multiplier = 2,
            reset = FALSE)

if(requireNamespace("vdiffr")){
  disp_cs <- function() plot(cs)
  disp_new <- function()  plot(cs$t, cs$ct, type = "l")
  vdiffr::expect_doppelganger("CUSUM plot", disp_cs)
  vdiffr::expect_doppelganger("Test plot", disp_new)
}
