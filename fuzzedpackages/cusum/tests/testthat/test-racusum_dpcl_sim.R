context("racusum_dpcl_sim")

risks <- c(0.001, 0.01, 0.1, 0.002, 0.02, 0.2)

set.seed(2046)
patient_risks <- sample(x = risks, size = 10, replace = TRUE)

test_that("Output of dpcl sim", {
          observe <- round(dpcl <- racusum_limit_dpcl(
            patient_risks = patient_risks,
            N = 1000,
            odds_multiplier = 2,
            alpha = 0.05,
            seed = 32423
            ),4)
          
          expect <- c(0,0,0,0,0,0,0.5108,0.5108,0.8393,0.5098)
          expect_equal(observe, expect)
}          
)
