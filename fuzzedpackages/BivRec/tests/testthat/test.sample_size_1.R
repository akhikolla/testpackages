context("Sample size 1")

library(BivRec)
bivrec_data <- simBivRec(nsize=1, beta1=c(0.5,0.5), beta2=c(0,-0.5),
                tau_c=63, set=1.1)

test_that("data check 1", {
  expect_is(bivrec_data, "data.frame")
  expect_equal(unique(bivrec_data$id), 1)
})

check_np <- function() {
  if (max(bivrec_data$epi)==1) {
    skip("np check")
  }
}

check_reg <- function() {
  if (max(bivrec_data$epi)==1) {
    expect_error(bivrecReg(bivrecSurv(id, epi, xij, yij, d1, d2) ~ a1 + a2,
                           bivrec_data, "Lee.et.al"))
    expect_error(bivrecReg(bivrecSurv(id, epi, xij, yij, d1, d2) ~ a1 + a2,
                           bivrec_data, "Chang"))
  }
}

### Lee and Chang methods for 1 subject with several episodes may work
####   or may lead to singular systems / no convergence

test_that("np check", {
  npresult <- bivrecNP(response = with(bivrec_data, bivrecSurv(id, epi, xij, yij, d1, d2)),
                       ai=1, u1 = seq(2, 15, 1), u2 = seq(1, 10, 1), conditional = FALSE)
  expect_is(npresult, "bivrecNP")
  expect_is(npresult$joint_cdf, "data.frame")
  expect_is(npresult$marginal_survival, "data.frame")
  expect_is(npresult$conditional_cdf, "NULL")
#note that conditional SE and CI's are bootstrap based so cannot run conditional for a sample of 1
})
