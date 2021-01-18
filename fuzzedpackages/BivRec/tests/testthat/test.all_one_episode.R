context("all subjects have only 1 episode")

library(BivRec)
n = round(runif(1, 80,100), digits = 0)
bivrec_data <- simBivRec(nsize=n, beta1=c(0.5,0.5), beta2=c(0,-0.5),
                tau_c=63, set=1.1)
firstepi <- which(bivrec_data$epi == 1)
bivrec_data2 = bivrec_data3 = bivrec_data[firstepi, -5]

#Idea of data where censoring happened on same day as event (rehospitalization)
#Fix for analysis with small amount (0.000000001) extra time was added to xij in a
#second episode with second episode (d1, d2) = (0,0).

for (i in unique(bivrec_data2$id)) {
  tmp_index <- which(bivrec_data2$id == i)
  if (bivrec_data2$d2[tmp_index]!=0) {
    tmp <- c(i, 2, 0.000000001, 0,
             0, 0, bivrec_data2$a1[tmp_index], bivrec_data2$a2[tmp_index])
    bivrec_data2 <- rbind(bivrec_data2, tmp)}
}

test_that("data check 2", {
  expect_equal(max(unique(bivrec_data2$epi)), 2)
  expect_equal(min(bivrec_data2$xij), 0.000000001)
})

test_that("lee check 2", {
  lee_reg <- bivrecReg(bivrecSurv(id, epi, xij, yij, d1, d2) ~ a1 + a2,
                       bivrec_data2, "Lee.et.al")
  expect_is(lee_reg, "bivrecReg")
  lee_coeffs <- coef.bivrecReg(lee_reg)
  expect_is(lee_coeffs, "matrix")
  expect_is(lee_coeffs[,1], "numeric")
})

test_that("np check 2", {
  npresult <- bivrecNP(response = with(bivrec_data2, bivrecSurv(id, epi, xij, yij, d1, d2)),
                       ai=1, u1 = seq(2, 15, 1), u2 = seq(1, 10, 1), conditional = TRUE,
                       given.interval = c(0, 10), level = 0.90)
  expect_is(npresult, "bivrecNP")
  expect_is(npresult$joint_cdf, "data.frame")
  expect_is(npresult$marginal_survival, "data.frame")
  expect_is(npresult$conditional_cdf$conditional, "data.frame")
})

#analysis for all d2=0 dataset
bivrec_data3$d2 <- rep(0, nrow(bivrec_data3))

test_that("data check 3", {
  expect_error(with(bivrec_data3, bivrecSurv(id, epi, xij, yij, d1, d2)))
})

test_that("lee check 3", {
  expect_error(bivrecReg(bivrecSurv(id, epi, xij, yij, d1, d2) ~ a1 + a2,
                       bivrec_data3, "Lee.et.al"))
})

test_that("chang check 3", {
  expect_error(bivrecReg(bivrecSurv(id, epi, xij, yij, d1, d2) ~ a1 + a2,
                         bivrec_data3, "Chang"))
})

test_that("np check 3", {
  expect_error(bivrecNP(response = with(bivrec_data3, bivrecSurv(id, epi, xij, yij, d1, d2)),
                       ai=1, u1 = seq(2, 15, 1), u2 = seq(1, 10, 1), conditional = TRUE,
                       given.interval = c(0, 10), level = 0.90))
})
