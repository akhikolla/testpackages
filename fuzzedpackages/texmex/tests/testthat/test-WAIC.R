context("WAIC")

test_that("Telling AIC whether to compute DIC or WAIC behaves as expected", {
  skip_on_cran()
  skip_on_travis()

  m <- evm(log(ALT.M / ALT.B), data = liver, qu=.7, method = "sim")
  aic <- AIC(m)

  expect_equal(length(aic), 3, label = "WAIC: AIC returns AIC, DIC and WAIC by default")

  aic <- AIC(m, DIC = FALSE)
  expect_equal(names(aic), c("AIC", "WAIC"), label="WAIC: AIC returns AIC and WAIC when DIC=FALSE")

  aic <- AIC(m, WAIC = FALSE)
  expect_equal(names(aic), c("AIC", "DIC"), label="WAIC: AIC returns AIC and DIC when WAIC=FALSE")

  aic <- AIC(m, DIC = FALSE, WAIC = FALSE)
  expect_equal(names(aic), "AIC", label = "WAIC: AIC returns only AIC when DIC=FALSE, WAIC=FALSE")
})

test_that("WAIC, AIC and DIC are broadly in alignment on well behaved data", {
  skip_on_cran()
  set.seed(1234)
  b <- function(){
    i <- sample(1:nrow(liver), size=nrow(liver), replace=TRUE)
    m <- evm(log(ALT.M / ALT.B), data = liver[i, ], qu=.7, eta = ~ as.numeric(dose),
             method = "sim", iter = 2000, thin = 1, family = cgpd)
    # Use cgpd to avoid xi < -0.5 whilst bootstrapping

    AIC(m)
  }

  aic <- t(replicate(100, b()))
  expect_gt(min(cor(aic)), .98, label="WAIC: AIC, DIC and WAIC are well correlated: cgpd liver")

  # We've analysed the liver data to death, so we "know" there is a dependence
  # on dose
  m0 <- evm(log(ALT.M / ALT.B), data = liver, qu=.7, method = "sim")
  m1 <- evm(log(ALT.M / ALT.B), data = liver, qu=.7, xi = ~ as.numeric(dose), method = "sim")

  expect_lt(AIC(m1)["WAIC"], AIC(m0)["WAIC"],
            label = "WAIC: 'correct' model selected for liver data: gpd")

  b <- function(){
    i <- sample(1:nrow(portpirie), size=nrow(portpirie), replace=TRUE)
    m <- evm(SeaLevel, data = portpirie[i, ], family = gev,
             method = "sim", iter = 2000, thin = 1)
    AIC(m)
  }

  aic <- t(replicate(100, b()))
  expect_gt(min(cor(aic)), .98, label="WAIC: AIC, DIC and WAIC are well correlated: gev portpirie")
})
