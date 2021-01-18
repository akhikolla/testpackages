

test_that("Estimates lambda2 and omega are correct", {

  data(asrm, package = "Bayesrel")
  set.seed(1234)
  ee <- Bayesrel::strel(asrm, estimates = c("lambda2", "omega"), n.iter = 100, n.boot = 100, n.chains = 2)

  expect_equal(ee$Bayes$est$Bayes_lambda2, 0.7948049, tolerance = 1e-3)
  expect_equal(ee$freq$est$freq_lambda2, 0.7960336, tolerance = 1e-3)
  expect_equal(ee$Bayes$est$Bayes_omega, 0.7708523, tolerance = 1e-3)
  expect_equal(ee$freq$est$freq_omega, 0.7919616, tolerance = 1e-3)
  expect_equal(ee$Bayes$cred$low$Bayes_omega, 0.6719616, tolerance = 1e-3)
  if (as.numeric(R.Version()$major >= 4)) {
    expect_equal(as.numeric(ee$freq$conf$up$freq_lambda2), 0.865121, tolerance = 1e-3)
  } # because of the change in the RNG brought by the new R version

})



test_that("Bayes glb is correct", {

  data(asrm, package = "Bayesrel")
  set.seed(1234)
  ee <- Bayesrel::strel(asrm, estimates = "glb", n.iter = 100, freq = F, n.chains = 1)

  expect_equal(ee$Bayes$est$Bayes_glb, 0.8542316, tolerance = 1e-3)
  expect_equal(ee$Bayes$cred$up$Bayes_glb, 0.8950283, tolerance = 1e-3)


})


test_that("Bayes Alpha if item deleted is correct", {

  data(asrm, package = "Bayesrel")
  set.seed(1234)
  ee <- Bayesrel::strel(asrm, estimates = "alpha", n.iter = 100, freq = F, item.dropped = T, n.chains = 2)
  expect_equal(as.numeric(ee$Bayes$ifitem$est$alpha[1:2]), c(0.7207363, 0.7245768),
               tolerance = 1e-3)
  expect_equal(as.numeric(ee$Bayes$ifitem$cred$alpha[c(1, 10)]), c(0.6450673, 0.8049170),
               tolerance = 1e-3)

})


test_that("Freq omega with PFA is correct", {

  data(asrm, package = "Bayesrel")
  set.seed(1234)
  tt <- Bayesrel::strel(asrm, estimates = "omega", n.boot = 100, Bayes = F, omega.freq.method = "pfa")
  expect_equal(as.numeric(tt$freq$est$freq_omega), c(0.7966209), tolerance = 1e-3)
  if (as.numeric(R.Version()$major >= 4)) {
    expect_equal(as.numeric(tt$freq$conf$up$freq_omega), 0.8617495, tolerance = 1e-3)
  } # because of the change in the RNG brought by the new R version

})


test_that("Bayes prior and posterior probability for Alpha >.8 is correct", {

  data(asrm, package = "Bayesrel")
  set.seed(1234)
  tt <- Bayesrel::strel(asrm, estimates = "alpha", n.iter = 100, freq = F, n.chains = 2)
  ee <- Bayesrel::p_strel(tt, "alpha", .8)

  expect_equal(as.numeric(ee), c(0.1552618, 0.3300000), tolerance = 1e-3)

})


test_that("Omega results with missing data are correct", {

  data(asrm_mis, package = "Bayesrel")
  set.seed(1234)
  ee <- Bayesrel::strel(asrm_mis, estimates = c("omega"), n.iter = 100, n.chains = 2, n.boot = 100)
  expect_equal(as.numeric(ee$Bayes$cred$low$Bayes_omega), c(0.6900274),
               tolerance = 1e-3)
  expect_equal(as.numeric(ee$freq$est$freq_omega), c(0.7943602),
               tolerance = 1e-3)

})

test_that("Frequentist Lambda6 results with missing data and parametric bootstrap are correct", {

  data(asrm_mis, package = "Bayesrel")
  set.seed(1234)
  ee <- Bayesrel::strel(asrm_mis, estimates = c("lambda6"), Bayes = F, n.boot = 100, para.boot = T)
  expect_equal(as.numeric(ee$freq$est$freq_lambda6), c(0.7927271),
               tolerance = 1e-3)
  if (as.numeric(R.Version()$major >= 4)) {
    es <- extSoftVersion()
    blas   <- as.character(es["BLAS"])
    lapack <- La_library()
    if (grepl("atlas", tolower(blas)) || grepl("atlas", tolower(lapack)))
      tol <- 1e-2
    else
      tol <- 1e-3
    expect_equal(as.numeric(ee$freq$conf$low$freq_lambda6), 0.7188984, tolerance = tol)
  } # because of the change in the RNG brought by the new R version


})

test_that("Results with input cov matrix are correct", {

  data(asrm, package = "Bayesrel")
  cc <- cov(asrm)
  set.seed(1234)
  ee <- Bayesrel::strel(cov.mat = cc, estimates = c("lambda2"), n.iter = 100, n.chains = 2, n.boot = 100, n.obs = 500)
  expect_equal(as.numeric(ee$Bayes$cred$up$Bayes_lambda2), c(0.8215358),
               tolerance = 1e-3)
  expect_equal(as.numeric(ee$freq$est$freq_lambda2), c(0.7960336),
               tolerance = 1e-3)
  if (as.numeric(R.Version()$major >= 4)) {
    es <- extSoftVersion()
    blas   <- as.character(es["BLAS"])
    lapack <- La_library()
    if (grepl("atlas", tolower(blas)) || grepl("atlas", tolower(lapack)))
      tol <- 1e-2
    else
      tol <- 1e-3
    expect_equal(as.numeric(ee$freq$conf$low$freq_lambda2), 0.7724344, tolerance = tol)
  } # because of the change in the RNG brought by the new R version

})

