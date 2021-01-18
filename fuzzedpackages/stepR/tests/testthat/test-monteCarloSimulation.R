context("monteCarloSimulation")

testrnorm1 <- function(testn, testr, testseed, testfamily, testintervalSystem, testlengths, testpenalty,
                      tolerance, allLengths = testlengths, save = TRUE) {
  retvector <- monteCarloSimulation(n = testn, r = testr, seed = testseed, family = testfamily,
                                    intervalSystem = testintervalSystem, lengths = NULL,
                                    penalty = NULL, output = "vector")
  retmax <- monteCarloSimulation(n = testn, r = testr, seed = testseed, family = testfamily,
                                 intervalSystem = testintervalSystem, lengths = testlengths,
                                 penalty = testpenalty, output = "maximum")
  set.seed(testseed)
  for (i in 1:testr) {
    y <- rnorm(testn)
    expect_equal(computeStat(y, family = testfamily, intervalSystem = testintervalSystem,
                             lengths = NULL, penalty = "none", output = "vector", sd = 1),
                 retvector[ ,i], tolerance = tolerance)
    expect_equal(computeStat(y, family = testfamily, intervalSystem = testintervalSystem,
                             lengths = testlengths, penalty = testpenalty, output = "maximum", sd = 1),
                 retmax[i], tolerance = tolerance)
  }
  testkeyList <- list(testn, digest::sha1(list("gauss"), digits = 6), testintervalSystem)
  testkey <- digest::digest(testkeyList)
  expect_identical(attr(retvector, "keyList"), testkeyList)
  expect_identical(attr(retvector, "key"), testkey)
  expect_identical(attr(retvector, "n"), testn)
  expect_identical(attr(retvector, "lengths"), allLengths)
  expect_identical(attr(retvector, "save"), save)
  expect_s3_class(retvector, "MCSimulationVector")
  expect_identical(names(attributes(retvector)), c("dim", "class", "keyList", "key", "n", "lengths", "save"))

  testkeyList <- list(testn, digest::sha1(list("gauss"), digits = 6),
                      testintervalSystem, as.integer(testlengths), testpenalty)
  testkey <- digest::digest(testkeyList)
  expect_identical(attr(retmax, "keyList"), testkeyList)
  expect_identical(attr(retmax, "key"), testkey)
  expect_identical(attr(retmax, "n"), testn)
  expect_identical(attr(retmax, "lengths"), testlengths)
  expect_identical(attr(retmax, "save"), save)
  expect_s3_class(retmax, "MCSimulationMaximum")
  expect_identical(names(attributes(retmax)), c("class", "keyList", "key", "n", "lengths", "save"))
}

testrnorm <- function(testn, testr, testseed, testfamily, testintervalSystem, testlengths, testpenalty,
                       tolerance, allLengths = testlengths, save = TRUE) {
  retvector <- monteCarloSimulation(n = testn, r = testr, seed = testseed, family = testfamily,
                                    intervalSystem = testintervalSystem, lengths = NULL,
                                    penalty = NULL, output = "vector")
  retmax <- monteCarloSimulation(n = testn, r = testr, seed = testseed, family = testfamily,
                                 intervalSystem = testintervalSystem, lengths = testlengths,
                                 penalty = testpenalty, output = "maximum")
  set.seed(testseed)
  for (i in 1:testr) {
    y <- rnorm(testn)
    expect_equal(computeStat(y, family = testfamily, intervalSystem = testintervalSystem,
                             lengths = NULL, penalty = "none", output = "vector"),
                 retvector[ ,i], tolerance = tolerance)
    expect_equal(computeStat(y, family = testfamily, intervalSystem = testintervalSystem,
                             lengths = testlengths, penalty = testpenalty, output = "maximum"),
                 retmax[i], tolerance = tolerance)
  }

  testkeyList <- list(testn, digest::sha1(list("hsmuce"), digits = 6),
                      testintervalSystem)
  testkey <- digest::digest(testkeyList)
  expect_identical(attr(retvector, "keyList"), testkeyList)
  expect_identical(attr(retvector, "key"), testkey)
  expect_identical(attr(retvector, "n"), testn)
  expect_identical(attr(retvector, "lengths"), allLengths)
  expect_identical(attr(retvector, "save"), save)
  expect_s3_class(retvector, "MCSimulationVector")
  expect_identical(names(attributes(retvector)), c("dim", "class", "keyList", "key", "n", "lengths", "save"))

  testkeyList <- list(testn, digest::sha1(list("hsmuce"), digits = 6),
                      testintervalSystem, as.integer(testlengths), testpenalty)
  testkey <- digest::digest(testkeyList)
  expect_identical(attr(retmax, "keyList"), testkeyList)
  expect_identical(attr(retmax, "key"), testkey)
  expect_identical(attr(retmax, "n"), testn)
  expect_identical(attr(retmax, "lengths"), testlengths)
  expect_identical(attr(retmax, "save"), save)
  expect_s3_class(retmax, "MCSimulationMaximum")
  expect_identical(names(attributes(retmax)), c("class", "keyList", "key", "n", "lengths", "save"))
}

testrand.genMDependent <- function(data) {
  ma <- testcomputeMA(data$covariances)
  as.numeric(arima.sim(n = data$n, list(ma = ma), sd = 1)) / sqrt(sum(c(1, ma)^2))
}

testrand.genJsmurf <- function(data) {
  ma <- testcomputeMA(data$filter$acf)
  as.numeric(arima.sim(n = data$n, list(ma = ma), sd = 1)) / sqrt(sum(c(1, ma)^2))
}

testcomputeMA <- function(cov) {
  N <- length(cov) - 1
  alpha <- matrix(rep(0, N * N), nrow = N, ncol = N)
  vs <- rep(cov[1], N + 1)

  for(i in 1:N){
    alpha[i, 1:i] <- rep(0, i)
    alpha[i, i] <- cov[i + 1] / vs[1]
    if(i > 1) {
      for(k in 1:(i - 1)) {
        js <- 0:(k - 1)
        alpha[i, i - k] <- (cov[i - k + 1] - sum(alpha[i, i - js] * alpha[k, k - js] * vs[js + 1])) / vs[k + 1]
      }
    }
    js <- 0:(i - 1)
    vs[i + 1] <- vs[i + 1] - sum(alpha[i, i - js]^2 * vs[js + 1])
  }

  alpha[N, ]
}

testbyfilter <- function(testn, testr, testseed, testfamily, testintervalSystem, testlengths, testpenalty,
                         tolerance, testcovariances, allLengths = testlengths, save = TRUE) {
  retvector <- monteCarloSimulation(n = testn, r = testr, seed = testseed, family = testfamily,
                                    intervalSystem = testintervalSystem, lengths = NULL,
                                    penalty = NULL, output = "vector", covariances = testcovariances)
  retmax <- monteCarloSimulation(n = testn, r = testr, seed = testseed, family = testfamily,
                                 intervalSystem = testintervalSystem, lengths = testlengths,
                                 penalty = testpenalty, output = "maximum", covariances = testcovariances)
  set.seed(testseed)
  ma <- testcomputeMA(testcovariances)
  for (i in 1:testr) {
    y <- as.numeric(arima.sim(n = testn, list(ma = ma), sd = 1)) / sqrt(sum(c(1, ma)^2))

    expect_equal(computeStat(y, family = testfamily, intervalSystem = testintervalSystem,
                             lengths = NULL, penalty = "none", output = "vector",
                             covariances = testcovariances),
                 retvector[ ,i], tolerance = tolerance)
    expect_equal(computeStat(y, family = testfamily, intervalSystem = testintervalSystem,
                             lengths = testlengths, penalty = testpenalty, output = "maximum",
                             covariances = testcovariances),
                 retmax[i], tolerance = tolerance)
  }

  testkeyList <- list(testn, digest::sha1(list("mDependentPS", ma), digits = 6),
                      testintervalSystem)
  testkey <- digest::digest(testkeyList)
  expect_identical(attr(retvector, "keyList"), testkeyList)
  expect_identical(attr(retvector, "key"), testkey)
  expect_identical(attr(retvector, "n"), testn)
  expect_identical(attr(retvector, "lengths"), allLengths)
  expect_identical(attr(retvector, "save"), save)
  expect_s3_class(retvector, "MCSimulationVector")
  expect_identical(names(attributes(retvector)), c("dim", "class", "keyList", "key", "n", "lengths", "save"))

  testkeyList <- list(testn, digest::sha1(list("mDependentPS", ma), digits = 6),
                      testintervalSystem, as.integer(testlengths), testpenalty)
  testkey <- digest::digest(testkeyList)
  expect_identical(attr(retmax, "keyList"), testkeyList)
  expect_identical(attr(retmax, "key"), testkey)
  expect_identical(attr(retmax, "n"), testn)
  expect_identical(attr(retmax, "lengths"), testlengths)
  expect_identical(attr(retmax, "save"), save)
  expect_s3_class(retmax, "MCSimulationMaximum")
  expect_identical(names(attributes(retmax)), c("class", "keyList", "key", "n", "lengths", "save"))
}

testJsmurf <- function(testn, testr, testseed, testfamily, testintervalSystem, testlengths, testpenalty,
                       tolerance, testfilter, allLengths = testlengths, save = TRUE) {
  retvector <- monteCarloSimulation(n = testn, r = testr, seed = testseed, family = testfamily,
                                    intervalSystem = testintervalSystem, lengths = NULL,
                                    penalty = NULL, output = "vector", filter = testfilter)
  retmax <- monteCarloSimulation(n = testn, r = testr, seed = testseed, family = testfamily,
                                 intervalSystem = testintervalSystem, lengths = testlengths,
                                 penalty = testpenalty, output = "maximum", filter = testfilter)
  set.seed(testseed)
  ma <- testcomputeMA(testfilter$acf)
  for (i in 1:testr) {
    y <- as.numeric(arima.sim(n = testn, list(ma = ma), sd = 1)) / sqrt(sum(c(1, ma)^2))

    expect_equal(computeStat(y, family = testfamily, intervalSystem = testintervalSystem,
                             lengths = NULL, penalty = "none", output = "vector",
                             filter = testfilter, sd = 1),
                 retvector[ ,i], tolerance = tolerance)
    expect_equal(computeStat(y, family = testfamily, intervalSystem = testintervalSystem,
                             lengths = testlengths, penalty = testpenalty, output = "maximum",
                             filter = testfilter, sd = 1),
                 retmax[i], tolerance = tolerance)
  }

  testkeyList <- list(testn, digest::sha1(list("jsmurf", ma), digits = 6),
                      testintervalSystem)
  testkey <- digest::digest(testkeyList)
  expect_identical(attr(retvector, "keyList"), testkeyList)
  expect_identical(attr(retvector, "key"), testkey)
  expect_identical(attr(retvector, "n"), testn)
  expect_identical(attr(retvector, "lengths"), allLengths)
  expect_identical(attr(retvector, "save"), save)
  expect_s3_class(retvector, "MCSimulationVector")
  expect_identical(names(attributes(retvector)), c("dim", "class", "keyList", "key", "n", "lengths", "save"))

  testkeyList <- list(testn, digest::sha1(list("jsmurf", ma), digits = 6),
                      testintervalSystem, as.integer(testlengths), testpenalty)
  testkey <- digest::digest(testkeyList)
  expect_identical(attr(retmax, "keyList"), testkeyList)
  expect_identical(attr(retmax, "key"), testkey)
  expect_identical(attr(retmax, "n"), testn)
  expect_identical(attr(retmax, "lengths"), testlengths)
  expect_identical(attr(retmax, "save"), save)
  expect_s3_class(retmax, "MCSimulationMaximum")
  expect_identical(names(attributes(retmax)), c("class", "keyList", "key", "n", "lengths", "save"))
}

testJsmurfPS <- function(testn, testr, testseed, testfamily, testintervalSystem, testlengths, testpenalty,
                         tolerance, testfilter, allLengths = testlengths, save = TRUE) {
  retvector <- monteCarloSimulation(n = testn, r = testr, seed = testseed, family = testfamily,
                                    intervalSystem = testintervalSystem, lengths = NULL,
                                    penalty = NULL, output = "vector", filter = testfilter)
  retmax <- monteCarloSimulation(n = testn, r = testr, seed = testseed, family = testfamily,
                                 intervalSystem = testintervalSystem, lengths = testlengths,
                                 penalty = testpenalty, output = "maximum", filter = testfilter)
  set.seed(testseed)
  ma <- testcomputeMA(testfilter$acf)
  for (i in 1:testr) {
    y <- as.numeric(arima.sim(n = testn, list(ma = ma), sd = 1)) / sqrt(sum(c(1, ma)^2))

    expect_equal(computeStat(y, family = testfamily, intervalSystem = testintervalSystem,
                             lengths = NULL, penalty = "none", output = "vector",
                             filter = testfilter, sd = 1),
                 retvector[ ,i], tolerance = tolerance)
    expect_equal(computeStat(y, family = testfamily, intervalSystem = testintervalSystem,
                             lengths = testlengths, penalty = testpenalty, output = "maximum",
                             filter = testfilter, sd = 1),
                 retmax[i], tolerance = tolerance)
  }

  testkeyList <- list(testn, digest::sha1(list("jsmurfPS", ma), digits = 6),
                      testintervalSystem)
  testkey <- digest::digest(testkeyList)
  expect_identical(attr(retvector, "keyList"), testkeyList)
  expect_identical(attr(retvector, "key"), testkey)
  expect_identical(attr(retvector, "n"), testn)
  expect_identical(attr(retvector, "lengths"), allLengths)
  expect_identical(attr(retvector, "save"), save)
  expect_s3_class(retvector, "MCSimulationVector")
  expect_identical(names(attributes(retvector)), c("dim", "class", "keyList", "key", "n", "lengths", "save"))

  testkeyList <- list(testn, digest::sha1(list("jsmurfPS", ma), digits = 6),
                      testintervalSystem, as.integer(testlengths), testpenalty)
  testkey <- digest::digest(testkeyList)
  expect_identical(attr(retmax, "keyList"), testkeyList)
  expect_identical(attr(retmax, "key"), testkey)
  expect_identical(attr(retmax, "n"), testn)
  expect_identical(attr(retmax, "lengths"), testlengths)
  expect_identical(attr(retmax, "save"), save)
  expect_s3_class(retmax, "MCSimulationMaximum")
  expect_identical(names(attributes(retmax)), c("class", "keyList", "key", "n", "lengths", "save"))
}

testJsmurfLR <- function(testn, testr, testseed, testfamily, testintervalSystem, testlengths, testpenalty,
                         tolerance, testfilter, allLengths = testlengths, save = TRUE) {
  retvector <- monteCarloSimulation(n = testn, r = testr, seed = testseed, family = testfamily,
                                    intervalSystem = testintervalSystem, lengths = NULL,
                                    penalty = NULL, output = "vector", filter = testfilter)
  retmax <- monteCarloSimulation(n = testn, r = testr, seed = testseed, family = testfamily,
                                 intervalSystem = testintervalSystem, lengths = testlengths,
                                 penalty = testpenalty, output = "maximum", filter = testfilter)
  set.seed(testseed)
  ma <- testcomputeMA(testfilter$acf)
  for (i in 1:testr) {
    y <- as.numeric(arima.sim(n = testn, list(ma = ma), sd = 1)) / sqrt(sum(c(1, ma)^2))

    expect_equal(computeStat(y, family = testfamily, intervalSystem = testintervalSystem,
                             lengths = NULL, penalty = "none", output = "vector",
                             filter = testfilter, sd = 1),
                 retvector[ ,i], tolerance = tolerance)
    expect_equal(computeStat(y, family = testfamily, intervalSystem = testintervalSystem,
                             lengths = testlengths, penalty = testpenalty, output = "maximum",
                             filter = testfilter, sd = 1),
                 retmax[i], tolerance = tolerance)
  }

  testkeyList <- list(testn, digest::sha1(list("jsmurfLR", ma), digits = 6),
                      testintervalSystem)
  testkey <- digest::digest(testkeyList)
  expect_identical(attr(retvector, "keyList"), testkeyList)
  expect_identical(attr(retvector, "key"), testkey)
  expect_identical(attr(retvector, "n"), testn)
  expect_identical(attr(retvector, "lengths"), allLengths)
  expect_identical(attr(retvector, "save"), save)
  expect_s3_class(retvector, "MCSimulationVector")
  expect_identical(names(attributes(retvector)), c("dim", "class", "keyList", "key", "n", "lengths", "save"))

  testkeyList <- list(testn, digest::sha1(list("jsmurfLR", ma), digits = 6),
                      testintervalSystem, as.integer(testlengths), testpenalty)
  testkey <- digest::digest(testkeyList)
  expect_identical(attr(retmax, "keyList"), testkeyList)
  expect_identical(attr(retmax, "key"), testkey)
  expect_identical(attr(retmax, "n"), testn)
  expect_identical(attr(retmax, "lengths"), testlengths)
  expect_identical(attr(retmax, "save"), save)
  expect_s3_class(retmax, "MCSimulationMaximum")
  expect_identical(names(attributes(retmax)), c("class", "keyList", "key", "n", "lengths", "save"))
}

testHjsmurf <- function(testn, testr, testseed, testfamily, testintervalSystem, testlengths, testpenalty,
                        tolerance, testfilter, allLengths = testlengths, save = TRUE) {
  retvector <- monteCarloSimulation(n = testn, r = testr, seed = testseed, family = testfamily,
                                    intervalSystem = testintervalSystem, lengths = NULL,
                                    penalty = NULL, output = "vector", filter = testfilter)
  retmax <- monteCarloSimulation(n = testn, r = testr, seed = testseed, family = testfamily,
                                 intervalSystem = testintervalSystem, lengths = testlengths,
                                 penalty = testpenalty, output = "maximum", filter = testfilter)
  set.seed(testseed)
  ma <- testcomputeMA(testfilter$acf)
  for (i in 1:testr) {
    y <- as.numeric(arima.sim(n = testn, list(ma = ma), sd = 1)) / sqrt(sum(c(1, ma)^2))

    expect_equal(computeStat(y, family = testfamily, intervalSystem = testintervalSystem,
                             lengths = NULL, penalty = "none", output = "vector",
                             filter = testfilter),
                 retvector[ ,i], tolerance = tolerance)
    expect_equal(computeStat(y, family = testfamily, intervalSystem = testintervalSystem,
                             lengths = testlengths, penalty = testpenalty, output = "maximum",
                             filter = testfilter),
                 retmax[i], tolerance = tolerance)
  }

  testkeyList <- list(testn, digest::sha1(list("hjsmurf", ma), digits = 6),
                      testintervalSystem)
  testkey <- digest::digest(testkeyList)
  expect_identical(attr(retvector, "keyList"), testkeyList)
  expect_identical(attr(retvector, "key"), testkey)
  expect_identical(attr(retvector, "n"), testn)
  expect_identical(attr(retvector, "lengths"), allLengths)
  expect_identical(attr(retvector, "save"), save)
  expect_s3_class(retvector, "MCSimulationVector")
  expect_identical(names(attributes(retvector)), c("dim", "class", "keyList", "key", "n", "lengths", "save"))

  testkeyList <- list(testn, digest::sha1(list("hjsmurf", ma), digits = 6),
                      testintervalSystem, as.integer(testlengths), testpenalty)
  testkey <- digest::digest(testkeyList)
  expect_identical(attr(retmax, "keyList"), testkeyList)
  expect_identical(attr(retmax, "key"), testkey)
  expect_identical(attr(retmax, "n"), testn)
  expect_identical(attr(retmax, "lengths"), testlengths)
  expect_identical(attr(retmax, "save"), save)
  expect_s3_class(retmax, "MCSimulationMaximum")
  expect_identical(names(attributes(retmax)), c("class", "keyList", "key", "n", "lengths", "save"))
}

testHjsmurfSPS <- function(testn, testr, testseed, testfamily, testintervalSystem, testlengths, testpenalty,
                           tolerance, testfilter, allLengths = testlengths, save = TRUE) {
  retvector <- monteCarloSimulation(n = testn, r = testr, seed = testseed, family = testfamily,
                                    intervalSystem = testintervalSystem, lengths = NULL,
                                    penalty = NULL, output = "vector", filter = testfilter)
  retmax <- monteCarloSimulation(n = testn, r = testr, seed = testseed, family = testfamily,
                                 intervalSystem = testintervalSystem, lengths = testlengths,
                                 penalty = testpenalty, output = "maximum", filter = testfilter)
  set.seed(testseed)
  ma <- testcomputeMA(testfilter$acf)
  for (i in 1:testr) {
    y <- as.numeric(arima.sim(n = testn, list(ma = ma), sd = 1)) / sqrt(sum(c(1, ma)^2))

    expect_equal(computeStat(y, family = testfamily, intervalSystem = testintervalSystem,
                             lengths = NULL, penalty = "none", output = "vector",
                             filter = testfilter),
                 retvector[ ,i], tolerance = tolerance)
    expect_equal(computeStat(y, family = testfamily, intervalSystem = testintervalSystem,
                             lengths = testlengths, penalty = testpenalty, output = "maximum",
                             filter = testfilter),
                 retmax[i], tolerance = tolerance)
  }

  testkeyList <- list(testn, digest::sha1(list("hjsmurfSPS", ma), digits = 6),
                      testintervalSystem)
  testkey <- digest::digest(testkeyList)
  expect_identical(attr(retvector, "keyList"), testkeyList)
  expect_identical(attr(retvector, "key"), testkey)
  expect_identical(attr(retvector, "n"), testn)
  expect_identical(attr(retvector, "lengths"), allLengths)
  expect_identical(attr(retvector, "save"), save)
  expect_s3_class(retvector, "MCSimulationVector")
  expect_identical(names(attributes(retvector)), c("dim", "class", "keyList", "key", "n", "lengths", "save"))

  testkeyList <- list(testn, digest::sha1(list("hjsmurfSPS", ma), digits = 6),
                      testintervalSystem, as.integer(testlengths), testpenalty)
  testkey <- digest::digest(testkeyList)
  expect_identical(attr(retmax, "keyList"), testkeyList)
  expect_identical(attr(retmax, "key"), testkey)
  expect_identical(attr(retmax, "n"), testn)
  expect_identical(attr(retmax, "lengths"), testlengths)
  expect_identical(attr(retmax, "save"), save)
  expect_s3_class(retmax, "MCSimulationMaximum")
  expect_identical(names(attributes(retmax)), c("class", "keyList", "key", "n", "lengths", "save"))
}

testHjsmurfLR <- function(testn, testr, testseed, testfamily, testintervalSystem, testlengths, testpenalty,
                          tolerance, testfilter, allLengths = testlengths, save = TRUE) {
  retvector <- monteCarloSimulation(n = testn, r = testr, seed = testseed, family = testfamily,
                                    intervalSystem = testintervalSystem, lengths = NULL,
                                    penalty = NULL, output = "vector", filter = testfilter)
  retmax <- monteCarloSimulation(n = testn, r = testr, seed = testseed, family = testfamily,
                                 intervalSystem = testintervalSystem, lengths = testlengths,
                                 penalty = testpenalty, output = "maximum", filter = testfilter)
  set.seed(testseed)
  ma <- testcomputeMA(testfilter$acf)
  for (i in 1:testr) {
    y <- as.numeric(arima.sim(n = testn, list(ma = ma), sd = 1)) / sqrt(sum(c(1, ma)^2))

    expect_equal(computeStat(y, family = testfamily, intervalSystem = testintervalSystem,
                             lengths = NULL, penalty = "none", output = "vector",
                             filter = testfilter),
                 retvector[ ,i], tolerance = tolerance)
    expect_equal(computeStat(y, family = testfamily, intervalSystem = testintervalSystem,
                             lengths = testlengths, penalty = testpenalty, output = "maximum",
                             filter = testfilter),
                 retmax[i], tolerance = tolerance)
  }

  testkeyList <- list(testn, digest::sha1(list("hjsmurfLR", ma), digits = 6),
                      testintervalSystem)
  testkey <- digest::digest(testkeyList)
  expect_identical(attr(retvector, "keyList"), testkeyList)
  expect_identical(attr(retvector, "key"), testkey)
  expect_identical(attr(retvector, "n"), testn)
  expect_identical(attr(retvector, "lengths"), allLengths)
  expect_identical(attr(retvector, "save"), save)
  expect_s3_class(retvector, "MCSimulationVector")
  expect_identical(names(attributes(retvector)), c("dim", "class", "keyList", "key", "n", "lengths", "save"))

  testkeyList <- list(testn, digest::sha1(list("hjsmurfLR", ma), digits = 6),
                      testintervalSystem, as.integer(testlengths), testpenalty)
  testkey <- digest::digest(testkeyList)
  expect_identical(attr(retmax, "keyList"), testkeyList)
  expect_identical(attr(retmax, "key"), testkey)
  expect_identical(attr(retmax, "n"), testn)
  expect_identical(attr(retmax, "lengths"), testlengths)
  expect_identical(attr(retmax, "save"), save)
  expect_s3_class(retmax, "MCSimulationMaximum")
  expect_identical(names(attributes(retmax)), c("class", "keyList", "key", "n", "lengths", "save"))
}

testLR <- function(testn, testr, testseed, testfamily, tolerance, testfilter, allLengths, save = TRUE, ...) {
  ret <- monteCarloSimulation(n = testn, r = testr, seed = testseed, family = testfamily,
                              intervalSystem = "all", output = "vector", filter = testfilter, ...)
  set.seed(testseed)
  ma <- testcomputeMA(testfilter$acf)
  for (i in 1:testr) {
    y <- as.numeric(arima.sim(n = testn, list(ma = ma), sd = 1)) / sqrt(sum(c(1, ma)^2))

    testfit <- stepblock(0, 0, testn / testfilter$sr, 0)
    expect_equal(computeStat(y, family = testfamily, penalty = "none", output = "vector",
                             filter = testfilter, fit = testfit, ...),
                 ret[ ,i], tolerance = tolerance)

  }

  testkeyList <- list(testn,
                      digest::sha1(list("LR", testfilter$type, unname(testfilter$param),
                                        testfilter$sr, testfilter$len), digits = 6),
                      "all")
  testkey <- digest::digest(testkeyList)
  expect_identical(attr(ret, "keyList"), testkeyList)
  expect_identical(attr(ret, "key"), testkey)
  expect_identical(attr(ret, "n"), testn)
  expect_identical(attr(ret, "lengths"), allLengths)
  expect_identical(attr(ret, "save"), save)
  expect_s3_class(ret, "MCSimulationVector")
  expect_identical(names(attributes(ret)), c("dim", "class", "keyList", "key", "n", "lengths", "save"))
}

test2Param <- function(testn, testr, testseed, testfamily, tolerance, testfilter, allLengths, save = TRUE, ...) {
  ret <- monteCarloSimulation(n = testn, r = testr, seed = testseed, family = testfamily,
                              intervalSystem = "all", output = "vector", filter = testfilter, ...)
  set.seed(testseed)
  ma <- testcomputeMA(testfilter$acf)
  for (i in 1:testr) {
    y <- as.numeric(arima.sim(n = testn, list(ma = ma), sd = 1)) / sqrt(sum(c(1, ma)^2))

    testfit <- stepblock(0, 0, testn / testfilter$sr, 0)
    expect_equal(computeStat(y, family = testfamily, penalty = "none", output = "vector",
                             filter = testfilter, fit = testfit, ...),
                 ret[ ,i], tolerance = tolerance)

  }

  testkeyList <- list(testn,
                      digest::sha1(list("2Param", testfilter$type, unname(testfilter$param),
                                        testfilter$sr, testfilter$len), digits = 6),
                      "all")
  testkey <- digest::digest(testkeyList)
  expect_identical(attr(ret, "keyList"), testkeyList)
  expect_identical(attr(ret, "key"), testkey)
  expect_identical(attr(ret, "n"), testn)
  expect_identical(attr(ret, "lengths"), allLengths)
  expect_identical(attr(ret, "save"), save)
  expect_s3_class(ret, "MCSimulationVector")
  expect_identical(names(attributes(ret)), c("dim", "class", "keyList", "key", "n", "lengths", "save"))
}

test_that("argument n is tested and works", {
  expect_error(monteCarloSimulation())
  expect_error(monteCarloSimulation(n = "s"))
  expect_error(monteCarloSimulation(n = c(1L, 2L)))
  expect_error(monteCarloSimulation(n = NA))
  expect_error(monteCarloSimulation(n = Inf))
  expect_error(monteCarloSimulation(n = 0L))
  expect_error(monteCarloSimulation(n = -1L))
  expect_identical(monteCarloSimulation(n = 32), monteCarloSimulation(n = 32L))
  expect_identical(monteCarloSimulation(n = 32.5), monteCarloSimulation(n = 32L))
  expect_equal(monteCarloSimulation(n = 32L, r = 10L),
                   monteCarloSimulation(n = 32L, r = 10L, seed = 32L,
                                        rand.gen = function(data) rnorm(data$n),
                                        family = "gauss", intervalSystem = "all", lengths = NULL,
                                        penalty = "none", output = "vector", sd = 1), check.attributes = FALSE)

  testrnorm1(testn = 32L, testr = 10L, testseed = 1L, testfamily = "gauss", testintervalSystem = "all",
            testlengths = 1:32, testpenalty = "none", tolerance = 1e-14)
})

test_that("argument r is tested and works", {
  expect_error(monteCarloSimulation(n = 32L, r = "s"))
  expect_error(monteCarloSimulation(n = 32L, r = c(1L, 2L)))
  expect_error(monteCarloSimulation(n = 32L, r = NA))
  expect_error(monteCarloSimulation(n = 32L, r = Inf))
  expect_error(monteCarloSimulation(n = 32L, r = 0L))
  expect_error(monteCarloSimulation(n = 32L, r = -1L))
  expect_identical(monteCarloSimulation(n = 32L), monteCarloSimulation(n = 32L, r = 1e4L))
  expect_identical(monteCarloSimulation(n = 32L, r = 100), monteCarloSimulation(n = 32L, r = 100L))
  expect_identical(monteCarloSimulation(n = 32L, r = 100.5), monteCarloSimulation(n = 32L, r = 100L))
})

test_that("argument seed is tested and works", {
  expect_error(suppressWarnings(monteCarloSimulation(n = 32L, seed = "s")))
  expect_error(suppressWarnings(monteCarloSimulation(n = 32L, seed = NA)))
  expect_error(suppressWarnings(monteCarloSimulation(n = 32L, seed = Inf)))
  set.seed(0L)
  expect_identical(monteCarloSimulation(n = 32L, r = 100L, seed = "no"),
                   monteCarloSimulation(n = 32L, r = 100L, seed = 0L))
  expect_identical(monteCarloSimulation(n = 32L, r = 100L),
                   monteCarloSimulation(n = 32L, r = 100L, seed = 32L))
  expect_identical(monteCarloSimulation(n = 32L, r = 100L, seed = 0),
                   monteCarloSimulation(n = 32L, r = 100L, seed = 0L))
  expect_identical(monteCarloSimulation(n = 32L, r = 100L, seed = 0.5),
                   monteCarloSimulation(n = 32L, r = 100L, seed = 0L))

  testrnorm1(testn = 32L, testr = 10L, testseed = -1L, testfamily = "gauss", testintervalSystem = "all",
            testlengths = 1:32, testpenalty = "none", tolerance = 1e-14)
})

test_that("argument rand.gen is tested and works", {
  expect_error(monteCarloSimulation(n = 32L, rand.gen = "s"))
  expect_error(monteCarloSimulation(n = 32L, rand.gen = function() {}))
  expect_error(monteCarloSimulation(n = 32L, rand.gen = function(test) {rnorm(32)}))
  expect_error(monteCarloSimulation(n = 32L, rand.gen = function(data, test) {rnorm(32)}))

  expect_error(monteCarloSimulation(n = 32L, rand.gen = function(data) {rnorm(31)}))
  expect_error(monteCarloSimulation(n = 32L, rand.gen = function(data) {c(rnorm(31), NA)}))
  expect_error(monteCarloSimulation(n = 32L, rand.gen = function(data) {c(rnorm(31), "NA")}))
  expect_error(monteCarloSimulation(n = 32L, rand.gen = function(data) {c(rnorm(31), Inf)}))

  ret <- monteCarloSimulation(n = 32L, r = 100L, rand.gen = function(data) {rnorm(32)})
  expect_equal(monteCarloSimulation(n = 32L, r = 100L), ret, check.attributes = FALSE)
  expect_false(attr(ret, "save"))
})

test_that("argument messages is tested and works", {
  expect_error(monteCarloSimulation(n = 32L, messages = "s"))
  expect_error(monteCarloSimulation(n = 32L, messages = c(1L, 2L)))
  expect_error(monteCarloSimulation(n = 32L, messages = NA))
  expect_error(monteCarloSimulation(n = 32L, messages = Inf))
  expect_error(monteCarloSimulation(n = 32L, messages = 0L))
  expect_error(monteCarloSimulation(n = 32L, messages = -1L))
  expect_identical(suppressMessages(monteCarloSimulation(n = 32L, r = 100L, messages = 10L)),
                   monteCarloSimulation(n = 32L, r = 100L))
  expect_identical(suppressMessages(monteCarloSimulation(n = 32L, r = 100L, messages = 3)),
                   suppressMessages(monteCarloSimulation(n = 32L, r = 100L, messages = 3.5)))
  expect_identical(suppressMessages(monteCarloSimulation(n = 32L, r = 100L, messages = 10L)),
                   monteCarloSimulation(n = 32L, r = 100L))
})

test_that("argument family is tested and works in the default case", { # other families below
  expect_error(monteCarloSimulation(n = 32L, family = ""))
  expect_error(monteCarloSimulation(n = 32L, family = c("gauss", "hsmuce")))
  expect_identical(monteCarloSimulation(n = 32L, r = 10L),
                   monteCarloSimulation(n = 32L, r = 10L, family = "gauss"))
})

test_that("argument intervalSystem is tested and works", {
  expect_error(monteCarloSimulation(n = 32L, intervalSystem = ""))
  expect_error(monteCarloSimulation(n = 32L, intervalSystem = "dya"))
  expect_error(monteCarloSimulation(n = 32L, intervalSystem = "dyalen"))
  expect_error(monteCarloSimulation(n = 32L, intervalSystem = "dyapar"))
  expect_identical(monteCarloSimulation(n = 34L, r = 10L),
                   monteCarloSimulation(n = 34L, r = 10L, intervalSystem = "all"))

  testrnorm1(testn = 28L, testr = 10L, testseed = 10L, testfamily = "gauss", testintervalSystem = "dyaLen",
            testlengths = as.integer(2^(0:4)), testpenalty = "sqrt", tolerance = 1e-14)
  testrnorm1(testn = 33L, testr = 10L, testseed = 10L, testfamily = "gauss", testintervalSystem = "dyaPar",
            testlengths = as.integer(2^(0:5)), testpenalty = "sqrt", tolerance = 1e-14)
})

test_that("argument lengths is tested and works", {
  expect_warning(ret <- monteCarloSimulation(n = 32L, r = 10L, lengths = 1:10))
  expect_identical(ret, monteCarloSimulation(n = 32L, r = 10L))

  expect_error(monteCarloSimulation(n = 32L, output = "maximum", lengths = "s"))
  expect_error(monteCarloSimulation(n = 32L, output = "maximum", lengths = c(1:10, NA)))
  expect_error(monteCarloSimulation(n = 32L, output = "maximum", lengths = c(1:10, Inf)))
  expect_error(monteCarloSimulation(n = 32L, output = "maximum", lengths = 0:10))
  expect_error(monteCarloSimulation(n = 32L, output = "maximum", lengths = -1L))
  expect_error(monteCarloSimulation(n = 32L, output = "maximum", lengths = c(1L, 46L)))
  expect_error(monteCarloSimulation(n = 32L, output = "maximum", intervalSystem = "dyaLen", lengths = 3L))
  expect_error(monteCarloSimulation(n = 32L, output = "maximum", intervalSystem = "dyaLen", lengths = 64L))
  expect_error(monteCarloSimulation(n = 32L, output = "maximum", intervalSystem = "dyaPar", lengths = 3L))
  expect_error(monteCarloSimulation(n = 32L, output = "maximum", intervalSystem = "dyaPar", lengths = 64L))

  expect_identical(monteCarloSimulation(n = 32L, output = "maximum", r = 10L, lengths = c(1, 4, 10)),
                   monteCarloSimulation(n = 32L, output = "maximum", r = 10L, lengths = c(1L, 4L, 10L)))
  expect_warning(ret <- monteCarloSimulation(n = 32L, output = "maximum", r = 10L, lengths = c(1:10, 10)))
  expect_identical(ret, monteCarloSimulation(n = 32L, output = "maximum", r = 10L, lengths = c(1:10)))
  expect_identical(monteCarloSimulation(n = 32L, output = "maximum", r = 10L, lengths = c(10:1)),
                   monteCarloSimulation(n = 32L, output = "maximum", r = 10L, lengths = c(1:10)))
  expect_identical(monteCarloSimulation(n = 32L, output = "maximum", r = 10L, lengths = c(1:10 + 0.5)),
                   monteCarloSimulation(n = 32L, output = "maximum", r = 10L, lengths = c(1:10)))

  expect_equal(monteCarloSimulation(n = 32L, output = "maximum", r = 10L, lengths = 2^(0:5)),
               monteCarloSimulation(n = 32L, output = "maximum", r = 10L, intervalSystem = "dyaLen"),
               check.attributes = FALSE)
  expect_equal(monteCarloSimulation(n = 32L, output = "maximum", r = 10L, lengths = 2^(3:4)),
               monteCarloSimulation(n = 32L, output = "maximum", r = 10L, intervalSystem = "dyaLen",
                                    lengths = 2^(3:4)), check.attributes = FALSE)

  testlengths <- c(1L, 3:5, 10:23)
  testrnorm1(testn = 28L, testr = 10L, testseed = 10L, testfamily = "gauss", testintervalSystem = "all",
            testlengths = testlengths, testpenalty = "sqrt", tolerance = 1e-14, allLengths = 1:28)

  testlengths <- c(1L, 4L, 32L)
  testrnorm1(testn = 32L, testr = 10L, testseed = 10L, testfamily = "gauss", testintervalSystem = "dyaLen",
            testlengths = testlengths, testpenalty = "sqrt", tolerance = 1e-14, allLengths = as.integer(2^(0:5)))

  testlengths <- c(2L, 4L, 16L)
  testrnorm1(testn = 24L, testr = 10L, testseed = 10L, testfamily = "gauss", testintervalSystem = "dyaPar",
            testlengths = testlengths, testpenalty = "sqrt", tolerance = 1e-14, allLengths = as.integer(2^(0:4)))

  expect_warning(ret <- monteCarloSimulation(10L, r = 10L, lengths = integer(0), output = "maximum"))
  compare <- rep(-Inf, 10)
  class(compare) <- "MCSimulationMaximum"
  expect_equal(ret, compare, check.attributes = FALSE)
})

test_that("argument penalty is tested and works in the default case", {
  expect_warning(ret <- monteCarloSimulation(n = 32L, r = 10L, penalty = "log"))
  expect_identical(ret, monteCarloSimulation(n = 32L, r = 10L))

  expect_error(monteCarloSimulation(n = 32L, output = "maximum", penalty = ""))
  expect_error(monteCarloSimulation(n = 32L, output = "maximum", penalty = c("gauss")))
  expect_identical(monteCarloSimulation(n = 32L, r = 10L, sd = 1),
                   monteCarloSimulation(n = 32L, r = 10L, penalty = "none", sd = 1))

  testlengths <- c(1L, 3:5, 10:23)
  testrnorm1(testn = 28L, testr = 10L, testseed = 10L, testfamily = "gauss", testintervalSystem = "all",
            testlengths = testlengths, testpenalty = "sqrt", tolerance = 1e-14, allLengths = 1:28)

  testlengths <- c(1L, 4L, 32L)
  testrnorm1(testn = 32L, testr = 10L, testseed = 10L, testfamily = "gauss", testintervalSystem = "dyaLen",
            testlengths = testlengths, testpenalty = "sqrt", tolerance = 1e-14, allLengths = as.integer(2^(0:5)))

  testlengths <- c(2L, 4L, 16L)
  testrnorm1(testn = 24L, testr = 10L, testseed = 10L, testfamily = "gauss", testintervalSystem = "dyaPar",
            testlengths = testlengths, testpenalty = "sqrt", tolerance = 1e-14, allLengths = as.integer(2^(0:4)))

  testlengths <- c(1L, 3:5, 10:23)
  testrnorm1(testn = 28L, testr = 10L, testseed = 10L, testfamily = "gauss", testintervalSystem = "all",
            testlengths = testlengths, testpenalty = "log", tolerance = 1e-14, allLengths = 1:28)

  testlengths <- c(1L, 4L, 32L)
  testrnorm1(testn = 32L, testr = 10L, testseed = 10L, testfamily = "gauss", testintervalSystem = "dyaLen",
            testlengths = testlengths, testpenalty = "log", tolerance = 1e-14, allLengths = as.integer(2^(0:5)))

  testlengths <- c(2L, 4L, 16L)
  testrnorm1(testn = 24L, testr = 10L, testseed = 10L, testfamily = "gauss", testintervalSystem = "dyaPar",
            testlengths = testlengths, testpenalty = "log", tolerance = 1e-14, allLengths = as.integer(2^(0:4)))
})

test_that("argument output is tested and works in the default case", {
  expect_error(monteCarloSimulation(n = 32L, output = 1))
  expect_error(monteCarloSimulation(n = 32L, output = "mu"))
  expect_identical(monteCarloSimulation(n = 32L, r = 10L),
                   monteCarloSimulation(n = 32L, r = 10L, output = "vector"))
})

test_that("... is tested and works", {
  expect_error(monteCarloSimulation(n = 32L, familty = "hsmuce"))
  expect_error(monteCarloSimulation(n = 32L, sd = "s"))
  expect_identical(monteCarloSimulation(n = 32L, r = 1e2),
                   monteCarloSimulation(n = 32L, r = 1e2, sd = 5))
})

test_that("family 'hsmuce' works", {
  testthat::skip_on_cran()
  expect_error(monteCarloSimulation(family = "hsmuce"))
  expect_error(monteCarloSimulation(family = "hsmuce", n = "s"))
  expect_identical(monteCarloSimulation(family = "hsmuce", n = 32),
                   monteCarloSimulation(family = "hsmuce", n = 32L))
  expect_identical(monteCarloSimulation(family = "hsmuce", n = 32.5),
                   monteCarloSimulation(family = "hsmuce", n = 32L))
  expect_equal(monteCarloSimulation(family = "hsmuce", n = 32L, r = 10L),
               monteCarloSimulation(n = 32L, r = 10L, seed = 32L,
                                    rand.gen = function(data) rnorm(data$n),
                                    family = "hsmuce", intervalSystem = "dyaPar", output = "vector"),
               check.attributes = FALSE)

  testrnorm(testn = 32L, testr = 10L, testseed = 1L, testfamily = "hsmuce", testintervalSystem = "all",
            testlengths = 2:32, testpenalty = "none", tolerance = 1e-14)

  testrnorm(testn = 32L, testr = 10L, testseed = -1L, testfamily = "hsmuce", testintervalSystem = "dyaPar",
            testlengths = as.integer(2^(1:5)), testpenalty = "none", tolerance = 1e-14)

  expect_error(monteCarloSimulation(family = "hsmuce", n = 32L, rand.gen = "s"))
  expect_error(monteCarloSimulation(family = "hsmuce", n = 32L, rand.gen = function() {}))
  expect_error(monteCarloSimulation(family = "hsmuce", n = 32L, rand.gen = function(test) {rnorm(32)}))
  expect_error(monteCarloSimulation(family = "hsmuce", n = 32L,
                                    rand.gen = function(data, test) {rnorm(32)}))

  expect_error(monteCarloSimulation(family = "hsmuce", n = 32L, rand.gen = function(data) {rnorm(31)}))
  expect_error(monteCarloSimulation(family = "hsmuce", n = 32L,
                                    rand.gen = function(data) {c(rnorm(31), NA)}))
  expect_error(monteCarloSimulation(family = "hsmuce", n = 32L,
                                    rand.gen = function(data) {c(rnorm(31), "NA")}))
  expect_error(monteCarloSimulation(family = "hsmuce", n = 32L,
                                    rand.gen = function(data) {c(rnorm(31), Inf)}))

  expect_equal(monteCarloSimulation(family = "hsmuce", n = 32L, r = 100L),
               monteCarloSimulation(family = "hsmuce", n = 32L, r = 100L,
                                    rand.gen = function(data) {rnorm(32)}), check.attributes = FALSE)

  expect_identical(monteCarloSimulation(family = "hsmuce", n = 34L, r = 10L),
                   monteCarloSimulation(family = "hsmuce", n = 34L, r = 10L, intervalSystem = "dyaPar"))

  testrnorm(testn = 28L, testr = 10L, testseed = 10L, testfamily = "hsmuce", testintervalSystem = "all",
            testlengths = 2:28, testpenalty = "none", tolerance = 1e-14)
  testrnorm(testn = 28L, testr = 10L, testseed = 10L, testfamily = "hsmuce", testintervalSystem = "dyaLen",
            testlengths = as.integer(2^(1:4)), testpenalty = "none", tolerance = 1e-14)

  expect_error(monteCarloSimulation(family = "hsmuce", output = "maximum", n = 32L, lengths = 1L))
  expect_error(monteCarloSimulation(family = "hsmuce", output = "maximum", n = 32L, lengths = c(33L)))
  expect_error(monteCarloSimulation(family = "hsmuce", output = "maximum", n = 32L, intervalSystem = "all",
                                    lengths = 1L))
  expect_error(monteCarloSimulation(family = "hsmuce", output = "maximum", n = 32L, intervalSystem = "dyaLen",
                                    lengths = 1L))

  expect_equal(monteCarloSimulation(family = "hsmuce", output = "maximum", n = 32L, r = 10L,
                                    intervalSystem = "all", lengths = 2^(1:5)),
               monteCarloSimulation(family = "hsmuce", output = "maximum", n = 32L, r = 10L,
                                    intervalSystem = "dyaLen"), check.attributes = FALSE)
  expect_equal(monteCarloSimulation(family = "hsmuce", output = "maximum", n = 32L, r = 10L,
                                    intervalSystem = "all", lengths = 2^(3:4)),
               monteCarloSimulation(family = "hsmuce", output = "maximum", n = 32L, r = 10L,
                                    intervalSystem = "dyaLen", lengths = 2^(3:4)), check.attributes = FALSE)

  testlengths <- c(2:5, 10:23)
  testrnorm(testn = 28L, testr = 10L, testseed = 10L, testfamily = "hsmuce", testintervalSystem = "all",
            testlengths = testlengths, testpenalty = "none", tolerance = 1e-14, allLengths = 2:28)

  testlengths <- c(2L, 4L, 8L)
  testrnorm(testn = 32L, testr = 10L, testseed = 10L, testfamily = "hsmuce", testintervalSystem = "dyaLen",
            testlengths = testlengths, testpenalty = "none", tolerance = 1e-14, allLengths = as.integer(2^(1:5)))

  testlengths <- c(2L, 4L, 16L)
  testrnorm(testn = 24L, testr = 10L, testseed = 10L, testfamily = "hsmuce", testintervalSystem = "dyaPar",
            testlengths = testlengths, testpenalty = "none", tolerance = 1e-14, allLengths = as.integer(2^(1:4)))

  expect_identical(monteCarloSimulation(family = "hsmuce", n = 32L, r = 10L),
                   monteCarloSimulation(family = "hsmuce", n = 32L, r = 10L, penalty = "none"))
  expect_identical(monteCarloSimulation(family = "hsmuce", output = "maximum", n = 32L, r = 10L),
                   monteCarloSimulation(family = "hsmuce", output = "maximum", n = 32L, r = 10L, penalty = "none"))

  testlengths <- c(2:5, 8L, 10:23)
  testrnorm(testn = 28L, testr = 10L, testseed = 10L, testfamily = "hsmuce", testintervalSystem = "all",
            testlengths = testlengths, testpenalty = "sqrt", tolerance = 1e-14, allLengths = 2:28)

  testlengths <- c(2L, 4L, 32L)
  testrnorm(testn = 32L, testr = 10L, testseed = 10L, testfamily = "hsmuce", testintervalSystem = "dyaLen",
            testlengths = testlengths, testpenalty = "sqrt", tolerance = 1e-14, allLengths = as.integer(2^(1:5)))

  testlengths <- c(2L, 4L, 16L)
  testrnorm(testn = 24L, testr = 10L, testseed = 10L, testfamily = "hsmuce", testintervalSystem = "dyaPar",
            testlengths = testlengths, testpenalty = "sqrt", tolerance = 1e-14, allLengths = as.integer(2^(1:4)))

  testlengths <- c(2:5, 8L, 10:23)
  testrnorm(testn = 28L, testr = 10L, testseed = 10L, testfamily = "hsmuce", testintervalSystem = "all",
            testlengths = testlengths, testpenalty = "log", tolerance = 1e-14, allLengths = 2:28)

  testlengths <- c(2L, 4L, 32L)
  testrnorm(testn = 32L, testr = 10L, testseed = 10L, testfamily = "hsmuce", testintervalSystem = "dyaLen",
            testlengths = testlengths, testpenalty = "log", tolerance = 1e-14, allLengths = as.integer(2^(1:5)))

  testlengths <- c(2L, 4L, 16L)
  testrnorm(testn = 24L, testr = 10L, testseed = 10L, testfamily = "hsmuce", testintervalSystem = "dyaPar",
            testlengths = testlengths, testpenalty = "log", tolerance = 1e-14, allLengths = as.integer(2^(1:4)))

  expect_error(monteCarloSimulation(family = "hsmuce", n = 32L, intervalsystem = "dyaPar"))
  expect_error(monteCarloSimulation(family = "hsmuce", n = 32L, sd = 1))

  ret <- monteCarloSimulation(1L, family = "hsmuce", r = 10L)
  compare <- matrix(ncol = 10, nrow = 0)
  class(compare) <- "MCSimulationVector"
  expect_equal(ret, compare, check.attributes = FALSE)
})

test_that("family 'mDependentPS' works", {
  testthat::skip_on_cran()
  testcovariances <- as.numeric(ARMAacf(ar = c(), ma = c(0.7, 0.4, 0.3, 0.1), lag.max = 4))
  expect_error(monteCarloSimulation(family = "mDependentPS", covariances = testcovariances))
  expect_error(monteCarloSimulation(family = "mDependentPS", n = "s", covariances = testcovariances))
  expect_identical(monteCarloSimulation(family = "mDependentPS", n = 32, covariances = testcovariances),
                   monteCarloSimulation(family = "mDependentPS", n = 32L, covariances = testcovariances))
  expect_identical(monteCarloSimulation(family = "mDependentPS", n = 32.5, covariances = testcovariances),
                   monteCarloSimulation(family = "mDependentPS", n = 32L, covariances = testcovariances))
  expect_equal(monteCarloSimulation(family = "mDependentPS", n = 32L, covariances = testcovariances, r = 10L),
               monteCarloSimulation(n = 32L, r = 10L, seed = 32L,
                                    rand.gen = testrand.genMDependent,
                                    family = "mDependentPS", intervalSystem = "dyaLen", output = "vector",
                                    covariances = testcovariances), tolerance = 1e-14, check.attributes = FALSE)

  testbyfilter(testn = 32L, testr = 10L, testseed = 1L, testfamily = "mDependentPS",
               testintervalSystem = "dyaLen", testlengths = as.integer(2^(0:5)), testpenalty = "none",
               tolerance = 1e-14, testcovariances = testcovariances)

  testbyfilter(testn = 17L, testr = 8L, testseed = -1L, testfamily = "mDependentPS",
               testintervalSystem = "dyaLen", testlengths = as.integer(2^(0:4)), testpenalty = "none",
               tolerance = 1e-14, testcovariances = testcovariances)

  expect_error(monteCarloSimulation(family = "mDependentPS", n = 32L, rand.gen = "s",
                                    covariances = testcovariances))
  expect_error(monteCarloSimulation(family = "mDependentPS", n = 32L, rand.gen = function() {},
                                    covariances = testcovariances))
  expect_error(monteCarloSimulation(family = "mDependentPS", n = 32L,
                                    rand.gen = function(test) {rnorm(32)}, covariances = testcovariances))
  expect_error(monteCarloSimulation(family = "mDependentPS", n = 32L,
                                    rand.gen = function(data, test) {rnorm(32)},
                                    covariances = testcovariances))

  expect_error(monteCarloSimulation(family = "mDependentPS", n = 32L,
                                    rand.gen = function(data) {rnorm(31)}, covariances = testcovariances))
  expect_error(monteCarloSimulation(family = "mDependentPS", n = 32L,
                                    rand.gen = function(data) {c(rnorm(31), NA)},
                                    covariances = testcovariances))
  expect_error(monteCarloSimulation(family = "mDependentPS", n = 32L,
                                    rand.gen = function(data) {c(rnorm(31), "NA")},
                                    covariances = testcovariances))
  expect_error(monteCarloSimulation(family = "mDependentPS", n = 32L,
                                    rand.gen = function(data) {c(rnorm(31), Inf)},
                                    covariances = testcovariances))

  expect_identical(monteCarloSimulation(family = "mDependentPS", n = 34L, r = 10L, covariances = testcovariances),
                   monteCarloSimulation(family = "mDependentPS", n = 34L, r = 10L, covariances = testcovariances,
                                        intervalSystem = "dyaLen"))

  testbyfilter(testn = 28L, testr = 7L, testseed = 1L, testfamily = "mDependentPS",
               testintervalSystem = "all", testlengths = 1:28, testpenalty = "none",
               tolerance = 1e-14, testcovariances = testcovariances)
  testbyfilter(testn = 28L, testr = 6L, testseed = 1L, testfamily = "mDependentPS",
               testintervalSystem = "dyaPar", testlengths = as.integer(2^(0:4)), testpenalty = "none",
               tolerance = 1e-14, testcovariances = testcovariances)

  expect_error(monteCarloSimulation(family = "mDependentPS", output = "maximum", n = 32L, lengths = c(33L),
                                    covariances = testcovariances))

  expect_equal(monteCarloSimulation(family = "mDependentPS", n = 32L, r = 10L, intervalSystem = "all",
                                    output = "maximum", lengths = 2^(0:5), covariances = testcovariances),
               monteCarloSimulation(family = "mDependentPS", n = 32L, r = 10L, intervalSystem = "dyaLen",
                                    output = "maximum", covariances = testcovariances),
               check.attributes = FALSE)
  expect_equal(monteCarloSimulation(family = "mDependentPS", n = 32L, r = 10L, intervalSystem = "all",
                                    output = "maximum", lengths = 2^(3:4), covariances = testcovariances),
               monteCarloSimulation(family = "mDependentPS", n = 32L, r = 10L, intervalSystem = "dyaLen",
                                    output = "maximum", lengths = 2^(3:4), covariances = testcovariances),
               check.attributes = FALSE)

  testlengths <- c(2:5, 10:23)
  testbyfilter(testn = 28L, testr = 10L, testseed = 10L, testfamily = "mDependentPS",
               testintervalSystem = "all", testlengths = testlengths, testpenalty = "none",
               tolerance = 1e-14, testcovariances = testcovariances, allLengths = 1:28)

  testlengths <- c(2L, 4L, 8L)
  testbyfilter(testn = 28L, testr = 10L, testseed = 10L, testfamily = "mDependentPS",
               testintervalSystem = "dyaLen", testlengths = testlengths, testpenalty = "none",
               tolerance = 1e-14, testcovariances = testcovariances, allLengths = as.integer(2^(0:4)))

  testlengths <- c(2L, 4L, 16L)
  testbyfilter(testn = 28L, testr = 10L, testseed = 10L, testfamily = "mDependentPS",
               testintervalSystem = "dyaPar", testlengths = testlengths, testpenalty = "none",
               tolerance = 1e-14, testcovariances = testcovariances, allLengths = as.integer(2^(0:4)))

  expect_identical(monteCarloSimulation(family = "mDependentPS", n = 32L, r = 10L, covariances = testcovariances),
                   monteCarloSimulation(family = "mDependentPS", n = 32L, r = 10L, penalty = "none",
                                        covariances = testcovariances))
  expect_identical(monteCarloSimulation(family = "mDependentPS", output = "maximum", n = 32L, r = 10L,
                                        covariances = testcovariances),
                   monteCarloSimulation(family = "mDependentPS", output = "maximum", n = 32L, r = 10L,
                                        penalty = "sqrt", covariances = testcovariances))

  testlengths <- c(2:5, 10:23)
  testbyfilter(testn = 28L, testr = 10L, testseed = 10L, testfamily = "mDependentPS",
               testintervalSystem = "all", testlengths = testlengths, testpenalty = "sqrt",
               tolerance = 1e-14, testcovariances = testcovariances, allLengths = 1:28)

  testlengths <- c(2L, 4L, 8L)
  testbyfilter(testn = 28L, testr = 10L, testseed = 10L, testfamily = "mDependentPS",
               testintervalSystem = "dyaLen", testlengths = testlengths, testpenalty = "sqrt",
               tolerance = 1e-14, testcovariances = testcovariances, allLengths = as.integer(2^(0:4)))

  testlengths <- c(2L, 4L, 16L)
  testbyfilter(testn = 28L, testr = 10L, testseed = 10L, testfamily = "mDependentPS",
               testintervalSystem = "dyaPar", testlengths = testlengths, testpenalty = "sqrt",
               tolerance = 1e-14, testcovariances = testcovariances, allLengths = as.integer(2^(0:4)))

  testlengths <- c(2:5, 10:23)
  testbyfilter(testn = 28L, testr = 10L, testseed = 10L, testfamily = "mDependentPS",
               testintervalSystem = "all", testlengths = testlengths, testpenalty = "log",
               tolerance = 1e-14, testcovariances = testcovariances, allLengths = 1:28)

  testlengths <- c(2L, 4L, 8L)
  testbyfilter(testn = 28L, testr = 10L, testseed = 10L, testfamily = "mDependentPS",
               testintervalSystem = "dyaLen", testlengths = testlengths, testpenalty = "log",
               tolerance = 1e-14, testcovariances = testcovariances, allLengths = as.integer(2^(0:4)))

  testlengths <- c(2L, 4L, 16L)
  testbyfilter(testn = 28L, testr = 10L, testseed = 10L, testfamily = "mDependentPS",
               testintervalSystem = "dyaPar", testlengths = testlengths, testpenalty = "log",
               tolerance = 1e-14, testcovariances = testcovariances, allLengths = as.integer(2^(0:4)))

  expect_error(monteCarloSimulation(family = "mDependentPS", n = 32L, intervalsystem = "dyaPar",
                                    covariances = testcovariances))

  expect_identical(monteCarloSimulation(10L, r = 100, family = "mDependentPS", covariances = testcovariances * 2),
                   monteCarloSimulation(10L, r = 100, family = "mDependentPS", covariances = testcovariances))
  expect_error(monteCarloSimulation(10L, family = "mDependentPS"))
  expect_error(monteCarloSimulation(10L, family = "mDependentPS", covariances = c(testcovariances, "s")))
  expect_error(monteCarloSimulation(10L, family = "mDependentPS", covariances = c(0.01, testcovariances)))
  expect_error(monteCarloSimulation(10L, family = "mDependentPS", covariances = c(-10, testcovariances)))
  expect_error(monteCarloSimulation(10L, family = "mDependentPS", covariances = c(testcovariances, 0)))

  expect_identical(monteCarloSimulation(10L, r = 100, family = "mDependentPS", covariances = testcovariances * 2),
                   monteCarloSimulation(10L, r = 100, family = "mDependentPS", correlations = testcovariances))
  expect_equal(monteCarloSimulation(10L, r = 100, family = "mDependentPS", covariances = testcovariances * 2),
               monteCarloSimulation(10L, r = 100, family = "mDependentPS",
                                    correlations = testcovariances, sd = 3))
})

test_that("family 'jsmurf' works", {
  testthat::skip_on_cran()
  testfilter <- lowpassFilter::lowpassFilter(type = "bessel", param = list(pole = 4L, cutoff = 0.1), sr = 1e4)
  expect_error(monteCarloSimulation(family = "jsmurf", filter = testfilter))
  expect_error(monteCarloSimulation(family = "jsmurf", n = "s", filter = testfilter))
  expect_identical(monteCarloSimulation(family = "jsmurf", n = 32, r = 100L, filter = testfilter),
                   monteCarloSimulation(family = "jsmurf", n = 32L, r = 100L, filter = testfilter))
  expect_identical(monteCarloSimulation(family = "jsmurf", n = 32.5, r = 100L, filter = testfilter),
                   monteCarloSimulation(family = "jsmurf", n = 32L, r = 100L, filter = testfilter))
  expect_equal(monteCarloSimulation(family = "jsmurf", n = 32L, filter = testfilter, r = 10L),
               monteCarloSimulation(n = 32L, r = 10L, seed = 32L,
                                    rand.gen = testrand.genJsmurf,
                                    family = "jsmurf", intervalSystem = "dyaLen", output = "vector",
                                    filter = testfilter), tolerance = 1e-14, check.attributes = FALSE)

  testJsmurf(testn = 64L, testr = 10L, testseed = 1L, testfamily = "jsmurf",
             testintervalSystem = "dyaLen", testlengths = as.integer(2^(4:6)), testpenalty = "sqrt",
             tolerance = 1e-14, testfilter = testfilter)
  testJsmurfPS(testn = 64L, testr = 10L, testseed = 1L, testfamily = "jsmurfPS",
               testintervalSystem = "dyaLen", testlengths = as.integer(2^(4:6)), testpenalty = "sqrt",
               tolerance = 1e-14, testfilter = testfilter)
  testJsmurfLR(testn = 64L, testr = 10L, testseed = 1L, testfamily = "jsmurfLR",
               testintervalSystem = "dyaLen", testlengths = as.integer(2^(4:6)), testpenalty = "sqrt",
               tolerance = 1e-14, testfilter = testfilter)

  expect_identical(monteCarloSimulation(family = "jsmurfPS", n = 34L, r = 10L, filter = testfilter),
                   monteCarloSimulation(family = "jsmurfPS", n = 34L, r = 10L, filter = testfilter,
                                        intervalSystem = "dyaLen"))

  testJsmurf(testn = 32L, testr = 10L, testseed = 1L, testfamily = "jsmurf",
             testintervalSystem = "all", testlengths = as.integer(c(12:32)), testpenalty = "log",
             tolerance = 1e-14, testfilter = testfilter, allLengths = as.integer(c(12:32)))
  testJsmurfLR(testn = 64L, testr = 10L, testseed = 1L, testfamily = "jsmurfLR",
               testintervalSystem = "dyaPar", testlengths = as.integer(2^(4:6)), testpenalty = "none",
               tolerance = 1e-14, testfilter = testfilter)

  testlengths <- c(12:15, 23L)
  testJsmurfPS(testn = 28L, testr = 10L, testseed = 10L, testfamily = "jsmurfPS",
               testintervalSystem = "all", testlengths = testlengths, testpenalty = "none",
               tolerance = 1e-14, testfilter = testfilter, allLengths = 12:28)
  expect_error(monteCarloSimulation(n = 32L, r = 10L, filter = testfilter, family = "jsmurf",
                                    output = "maximum", lengths = c(8, 16)))

  ret <- monteCarloSimulation(10L, family = "jsmurfPS", filter = testfilter, r = 10L)
  compare <- matrix(ncol = 10, nrow = 0)
  class(compare) <- "MCSimulationVector"
  expect_equal(ret, compare, check.attributes = FALSE)

  testfilter2 <- lowpassFilter::lowpassFilter(param = list(pole = 4, cutoff = 0.1), sr = 10, shift = 0.8)
  expect_identical(monteCarloSimulation(10L, r = 100, family = "jsmurfPS", filter = testfilter),
                   monteCarloSimulation(10L, r = 100, family = "jsmurfPS", filter = testfilter2))

  expect_identical(monteCarloSimulation(10L, r = 100, family = "jsmurf", filter = testfilter),
                   monteCarloSimulation(10L, r = 100, family = "jsmurf", filter = testfilter, sd = 3))

  testfilter <- lowpassFilter::lowpassFilter(param = list(pole = 3, cutoff = 0.025), len = 15, shift = 0.12)
  testJsmurf(testn = 64L, testr = 10L, testseed = 1L, testfamily = "jsmurf",
             testintervalSystem = "dyaLen", testlengths = as.integer(2^(4:6)), testpenalty = "none",
             tolerance = 1e-14, testfilter = testfilter)
})

test_that("family 'hjsmurf' works", {
  testthat::skip_on_cran()
  testfilter <- lowpassFilter::lowpassFilter(type = "bessel", param = list(pole = 4L, cutoff = 0.1), sr = 1e4)
  expect_error(monteCarloSimulation(family = "hjsmurf", filter = testfilter))
  expect_error(monteCarloSimulation(family = "hjsmurf", n = "s", filter = testfilter))
  expect_identical(monteCarloSimulation(family = "hjsmurf", n = 32, r = 100L, filter = testfilter),
                   monteCarloSimulation(family = "hjsmurf", n = 32L, r = 100L, filter = testfilter))
  expect_identical(monteCarloSimulation(family = "hjsmurf", n = 32.5, r = 100L, filter = testfilter),
                   monteCarloSimulation(family = "hjsmurf", n = 32L, r = 100L, filter = testfilter))
  expect_equal(monteCarloSimulation(family = "hjsmurf", n = 32L, filter = testfilter, r = 10L),
               monteCarloSimulation(n = 32L, r = 10L, seed = 32L,
                                    rand.gen = testrand.genJsmurf,
                                    family = "hjsmurf", intervalSystem = "dyaLen", output = "vector",
                                    filter = testfilter), tolerance = 1e-14, check.attributes = FALSE)

  testHjsmurf(testn = 64L, testr = 10L, testseed = 1L, testfamily = "hjsmurf",
              testintervalSystem = "dyaLen", testlengths = as.integer(2^(4:6)), testpenalty = "none",
              tolerance = 1e-14, testfilter = testfilter)
  testHjsmurfSPS(testn = 64L, testr = 10L, testseed = 1L, testfamily = "hjsmurfSPS",
                 testintervalSystem = "dyaLen", testlengths = as.integer(2^(4:6)), testpenalty = "none",
                 tolerance = 1e-14, testfilter = testfilter)
  testHjsmurfLR(testn = 64L, testr = 10L, testseed = 1L, testfamily = "hjsmurfLR",
                testintervalSystem = "dyaLen", testlengths = as.integer(2^(4:6)), testpenalty = "none",
                tolerance = 1e-14, testfilter = testfilter)

  expect_identical(monteCarloSimulation(family = "hjsmurfSPS", n = 34L, r = 10L, filter = testfilter),
                   monteCarloSimulation(family = "hjsmurfSPS", n = 34L, r = 10L, filter = testfilter,
                                        intervalSystem = "dyaLen"))

  testHjsmurf(testn = 32L, testr = 10L, testseed = 1L, testfamily = "hjsmurf",
              testintervalSystem = "all", testlengths = as.integer(c(13:32)), testpenalty = "log",
              tolerance = 1e-14, testfilter = testfilter, allLengths = as.integer(c(13:32)))
  testHjsmurfLR(testn = 64L, testr = 10L, testseed = 1L, testfamily = "hjsmurfLR",
                testintervalSystem = "dyaPar", testlengths = as.integer(2^(4:6)), testpenalty = "none",
                tolerance = 1e-14, testfilter = testfilter)

  testlengths <- c(13:15, 23L)
  testHjsmurfSPS(testn = 28L, testr = 10L, testseed = 10L, testfamily = "hjsmurfSPS",
                 testintervalSystem = "all", testlengths = testlengths, testpenalty = "none",
                 tolerance = 1e-14, testfilter = testfilter, allLengths = 13:28)
  expect_error(monteCarloSimulation(n = 32L, r = 10L, filter = testfilter, family = "hjsmurf",
                                    output = "maximum", lengths = c(8, 16)))

  ret <- monteCarloSimulation(10L, family = "hjsmurfSPS", filter = testfilter, r = 10L)
  compare <- matrix(ncol = 10, nrow = 0)
  class(compare) <- "MCSimulationVector"
  expect_equal(ret, compare, check.attributes = FALSE)

  testfilter2 <- lowpassFilter::lowpassFilter(param = list(pole = 4, cutoff = 0.1), sr = 10, shift = 0.8)
  expect_identical(monteCarloSimulation(10L, r = 100, family = "hjsmurfSPS", filter = testfilter),
                   monteCarloSimulation(10L, r = 100, family = "hjsmurfSPS", filter = testfilter2))

  testfilter <- lowpassFilter::lowpassFilter(param = list(pole = 3, cutoff = 0.025), len = 14, shift = 0.12)
  testHjsmurf(testn = 64L, testr = 10L, testseed = 1L, testfamily = "hjsmurf",
             testintervalSystem = "dyaLen", testlengths = as.integer(2^(4:6)), testpenalty = "none",
             tolerance = 1e-14, testfilter = testfilter)
})

test_that("family 'LR' works", {
  testthat::skip_on_cran()
  testfilter <- lowpassFilter::lowpassFilter(type = "bessel", param = list(pole = 4L, cutoff = 0.1), sr = 1e4)
  expect_error(monteCarloSimulation(family = "LR", filter = testfilter))
  expect_error(monteCarloSimulation(family = "LR", n = "s", filter = testfilter))
  expect_identical(monteCarloSimulation(family = "LR", n = 32, r = 10L, filter = testfilter),
                   monteCarloSimulation(family = "LR", n = 32L, r = 10L, filter = testfilter))
  expect_identical(monteCarloSimulation(family = "LR", n = 32.5, r = 10L, filter = testfilter),
                   monteCarloSimulation(family = "LR", n = 32L, r = 10L, filter = testfilter))
  # takes very long
  # expect_equal(monteCarloSimulation(family = "LR", n = 32L, filter = testfilter),
  #              monteCarloSimulation(n = 32L, r = 1e4L, seed = 32L, rand.gen = testrand.genJsmurf,
  #                                   family = "LR", intervalSystem = "all", output = "vector",
  #                                   filter = testfilter), tolerance = 1e-14, check.attributes = FALSE)

  testLR(testn = 32L, testr = 10L, testseed = 1L, testfamily = "LR", tolerance = 1e-12,
         testfilter = testfilter, allLengths = 1:20)

  expect_error(monteCarloSimulation(family = "LR", n = 32L, filter = testfilter,
                                    rand.gen = function(data, test) {rnorm(32)}))
  expect_error(monteCarloSimulation(family = "LR", n = 32L,
                                    rand.gen = function(data) {rnorm(31)}, filter = testfilter))

  expect_identical(monteCarloSimulation(family = "LR", n = 34L, r = 10L, filter = testfilter),
                   monteCarloSimulation(family = "LR", n = 34L, r = 10L, filter = testfilter,
                                        intervalSystem = "all"))

  expect_error(monteCarloSimulation(family = "LR", output = "maximum", n = 32L, lengths = c(33L),
                                    filter = testfilter))

  expect_identical(monteCarloSimulation(family = "LR", n = 32L, r = 10L, filter = testfilter),
                   monteCarloSimulation(family = "LR", n = 32L, r = 10L, penalty = "none",
                                        filter = testfilter))
  expect_identical(monteCarloSimulation(family = "LR", output = "maximum", n = 32L, r = 10L,
                                        filter = testfilter),
                   monteCarloSimulation(family = "LR", output = "maximum", n = 32L, r = 10L,
                                        penalty = "weights", filter = testfilter))

  expect_error(monteCarloSimulation(family = "LR", n = 32L, intervalsystem = "dyaPar",
                                    filter = testfilter))
  expect_error(monteCarloSimulation(family = "LR", n = 32L, intervalSystem = "dyaPar",
                                    filter = testfilter))

  ret <- monteCarloSimulation(20L, family = "LR", filter = testfilter, r = 10L)
  compare <- matrix(-Inf, ncol = 10, nrow = 20)
  class(compare) <- "MCSimulationVector"
  expect_equal(ret, compare, check.attributes = FALSE)


  expect_error(monteCarloSimulation(family = "LR", n = 32L, filter = unclass(testfilter)))

  expect_identical(monteCarloSimulation(family = "LR", n = 32L, r = 10L, filter = testfilter, fit = stepblock(1)),
                   monteCarloSimulation(family = "LR", n = 32L, r = 10L, filter = testfilter))

  expect_identical(monteCarloSimulation(family = "LR", n = 32L, r = 10L, filter = testfilter,
                                        startTime = 0, thresholdLongSegment = 10L, localValue = stats::median,
                                        regularization = 1, suppressWarningNoDeconvolution = FALSE, localList = NULL),
                   monteCarloSimulation(family = "LR", n = 32L, r = 10L, filter = testfilter))


  expect_error(monteCarloSimulation(family = "LR", n = 32L, r = 10L, filter = testfilter, startTime = NA))
  expect_error(monteCarloSimulation(family = "LR", n = 32L, r = 10L, filter = testfilter, startTime = c(0, 1)))

  expect_identical(monteCarloSimulation(family = "LR", n = 32L, r = 10L, filter = testfilter, startTime = -100),
                   monteCarloSimulation(family = "LR", n = 32L, r = 10L, filter = testfilter, startTime = 0))
  expect_identical(monteCarloSimulation(family = "LR", n = 32L, r = 10L, filter = testfilter, startTime = 0.001),
                   monteCarloSimulation(family = "LR", n = 32L, r = 10L, filter = testfilter, startTime = 0))

  expect_error(monteCarloSimulation(family = "LR", n = 32L, r = 10L, filter = testfilter, thresholdLongSegment = NA))
  expect_error(monteCarloSimulation(family = "LR", n = 32L, r = 10L, filter = testfilter, thresholdLongSegment = -1L))
  expect_error(monteCarloSimulation(family = "LR", n = 32L, r = 10L, filter = testfilter,
                           thresholdLongSegment = c(10L, 23L)))

  expect_identical(monteCarloSimulation(family = "LR", n = 32L, r = 10L, filter = testfilter,
                                        thresholdLongSegment = 10.5),
                   monteCarloSimulation(family = "LR", n = 32L, r = 10L, filter = testfilter))


  testLR(testn = 32L, testr = 10L, testseed = 1L, testfamily = "LR", tolerance = 1e-12,
         testfilter = testfilter, allLengths = 1:20, thresholdLongSegment = 8L, save = FALSE)
  testLR(testn = 32L, testr = 10L, testseed = 1L, testfamily = "LR", tolerance = 1e-12,
         testfilter = testfilter, allLengths = 1:20, thresholdLongSegment = 12L, save = FALSE)
  testLR(testn = 30L, testr = 10L, testseed = 1L, testfamily = "LR", tolerance = 1e-12,
         testfilter = testfilter, allLengths = 1:20)

  expect_error(monteCarloSimulation(family = "LR", n = 32L, r = 10L, filter = testfilter, localValue = 1))
  expect_error(monteCarloSimulation(family = "LR", n = 32L, r = 10L, filter = testfilter, localValue = function() 1))
  expect_error(monteCarloSimulation(family = "LR", n = 32L, r = 10L, filter = testfilter, localValue = function(a) NA))
  expect_error(monteCarloSimulation(family = "LR", n = 32L, r = 10L, filter = testfilter,
                                    localValue = function(a) c(1, 2)))
  testLR(testn = 32L, testr = 10L, testseed = 1L, testfamily = "LR", tolerance = 1e-12,
         testfilter = testfilter, allLengths = 1:20, localValue = mean, save = FALSE)

  expect_error(monteCarloSimulation(family = "LR", n = 32L, r = 10L, filter = testfilter, regularization = NULL))
  expect_error(monteCarloSimulation(family = "LR", n = 32L, r = 10L, filter = testfilter, regularization = c(1, NA)))
  testLR(testn = 32L, testr = 10L, testseed = 1L, testfamily = "LR", tolerance = 1e-12,
         testfilter = testfilter, allLengths = 1:20, regularization = 2, save = TRUE)
  testLR(testn = 32L, testr = 10L, testseed = 1L, testfamily = "LR", tolerance = 1e-12,
         testfilter = testfilter, allLengths = 1:20, regularization = c(2, 1), save = TRUE)

  expect_error(monteCarloSimulation(family = "LR", n = 32L, r = 10L, filter = testfilter, correlations = "s"))
  expect_error(monteCarloSimulation(family = "LR", n = 32L, r = 10L, filter = testfilter, correlations = c(Inf, 3)))
  expect_error(monteCarloSimulation(family = "LR", n = 32L, r = 10L, filter = testfilter, correlations = c(1, 0)))
  testcor <- testfilter$acf
  testcor[1] <- testcor[1] * 3
  expect_identical(monteCarloSimulation(family = "LR", n = 32L, r = 10L, filter = testfilter,
                                        correlations = testcor),
                   monteCarloSimulation(family = "LR", n = 32L, r = 10L, filter = testfilter,
                                        regularization = 2))

  testcor <- testfilter$acf
  testcor[1] <- testcor[1] * 2
  ret <- monteCarloSimulation(family = "LR", n = 32L, r = 10L, filter = testfilter, correlations = testcor)
  expect_true(attr(ret, "save"))

  expect_error(monteCarloSimulation(family = "LR", n = 32L, r = 10L, filter = testfilter,
                           suppressWarningNoDeconvolution = 1))
  expect_error(monteCarloSimulation(family = "LR", n = 32L, r = 10L, filter = testfilter,
                           suppressWarningNoDeconvolution = as.logical(NA)))
  expect_error(monteCarloSimulation(family = "LR", n = 32L, r = 10L, filter = testfilter,
                           suppressWarningNoDeconvolution = c(TRUE, TRUE)))

  expect_identical(monteCarloSimulation(family = "LR", n = 32L, r = 10L, filter = testfilter,
                                        suppressWarningNoDeconvolution = TRUE),
                   monteCarloSimulation(family = "LR", n = 32L, r = 10L, filter = testfilter))
})

test_that("family '2Param' works", {
  testthat::skip_on_cran()
  testfilter <- lowpassFilter::lowpassFilter(type = "bessel", param = list(pole = 4L, cutoff = 0.1), sr = 1e4)
  expect_error(monteCarloSimulation(family = "2Param", filter = testfilter))
  expect_error(monteCarloSimulation(family = "2Param", n = "s", filter = testfilter))
  expect_identical(monteCarloSimulation(family = "2Param", n = 64, r = 3L, filter = testfilter),
                   monteCarloSimulation(family = "2Param", n = 64L, r = 3L, filter = testfilter))
  expect_identical(monteCarloSimulation(family = "2Param", n = 64.5, r = 3L, filter = testfilter),
                   monteCarloSimulation(family = "2Param", n = 64L, r = 3L, filter = testfilter))
  # only r = 2L, since it takes very long
  # expect_equal(monteCarloSimulation(family = "2Param", n = 64L, filter = testfilter, r = 2L),
  #              monteCarloSimulation(n = 64L, r = 2L, seed = 64L, rand.gen = testrand.genJsmurf,
  #                                   family = "2Param", intervalSystem = "all", output = "vector",
  #                                   filter = testfilter), tolerance = 1e-14, check.attributes = FALSE)


  test2Param(testn = 64L, testr = 2L, testseed = 1L, testfamily = "2Param", tolerance = 1e-12,
             testfilter = testfilter, allLengths = 1:64)

  expect_error(monteCarloSimulation(family = "2Param", n = 32L, filter = testfilter,
                                    rand.gen = function(data, test) {rnorm(32)}))
  expect_error(monteCarloSimulation(family = "2Param", n = 32L,
                                    rand.gen = function(data) {rnorm(31)}, filter = testfilter))

  expect_error(monteCarloSimulation(family = "2Param", output = "maximum", n = 32L, lengths = c(33L),
                                    filter = testfilter))

  expect_error(monteCarloSimulation(family = "2Param", n = 32L, intervalsystem = "dyaPar",
                                    filter = testfilter))
  expect_error(monteCarloSimulation(family = "2Param", n = 32L, intervalSystem = "dyaPar",
                                    filter = testfilter))

  ret <- monteCarloSimulation(20L, family = "2Param", filter = testfilter, r = 2L)
  compare <- matrix(-Inf, ncol = 2, nrow = 20)
  class(compare) <- "MCSimulationVector"
  expect_equal(ret, compare, check.attributes = FALSE)


  expect_error(monteCarloSimulation(family = "2Param", n = 32L, filter = unclass(testfilter)))

  expect_identical(monteCarloSimulation(family = "2Param", n = 64L, r = 1L, filter = testfilter, fit = stepblock(1)),
                   monteCarloSimulation(family = "2Param", n = 64L, r = 1L, filter = testfilter))

  filter <- testfilter
  expect_identical(monteCarloSimulation(family = "2Param", n = 64L, r = 2L, filter = testfilter,
                                        startTime = 0, thresholdLongSegment = 25L, localValue = stats::median,
                                        localVar = function(data) sdrobnorm(data, lag = filter$len + 1L)^2,
                                        regularization = 1, suppressWarningNoDeconvolution = FALSE, localList = NULL),
                   monteCarloSimulation(family = "2Param", n = 64L, r = 2L, filter = testfilter))


  expect_error(monteCarloSimulation(family = "2Param", n = 32L, r = 10L, filter = testfilter, startTime = NA))
  expect_error(monteCarloSimulation(family = "2Param", n = 32L, r = 10L, filter = testfilter, startTime = c(0, 1)))

  expect_identical(monteCarloSimulation(family = "2Param", n = 32L, r = 10L, filter = testfilter, startTime = -0.001),
                   monteCarloSimulation(family = "2Param", n = 32L, r = 10L, filter = testfilter, startTime = 0))
  expect_identical(monteCarloSimulation(family = "2Param", n = 32L, r = 10L, filter = testfilter, startTime = 100),
                   monteCarloSimulation(family = "2Param", n = 32L, r = 10L, filter = testfilter, startTime = 0))

  expect_error(monteCarloSimulation(family = "2Param", n = 32L, r = 10L, filter = testfilter,
                                    thresholdLongSegment = NA))
  expect_error(monteCarloSimulation(family = "2Param", n = 32L, r = 10L, filter = testfilter,
                                    thresholdLongSegment = -1L))
  expect_error(monteCarloSimulation(family = "2Param", n = 32L, r = 10L, filter = testfilter,
                                    thresholdLongSegment = c(10L, 23L)))

  expect_identical(monteCarloSimulation(family = "2Param", n = 64L, r = 1L, filter = testfilter,
                                        thresholdLongSegment = 25.5),
                   monteCarloSimulation(family = "2Param", n = 64L, r = 1L, filter = testfilter))

  test2Param(testn = 64L, testr = 2L, testseed = 1L, testfamily = "2Param", tolerance = 1e-12,
         testfilter = testfilter, allLengths = 1:64, thresholdLongSegment = 20L, save = FALSE)

  expect_error(monteCarloSimulation(family = "2Param", n = 32L, r = 10L, filter = testfilter, localValue = 1))
  expect_error(monteCarloSimulation(family = "2Param", n = 32L, r = 10L, filter = testfilter,
                                    localValue = function() 1))
  expect_error(monteCarloSimulation(family = "2Param", n = 64L, r = 10L, filter = testfilter,
                                    localValue = function(a) NA))
  expect_error(monteCarloSimulation(family = "2Param", n = 64L, r = 10L, filter = testfilter,
                                    localValue = function(a) c(1, 2)))

  expect_error(monteCarloSimulation(family = "2Param", n = 32L, r = 10L, filter = testfilter, localVar = 1))
  expect_error(monteCarloSimulation(family = "2Param", n = 32L, r = 10L, filter = testfilter,
                                    localVar = function() 1))
  expect_error(monteCarloSimulation(family = "2Param", n = 64L, r = 10L, filter = testfilter,
                                    localVar = function(a) NA))
  expect_error(monteCarloSimulation(family = "2Param", n = 64L, r = 10L, filter = testfilter,
                                    localVar = function(a) c(1, 2)))
  expect_error(monteCarloSimulation(family = "2Param", n = 64L, r = 10L, filter = testfilter,
                                    localVar = function(a) -1))

  test2Param(testn = 64L, testr = 1L, testseed = 1L, testfamily = "2Param", tolerance = 1e-12,
         testfilter = testfilter, allLengths = 1:64, localValue = mean, localVar = function(a) 1, save = FALSE)

  expect_error(monteCarloSimulation(family = "2Param", n = 32L, r = 10L, filter = testfilter, regularization = NULL))
  expect_error(monteCarloSimulation(family = "2Param", n = 32L, r = 10L, filter = testfilter,
                                    regularization = c(1, NA)))
  test2Param(testn = 64L, testr = 1L, testseed = 1L, testfamily = "2Param", tolerance = 1e-12,
         testfilter = testfilter, allLengths = 1:64, regularization = c(2, 1), save = TRUE)

  expect_error(monteCarloSimulation(family = "2Param", n = 32L, r = 10L, filter = testfilter, correlations = "s"))
  expect_error(monteCarloSimulation(family = "2Param", n = 32L, r = 10L, filter = testfilter, correlations = c(Inf, 3)))
  expect_error(monteCarloSimulation(family = "2Param", n = 32L, r = 10L, filter = testfilter, correlations = c(1, 0)))
  testcor <- testfilter$acf
  testcor[1] <- testcor[1] * 3
  expect_identical(monteCarloSimulation(family = "2Param", n = 64L, r = 1L, filter = testfilter,
                                        correlations = testcor),
                   monteCarloSimulation(family = "2Param", n = 64L, r = 1L, filter = testfilter,
                                        regularization = 2))

  testcor <- testfilter$acf
  testcor[1] <- testcor[1] * 2
  ret <- monteCarloSimulation(family = "2Param", n = 64L, r = 1L, filter = testfilter, correlations = testcor)
  expect_true(attr(ret, "save"))

  expect_error(monteCarloSimulation(family = "2Param", n = 64L, r = 10L, filter = testfilter,
                                    suppressWarningNoDeconvolution = 1))
  expect_error(monteCarloSimulation(family = "2Param", n = 32L, r = 10L, filter = testfilter,
                                    suppressWarningNoDeconvolution = as.logical(NA)))
  expect_error(monteCarloSimulation(family = "2Param", n = 32L, r = 10L, filter = testfilter,
                                    suppressWarningNoDeconvolution = c(TRUE, TRUE)))
})
