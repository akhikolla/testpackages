context("computeStat")

source(system.file("tests/comparisons/indices.R", package = "stepR"))
source(system.file("tests/comparisons/singleStat.R", package = "stepR"))

statAll <- function(y, left, right, value, indices, singleStat, lengths, ...) {
  indices <- indices(length(y), lengths = lengths, ...)

  stat <- rep(-Inf, length(y))
  for (i in 1:length(indices$li)) {
    index <- sum(left <= indices$li[i])
    if (index == length(right) - sum(right >= indices$ri[i]) + 1) {
      stat[indices$ri[i] - indices$li[i] + 1] <- max(stat[indices$ri[i] - indices$li[i] + 1],
                                                     singleStat(y, value[index],
                                                                indices$li[i], indices$ri[i], ...))
    }
  }
  stat
}

test <- function(testy = NULL, testsignal = 0, left = NULL, right = NULL, value = NULL,
                 testfamily = NULL, testintervalSystem = NULL, testlengths = NULL,
                 finallengths = NULL, testnq = length(testy), defaultPenalty = NULL, indices = NULL,
                 singleStat = NULL, tolerance = 1e-12, shift = 0, ...) {
  # none
  compare <- statAll(testy, left = left, right = right, value = value, indices = indices,
                     singleStat = singleStat, lengths = finallengths, ...)[finallengths]

  expect_equal(computeStat(testy, signal = testsignal, family = testfamily,
                           intervalSystem = testintervalSystem, lengths = testlengths, penalty = "none",
                           nq = testnq, output = "maximum", ...),
               max(compare), info = "none", tolerance = tolerance)
  expect_equal(computeStat(testy, signal = testsignal, family = testfamily,
                           intervalSystem = testintervalSystem, lengths = testlengths, penalty = "none",
                           nq = testnq, output = "vector", ...),
               compare, info = "none", tolerance = tolerance)
  ret <- computeStat(testy, signal = testsignal, family = testfamily, intervalSystem = testintervalSystem,
                     lengths = testlengths, penalty = "none", nq = testnq, output = "list", ...)
  expect_equal(ret, list(maximum = max(compare), stat = compare, lengths = finallengths),
               info = "none", tolerance = tolerance)
  expect_identical(computeStat(testy, signal = testsignal, family = testfamily,
                               intervalSystem = testintervalSystem, lengths = testlengths, penalty = "none",
                               nq = testnq, ...), ret, info = "none")

  # sqrt
  compare <- statAll(testy, left = left, right = right, value = value, indices = indices,
                     singleStat = singleStat, lengths = finallengths, ...)[finallengths]
  finite <- compare != -Inf
  compare[finite] <- sqrt(2 * compare[finite]) - sqrt(2 * log(exp(1) * testnq / (finallengths[finite] + shift)))

  expect_equal(computeStat(testy, signal = testsignal, family = testfamily,
                           intervalSystem = testintervalSystem, lengths = testlengths, penalty = "sqrt",
                           nq = testnq, output = "maximum", ...),
               max(compare), info = "sqrt", tolerance = tolerance)
  expect_equal(computeStat(testy, signal = testsignal, family = testfamily,
                           intervalSystem = testintervalSystem, lengths = testlengths, penalty = "sqrt",
                           nq = testnq, output = "vector", ...),
               compare, info = "sqrt", tolerance = tolerance)
  ret <- computeStat(testy, signal = testsignal, family = testfamily, intervalSystem = testintervalSystem,
                     lengths = testlengths, penalty = "sqrt", nq = testnq, output = "list", ...)
  expect_equal(ret, list(maximum = max(compare), stat = compare, lengths = finallengths),
               info = "sqrt", tolerance = tolerance)
  expect_identical(computeStat(testy, signal = testsignal, family = testfamily, penalty = "sqrt",
                               intervalSystem = testintervalSystem, lengths = testlengths, nq = testnq, ...),
                   ret, info = "sqrt")

  # log
  compare <- statAll(testy, left = left, right = right, value = value, indices = indices,
                     singleStat = singleStat, lengths = finallengths, ...)[finallengths] -
    log(exp(1) * testnq / (finallengths + shift))

  expect_equal(computeStat(testy, signal = testsignal, family = testfamily,
                           intervalSystem = testintervalSystem, lengths = testlengths, penalty = "log",
                           nq = testnq, output = "maximum", ...),
               max(compare), info = "log", tolerance = tolerance)
  expect_equal(computeStat(testy, signal = testsignal, family = testfamily,
                           intervalSystem = testintervalSystem, lengths = testlengths, penalty = "log",
                           nq = testnq, output = "vector", ...),
               compare, info = "log", tolerance = tolerance)
  ret <- computeStat(testy, signal = testsignal, family = testfamily, intervalSystem = testintervalSystem,
                     lengths = testlengths, penalty = "log", nq = testnq, output = "list", ...)
  expect_equal(ret, list(maximum = max(compare), stat = compare, lengths = finallengths),
               info = "log", tolerance = tolerance)
  expect_identical(computeStat(testy, signal = testsignal, family = testfamily,
                               intervalSystem = testintervalSystem, lengths = testlengths, penalty = "log",
                               nq = testnq, ...), ret, info = "log")

  # default penalty
  expect_identical(computeStat(testy, signal = testsignal, family = testfamily,
                               intervalSystem = testintervalSystem, lengths = testlengths,
                               output = "maximum", nq = testnq, ...),
                   computeStat(testy, signal = testsignal, family = testfamily,
                               intervalSystem = testintervalSystem, lengths = testlengths,
                               output = "maximum", nq = testnq, penalty = defaultPenalty, ...))
  expect_identical(computeStat(testy, signal = testsignal, family = testfamily,
                               intervalSystem = testintervalSystem, lengths = testlengths,
                               output = "vector", nq = testnq, ...),
                   computeStat(testy, signal = testsignal, family = testfamily,
                               intervalSystem = testintervalSystem, lengths = testlengths,
                               output = "vector", nq = testnq, penalty = defaultPenalty, ...))
  expect_identical(computeStat(testy, signal = testsignal, family = testfamily,
                               intervalSystem = testintervalSystem, lengths = testlengths,
                               output = "list", nq = testnq, ...),
                   computeStat(testy, signal = testsignal, family = testfamily,
                               intervalSystem = testintervalSystem, lengths = testlengths,
                               output = "list", nq = testnq, penalty = defaultPenalty, ...))
  expect_identical(computeStat(testy, signal = testsignal, family = testfamily,
                               intervalSystem = testintervalSystem, lengths = testlengths,
                               nq = testnq, ...),
                   computeStat(testy, signal = testsignal, family = testfamily,
                               intervalSystem = testintervalSystem, lengths = testlengths,
                               output = "list", nq = testnq, penalty = defaultPenalty, ...))
}

test_that("argument y is tested", {
  testn <- 100L
  testy <- rnorm(testn)

  expect_error(computeStat())
  expect_error(computeStat(numeric(0)))
  expect_identical(computeStat(testy), computeStat(testy, family = "gauss", intervalSystem = "all",
                                                   lengths = 1:testn, penalty = "sqrt", nq = testn,
                                                   output = "list", sd = sdrobnorm(testy)))
  expect_error(computeStat(as.integer(testy)))
  expect_error(computeStat(c(testy, "s")))
  expect_error(computeStat(c(rnorm(10), NA)))

  test(testy, left = 1L, right = testn, value = 0, finallengths = 1:testn,
       indices = indicesAll, singleStat = singleStatGauss, sd = 1, defaultPenalty = "sqrt")
})

test_that("argument signal is tested and works in the default case", {
  testn <- 34L
  testy <- c(rnorm(testn / 2), rnorm(testn / 2, 1))
  testsignal <- list(leftIndex = c(1L, 12L), rightIndex = c(11L, 34L), value = c(0, 1))

  expect_error(computeStat(testy, signal = NULL))
  expect_error(computeStat(testy, signal = "s"))
  expect_identical(computeStat(testy, signal = 0L), computeStat(testy, signal = 0))
  expect_error(computeStat(testy, signal = c(1, 2)))
  expect_error(computeStat(testy, signal = NA))
  expect_error(computeStat(testy, signal = Inf))
  expect_identical(computeStat(testy, signal = 1),
                   computeStat(testy, signal = list(leftIndex = 1L, rightIndex = 34L, value = 1)))
  expect_equal(computeStat(testy, signal = 0),
               computeStat(testy, signal = list(leftIndex = 1L, rightIndex = 34L, value = 0)),
               tolerance = 1e-14)

  expect_error(computeStat(testy, signal = list()))
  expect_error(computeStat(testy, signal = list(leftIndex = 1L,
                                                rightIndex = 34L, value = NULL)))
  expect_error(computeStat(testy, signal = list(leftIndex = 1L,
                                                rightIndex = c(17L, 34L), value = numeric(0))))
  expect_error(computeStat(testy, signal = list(leftIndex = c(1L, 18L),
                                                rightIndex = c(17L, 34L), value = numeric(0))))
  expect_error(computeStat(testy, signal = list(leftIndex = c(1L, Inf),
                                                rightIndex = c(17L, 34L), value = c(1, 2))))
  expect_error(computeStat(testy, signal = list(leftIndex = c(1L, 18L),
                                                rightIndex = c(17L, 34L), value = c(1, NA))))
  expect_identical(computeStat(testy, signal = list(leftIndex = 1L, rightIndex = 34L, value = 1L)),
                   computeStat(testy, signal = list(leftIndex = 1L, rightIndex = 34L, value = 1)))
  expect_identical(computeStat(testy, signal = list(leftIndex = 1, rightIndex = 34, value = 1)),
                   computeStat(testy, signal = list(leftIndex = 1L, rightIndex = 34L, value = 1)))

  expect_error(computeStat(testy, signal = list(leftIndex = 2L, rightIndex = 35L, value = 1)))
  expect_error(computeStat(testy, signal = list(leftIndex = c(1L, 18L),
                                                rightIndex = c(17L, 35L), value = c(1, 2))))

  test(testy, left = 1L, right = testn, value = 0, testfamily = "gauss", finallengths = 1:testn,
       indices = indicesAll, singleStat = singleStatGauss, sd = 1.23,
       defaultPenalty = "sqrt")
  test(testy, testsignal = 1, left = 1L, right = testn, value = 1, testfamily = "gauss",
       finallengths = 1:testn, indices = indicesAll, singleStat = singleStatGauss,
       sd = 1.23, defaultPenalty = "sqrt")
  test(testy, testsignal = testsignal, left = c(1L, 12L), right = c(11L, 34L), value = c(0, 1),
       testfamily = "gauss", finallengths = 1:testn,
       indices = indicesAll, singleStat = singleStatGauss, sd = 1.23, defaultPenalty = "sqrt")

  testsignal <- list(leftIndex = c(1L, 12L), rightIndex = c(11L, 34L), value = c(0, 1))
  testsignal2 <- list(lefxx = "s", rtg = 1:10, leftIndex = c(1L, 12L), rightIndex = c(11L, 34L), value = c(0, 1))
  expect_identical(computeStat(testy, signal = testsignal2), computeStat(testy, signal = testsignal))

  teststepfit <- stepfit(23, family = "gauss", value = c(0, 1), param = 2, leftEnd = c(1L, 12L) / 34,
                         rightEnd = c(11L, 34L) / 34, x0 = -1, leftIndex = c(1L, 12L), rightIndex = c(11L, 34L))
  expect_identical(computeStat(testy, signal = teststepfit), computeStat(testy, signal = testsignal))
})

test_that("argument family is tested and works in the default case", { # other families below
  testn <- 34L
  testy <- rnorm(testn)

  expect_error(computeStat(testy, family = ""))
  expect_error(computeStat(testy, family = "bisom"))
  expect_error(computeStat(testy, family = c("gauss", "hsmuce")))
  expect_identical(computeStat(testy, signal = 0), computeStat(testy, signal = 0, family = "gauss"))
})

test_that("argument intervalSystem is tested and works", {
  testn <- 34L
  testy <- c(rnorm(testn / 2), rnorm(testn / 2, 1))
  testsignal <- list(leftIndex = c(1L, 12L), rightIndex = c(11L, 34L), value = c(0, 1))

  expect_error(computeStat(testy, intervalSystem = ""))
  expect_error(computeStat(testy, intervalSystem = "dya"))
  expect_error(computeStat(testy, intervalSystem = "dyalen"))
  expect_error(computeStat(testy, intervalSystem = "dyapar"))
  expect_error(computeStat(testy, intervalSystem = c("dyaPar", "all")))
  expect_identical(computeStat(testy, signal = 0), computeStat(testy, signal = 0, intervalSystem = "all"))
  expect_identical(computeStat(rep(testy, 30), signal = 0),
                   computeStat(rep(testy, 30), signal = 0, intervalSystem = "dyaLen"))

  test(testy, left = 1L, right = testn, value = 0, testfamily = "gauss", finallengths = 2^(0:5),
       testintervalSystem = "dyaLen", indices = indicesDyaLen, singleStat = singleStatGauss,
       sd = 1.23, defaultPenalty = "sqrt")
  test(testy, testsignal = 1, left = 1L, right = testn, value = 1, testfamily = "gauss",
       finallengths = 2^(0:5), testintervalSystem = "dyaLen",
       indices = indicesDyaLen, singleStat = singleStatGauss, sd = 1.23, defaultPenalty = "sqrt")
  test(testy, testsignal = testsignal, left = c(1L, 12L), right = c(11L, 34L), value = c(0, 1),
       testfamily = "gauss", finallengths = 2^(0:5), testintervalSystem = "dyaLen",
       indices = indicesDyaLen, singleStat = singleStatGauss, sd = 1.23, defaultPenalty = "sqrt")

  test(testy, left = 1L, right = testn, value = 0, testfamily = "gauss", finallengths = 2^(0:5),
       testintervalSystem = "dyaPar", indices = indicesDyaPar, singleStat = singleStatGauss,
       sd = 1.23, defaultPenalty = "sqrt")
  test(testy, testsignal = 1, left = 1L, right = testn, value = 1, testfamily = "gauss",
       finallengths = 2^(0:5), testintervalSystem = "dyaPar",
       indices = indicesDyaPar, singleStat = singleStatGauss, sd = 1.23, defaultPenalty = "sqrt")
  test(testy, testsignal = testsignal, left = c(1L, 12L), right = c(11L, 34L), value = c(0, 1),
       testfamily = "gauss", finallengths = 2^(0:5), testintervalSystem = "dyaPar",
       indices = indicesDyaPar, singleStat = singleStatGauss, sd = 1.23, defaultPenalty = "sqrt")
})

test_that("argument lengths is tested and works", {
  testn <- 45L
  testy <- c(rnorm(20, 1, 0.23), rnorm(testn - 20, -1, 0.34))
  testsignal <- list(leftIndex = c(1L, 23L), rightIndex = c(22L, 45L), value = c(0, -1))

  expect_error(computeStat(testy, lengths = "s"))
  expect_error(computeStat(testy, lengths = c(1:10, NA)))
  expect_error(computeStat(testy, lengths = c(1:10, Inf)))
  expect_error(computeStat(testy, lengths = 0:10))
  expect_error(computeStat(testy, lengths = -1L))
  expect_error(computeStat(testy, lengths = c(1L, 46L)))
  expect_error(computeStat(testy, intervalSystem = "dyaLen", lengths = 3L))
  expect_error(computeStat(testy, intervalSystem = "dyaLen", lengths = 64L))
  expect_error(computeStat(testy, intervalSystem = "dyaPar", lengths = 3L))
  expect_error(computeStat(testy, intervalSystem = "dyaPar", lengths = 64L))

  expect_warning(ret <- computeStat(testy, lengths = c(1:10, 10)))
  expect_identical(ret, computeStat(testy, lengths = c(1:10)))
  expect_identical(computeStat(testy, lengths = c(10:1)),
                   computeStat(testy, lengths = c(1:10)))
  expect_identical(computeStat(testy, lengths = c(1:10 + 0.5)),
                   computeStat(testy, lengths = c(1:10)))

  expect_equal(computeStat(testy, lengths = 2^(0:5)), computeStat(testy, intervalSystem = "dyaLen"))
  expect_equal(computeStat(testy, lengths = 2^(3:4)),
               computeStat(testy, intervalSystem = "dyaLen", lengths = 2^(3:4)))

  testlengths <- c(1, 3:5, 10:23)
  test(testy, left = 1L, right = testn, value = 0, testfamily = "gauss", testlengths = testlengths,
       finallengths = testlengths, indices = indicesAllLengths,
       singleStat = singleStatGauss, sd = 1.23, defaultPenalty = "sqrt")
  test(testy, testsignal = 1, left = 1L, right = testn, value = 1, testfamily = "gauss",
       testlengths = testlengths, finallengths = testlengths,
       indices = indicesAllLengths, singleStat = singleStatGauss, sd = 1.23, defaultPenalty = "sqrt")
  test(testy, testsignal = testsignal, left = c(1L, 23L), right = c(22L, 45L), value = c(0, -1),
       testfamily = "gauss", testlengths = testlengths, finallengths = testlengths,
       indices = indicesAllLengths, singleStat = singleStatGauss, sd = 1.23, defaultPenalty = "sqrt")

  testlengths <- c(1, 4, 32)
  test(testy, left = 1L, right = testn, value = 0, testfamily = "gauss", testlengths = testlengths,
       finallengths = testlengths, indices = indicesDyaLenLengths,
       singleStat = singleStatGauss, sd = 1.23, testintervalSystem = "dyaLen", defaultPenalty = "sqrt")
  test(testy, testsignal = 1, left = 1L, right = testn, value = 1, testfamily = "gauss",
       testlengths = testlengths, finallengths = testlengths,
       indices = indicesDyaLenLengths, singleStat = singleStatGauss, sd = 1.23,
       testintervalSystem = "dyaLen", defaultPenalty = "sqrt")
  test(testy, testsignal = testsignal, left = c(1L, 23L), right = c(22L, 45L), value = c(0, -1),
       testfamily = "gauss", testlengths = testlengths, finallengths = testlengths,
       indices = indicesDyaLenLengths, singleStat = singleStatGauss, sd = 1.23,
       testintervalSystem = "dyaLen", defaultPenalty = "sqrt")

  testlengths <- c(2, 4, 16)
  test(testy, left = 1L, right = testn, value = 0, testfamily = "gauss", testlengths = testlengths,
       finallengths = testlengths, indices = indicesDyaParLengths,
       singleStat = singleStatGauss, sd = 1.23, testintervalSystem = "dyaPar", defaultPenalty = "sqrt")
  test(testy, testsignal = 1, left = 1L, right = testn, value = 1, testfamily = "gauss",
       testlengths = testlengths, finallengths = testlengths,
       indices = indicesDyaParLengths, singleStat = singleStatGauss, sd = 1.23,
       testintervalSystem = "dyaPar", defaultPenalty = "sqrt")
  test(testy, testsignal = testsignal, left = c(1L, 23L), right = c(22L, 45L), value = c(0, -1),
       testfamily = "gauss", testlengths = testlengths, finallengths = testlengths,
       indices = indicesDyaParLengths, singleStat = singleStatGauss, sd = 1.23,
       testintervalSystem = "dyaPar", defaultPenalty = "sqrt")

  suppressWarnings(ret <- computeStat(testy, lengths = integer(0)))
  expect_identical(ret$maximum, -Inf)
  expect_identical(ret$stat, numeric(0))
  expect_identical(ret$lengths, integer(0))
})

test_that("argument penalty is tested and works in the default case", {
  testn <- 34L
  testy <- rnorm(testn)

  expect_error(computeStat(testy, penalty = ""))
  expect_error(computeStat(testy, penalty = c("gauss")))
  expect_error(computeStat(testy, penalty = c("sqrt", "log")))
  expect_identical(computeStat(testy, signal = 0), computeStat(testy, signal = 0, penalty = "sqrt"))
})

test_that("argument nq is tested and works", {
  testn <- 34L
  testy <- c(rnorm(testn / 2), rnorm(testn / 2, 1))
  testsignal <- list(leftIndex = c(1L, 12L), rightIndex = c(11L, 34L), value = c(0, 1))

  expect_error(computeStat(testy, nq = "s"))
  expect_error(computeStat(testy, nq = c(1L, 2L)))
  expect_error(computeStat(testy, nq = NA))
  expect_error(computeStat(testy, nq = Inf))
  expect_error(computeStat(testy, nq = 0L))
  expect_error(computeStat(testy, nq = -1L))
  expect_error(computeStat(testy, nq = 33L))
  expect_error(computeStat(testy, intervalSystem = "dyaPar", lengths = 64L, nq = 64L))
  expect_identical(computeStat(testy), computeStat(testy, nq = testn))
  expect_identical(computeStat(testy, nq = 100), computeStat(testy, nq = 100L))
  expect_identical(computeStat(testy, nq = 100.5), computeStat(testy, nq = 100L))

  test(testy, left = 1L, right = testn, value = 0, testfamily = "gauss", finallengths = 1:testn,
       testnq = 100L, indices = indicesAll, singleStat = singleStatGauss, sd = 1.23, defaultPenalty = "sqrt")
  test(testy, testsignal = 1, left = 1L, right = testn, value = 1, testfamily = "gauss",
       finallengths = 1:testn, testnq = 64L, indices = indicesAll, singleStat = singleStatGauss,
       sd = 1.23, defaultPenalty = "sqrt")
  test(testy, testsignal = testsignal, left = c(1L, 12L), right = c(11L, 34L), value = c(0, 1),
       testfamily = "gauss", finallengths = 1:testn, testnq = 55L,
       indices = indicesAll, singleStat = singleStatGauss, sd = 1.23, defaultPenalty = "sqrt")

  testlengths <- c(1, 3:5, 10:23)
  test(testy, left = 1L, right = testn, value = 0, testfamily = "gauss", testlengths = testlengths,
       finallengths = testlengths, testnq = 100L, indices = indicesAllLengths,
       singleStat = singleStatGauss, sd = 1.23, defaultPenalty = "sqrt")
  test(testy, testsignal = 1, left = 1L, right = testn, value = 1, testfamily = "gauss",
       testlengths = testlengths, finallengths = testlengths, testnq = 64L,
       indices = indicesAllLengths, singleStat = singleStatGauss, sd = 1.23, defaultPenalty = "sqrt")
  test(testy, testsignal = testsignal, left = c(1L, 12L), right = c(11L, 34L), value = c(0, 1),
       testfamily = "gauss", testlengths = testlengths, finallengths = testlengths, testnq = 45L,
       indices = indicesAllLengths, singleStat = singleStatGauss, sd = 1.23, defaultPenalty = "sqrt")
})

test_that("argument output is tested and works in the default case", {
  testn <- 34L
  testy <- rnorm(testn)

  expect_error(computeStat(testy, output = 1))
  expect_error(computeStat(testy, output = "mu"))
  expect_identical(computeStat(testy, signal = 0), computeStat(testy, signal = 0, output = "l"))
})

test_that("... is tested and works", {
  testn <- 34L
  testy <- rnorm(testn)

  expect_error(computeStat(testy, std = 1))
  expect_error(computeStat(testy, familty = "hsmuce"))
  expect_error(computeStat(testy, covariances = c(1, 0.6)))
  expect_error(computeStat(testy, sd = "s"))
  expect_error(computeStat(testy, sd = c(1, 2)))
  expect_error(computeStat(testy, sd = NA))
  expect_error(computeStat(c(1, 2)))
  expect_error(computeStat(testy, sd = Inf))
  expect_error(computeStat(testy, sd = 0))
  expect_error(computeStat(testy, sd = -0.1))
  expect_identical(computeStat(testy), computeStat(testy, sd = sdrobnorm(testy)))
  expect_identical(computeStat(testy, sd = 1L), computeStat(testy, sd = 1))
})

test_that("family hsmuce works", {
  testn <- 100L
  testy <- rnorm(testn)

  expect_error(computeStat(family = "hsmuce"))
  expect_error(computeStat(numeric(0), family = "hsmuce"))
  expect_identical(computeStat(testy, family = "hsmuce"),
                   computeStat(testy, family = "hsmuce", intervalSystem = "dyaPar", lengths = 2^(1:6),
                               penalty = "none", nq = testn, output = "list"))
  expect_error(computeStat(as.integer(testy), family = "hsmuce"))
  expect_error(computeStat(c(testy, "s"), family = "hsmuce"))
  expect_error(computeStat(c(rnorm(10), NA), family = "hsmuce"))
  expect_error(computeStat(c(rnorm(10), Inf), family = "hsmuce"))

  test(testy, left = 1L, right = testn, value = 0, testfamily = "hsmuce", finallengths = 2^(1:6),
       indices = indicesDyaPar, singleStat = singleStatHsmuce, defaultPenalty = "none", tolerance = 1e-6)

  testn <- 34L
  testy <- c(rnorm(testn / 2), rnorm(testn / 2, 1))
  testsignal <- list(leftIndex = c(1L, 12L), rightIndex = c(11L, 34L), value = c(0, 1))

  expect_equal(computeStat(testy, signal = 0, family = "hsmuce"),
               computeStat(testy, signal = list(leftIndex = 1L, rightIndex = 34L, value = 0),
                           family = "hsmuce"), tolerance = 1e-6)

  test(testy, left = 1L, right = testn, value = 0, testfamily = "hsmuce", finallengths = 2^(1:5),
       indices = indicesDyaPar, singleStat = singleStatHsmuce, defaultPenalty = "none", tolerance = 1e-6)
  test(testy, testsignal = 1, left = 1L, right = testn, value = 1, testfamily = "hsmuce",
       finallengths = 2^(1:5), indices = indicesDyaPar, singleStat = singleStatHsmuce,
       defaultPenalty = "none", tolerance = 1e-6)
  test(testy, testsignal = testsignal, left = c(1L, 12L), right = c(11L, 34L), value = c(0, 1),
       testfamily = "hsmuce", finallengths = 2^(1:5), indices = indicesDyaPar, singleStat = singleStatHsmuce,
       defaultPenalty = "none", tolerance = 1e-6)

  expect_identical(computeStat(testy, family = "hsm"), computeStat(testy, family = "hsmuce"))

  expect_identical(computeStat(testy, signal = 0, family = "hsmuce"),
                   computeStat(testy, signal = 0, family = "hsmuce", intervalSystem = "dyaPar"))

  test(testy, left = 1L, right = testn, value = 0, testfamily = "hsmuce", finallengths = 2:testn,
       testintervalSystem = "all", indices = indicesAll, singleStat = singleStatHsmuce,
       defaultPenalty = "none", tolerance = 1e-6)
  test(testy, testsignal = 1, left = 1L, right = testn, value = 1, testfamily = "hsmuce",
       finallengths = 2:testn, testintervalSystem = "all", indices = indicesAll, singleStat = singleStatHsmuce,
       defaultPenalty = "none", tolerance = 1e-6)
  test(testy, testsignal = testsignal, left = c(1L, 12L), right = c(11L, 34L), value = c(0, 1),
       testfamily = "hsmuce", finallengths = 2:testn, testintervalSystem = "all",
       indices = indicesAll, singleStat = singleStatHsmuce, defaultPenalty = "none", tolerance = 1e-6)

  test(testy, left = 1L, right = testn, value = 0, testfamily = "hsmuce", finallengths = 2^(1:5),
       testintervalSystem = "dyaLen", indices = indicesDyaLen, singleStat = singleStatHsmuce,
       defaultPenalty = "none", tolerance = 1e-6)
  test(testy, testsignal = 1, left = 1L, right = testn, value = 1, testfamily = "hsmuce",
       finallengths = 2^(1:5), testintervalSystem = "dyaLen", indices = indicesDyaLen,
       singleStat = singleStatHsmuce, defaultPenalty = "none", tolerance = 1e-6)
  test(testy, testsignal = testsignal, left = c(1L, 12L), right = c(11L, 34L), value = c(0, 1),
       testfamily = "hsmuce", finallengths = 2^(1:5), testintervalSystem = "dyaLen", indices = indicesDyaLen,
       singleStat = singleStatHsmuce, defaultPenalty = "none", tolerance = 1e-6)

  testn <- 45L
  testy <- c(rnorm(20, 1, 0.23), rnorm(testn - 20, -1, 0.34))
  testsignal <- list(leftIndex = c(1L, 23L), rightIndex = c(22L, 45L), value = c(0, -1))

  expect_error(computeStat(testy, family = "hsmuce", lengths = 1:3))
  expect_error(computeStat(testy, family = "hsmuce", intervalSystem = "all", lengths = 1:45))
  expect_error(computeStat(testy, family = "hsmuce", intervalSystem = "dyaLen", lengths = 1:2))

  testlengths <- c(2, 6:10, 23:45)
  test(testy, left = 1L, right = testn, value = 0, testfamily = "hsmuce", testintervalSystem = "all",
       testlengths = testlengths, finallengths = testlengths, indices = indicesAllLengths,
       singleStat = singleStatHsmuce, defaultPenalty = "none", tolerance = 1e-6)
  test(testy, testsignal = 1, left = 1L, right = testn, value = 1, testfamily = "hsmuce",
       testintervalSystem = "all", testlengths = testlengths, finallengths = testlengths,
       indices = indicesAllLengths, singleStat = singleStatHsmuce, defaultPenalty = "none",
       tolerance = 1e-6)
  test(testy, testsignal = testsignal, left = c(1L, 23L), right = c(22L, 45L), value = c(0, -1),
       testfamily = "hsmuce", testintervalSystem = "all", testlengths = testlengths,
       finallengths = testlengths, indices = indicesAllLengths,
       singleStat = singleStatHsmuce, defaultPenalty = "none", tolerance = 1e-6)

  testlengths <- c(2, 4, 32)
  test(testy, left = 1L, right = testn, value = 0, testfamily = "hsmuce", testlengths = testlengths,
       finallengths = testlengths, indices = indicesDyaLenLengths, singleStat = singleStatHsmuce,
       defaultPenalty = "none", tolerance = 1e-6, testintervalSystem = "dyaLen")
  test(testy, testsignal = 1, left = 1L, right = testn, value = 1, testfamily = "hsmuce",
       testlengths = testlengths, finallengths = testlengths, indices = indicesDyaLenLengths,
       singleStat = singleStatHsmuce, defaultPenalty = "none", tolerance = 1e-6, testintervalSystem = "dyaLen")
  test(testy, testsignal = testsignal, left = c(1L, 23L), right = c(22L, 45L), value = c(0, -1),
       testfamily = "hsmuce", testlengths = testlengths, finallengths = testlengths,
       indices = indicesDyaLenLengths, singleStat = singleStatHsmuce, defaultPenalty = "none",
       tolerance = 1e-6, testintervalSystem = "dyaLen")

  testlengths <- c(4, 8, 16)
  test(testy, left = 1L, right = testn, value = 0, testfamily = "hsmuce", testlengths = testlengths,
       finallengths = testlengths, indices = indicesDyaParLengths,
       singleStat = singleStatHsmuce, defaultPenalty = "none", tolerance = 1e-6)
  test(testy, testsignal = 1, left = 1L, right = testn, value = 1, testfamily = "hsmuce",
       testlengths = testlengths, finallengths = testlengths, indices = indicesDyaParLengths,
       singleStat = singleStatHsmuce, defaultPenalty = "none", tolerance = 1e-6)
  test(testy, testsignal = testsignal, left = c(1L, 23L), right = c(22L, 45L), value = c(0, -1),
       testfamily = "hsmuce", testlengths = testlengths, finallengths = testlengths,
       indices = indicesDyaParLengths, singleStat = singleStatHsmuce, defaultPenalty = "none",
       tolerance = 1e-6)

    expect_identical(computeStat(testy, signal = 0, family = "hsmuce"),
                     computeStat(testy, signal = 0, family = "hsmuce", penalty = "none"))
    expect_identical(computeStat(testy, signal = 0, family = "hsmuce", penalty = "weights"),
                     computeStat(testy, signal = 0, family = "hsmuce", penalty = "none"))

  test(testy, left = 1L, right = testn, value = 0, testfamily = "hsmuce", finallengths = 2^(1:5),
       testnq = 100L, indices = indicesDyaPar, singleStat = singleStatHsmuce, defaultPenalty = "none",
       tolerance = 1e-6)

  expect_error(computeStat(testy, family = "hsmuce", sd = 1))
  expect_error(computeStat(testy, family = "hsmuce", intervalsystem = "all"))

  suppressWarnings(ret <- computeStat(c(1), family = "hsmuce"))
  expect_identical(ret$maximum, -Inf)
  expect_identical(ret$stat, numeric(0))
  expect_identical(ret$lengths, integer(0))
})

test_that("family mDependentPS works", {
  testn <- 30L
  testy <- rnorm(testn)
  testCovariances <- seq(1, 0.2, -0.2)

  expect_error(computeStat(family = "mDependentPS", covariances = testCovariances))
  expect_error(computeStat(numeric(0), family = "mDependentPS", covariances = testCovariances))
  expect_identical(computeStat(testy, family = "mDependentPS", covariances = testCovariances),
                   computeStat(testy, family = "mDependentPS", intervalSystem = "dyaLen", lengths = 2^(0:4),
                               penalty = "sqrt", nq = testn, output = "list", covariances = testCovariances))
  expect_error(computeStat(as.integer(testy), family = "mDependentPS", covariances = testCovariances))
  expect_error(computeStat(c(testy, "s"), family = "mDependentPS", covariances = testCovariances))
  expect_error(computeStat(c(rnorm(10), NA), family = "mDependentPS", covariances = testCovariances))

  test(testy, left = 1L, right = testn, value = 0, testfamily = "mDependentPS", finallengths = 2^(0:4),
       indices = indicesDyaLen, singleStat = singleStatmDependentPS,
       covariances = testCovariances, defaultPenalty = "sqrt")

  testn <- 18L
  testy <- c(rnorm(testn / 2), rnorm(testn / 2, 1))
  testsignal <- list(leftIndex = c(1L, 12L), rightIndex = c(11L, 18L), value = c(0, 1))

  test(testy, left = 1L, right = testn, value = 0, testfamily = "mDependentPS", finallengths = 2^(0:4),
       indices = indicesDyaLen, singleStat = singleStatmDependentPS,
       covariances = testCovariances, defaultPenalty = "sqrt")
  test(testy, testsignal = 1, left = 1L, right = testn, value = 1, testfamily = "mDependentPS",
       finallengths = 2^(0:4), indices = indicesDyaLen, singleStat = singleStatmDependentPS,
       covariances = testCovariances, defaultPenalty = "sqrt")
  test(testy, testsignal = testsignal, left = c(1L, 12L), right = c(11L, testn), value = c(0, 1),
       testfamily = "mDependentPS", finallengths = 2^(0:4), indices = indicesDyaLen,
       singleStat = singleStatmDependentPS, covariances = testCovariances, defaultPenalty = "sqrt")

  expect_identical(computeStat(testy, signal = 0, family = "mDependentPS", covariances = testCovariances),
                   computeStat(testy, signal = 0, family = "mDependentPS", intervalSystem = "dyaLen",
                               covariances = testCovariances))

  test(testy, left = 1L, right = testn, value = 0, testfamily = "mDependentPS", finallengths = 1:testn,
       testintervalSystem = "all", indices = indicesAll, singleStat = singleStatmDependentPS,
       covariances = testCovariances, defaultPenalty = "sqrt")
  test(testy, testsignal = 1, left = 1L, right = testn, value = 1, testfamily = "mDependentPS",
       finallengths = 1:testn, testintervalSystem = "all", indices = indicesAll,
       singleStat = singleStatmDependentPS, covariances = testCovariances, defaultPenalty = "sqrt")
  test(testy, testsignal = testsignal, left = c(1L, 12L), right = c(11L, testn), value = c(0, 1),
       testfamily = "mDependentPS", finallengths = 1:testn, testintervalSystem = "all", indices = indicesAll,
       singleStat = singleStatmDependentPS, covariances = testCovariances, defaultPenalty = "sqrt")

  test(testy, left = 1L, right = testn, value = 0, testfamily = "mDependentPS", finallengths = 2^(0:4),
       testintervalSystem = "dyaPar", indices = indicesDyaPar,
       singleStat = singleStatmDependentPS, covariances = testCovariances, defaultPenalty = "sqrt")
  test(testy, testsignal = 1, left = 1L, right = testn, value = 1, testfamily = "mDependentPS",
       finallengths = 2^(0:4), testintervalSystem = "dyaPar", indices = indicesDyaPar,
       singleStat = singleStatmDependentPS, covariances = testCovariances, defaultPenalty = "sqrt")
  test(testy, testsignal = testsignal, left = c(1L, 12L), right = c(11L, testn), value = c(0, 1),
       testfamily = "mDependentPS", finallengths = 2^(0:4), testintervalSystem = "dyaPar",
       indices = indicesDyaPar, singleStat = singleStatmDependentPS, covariances = testCovariances,
       defaultPenalty = "sqrt")

  testn <- 17L
  testcovariances <- as.numeric(ARMAacf(ar = c(), ma = c(0.8, 0.5, 0.3), lag.max = 3))
  testsignal <- list(leftIndex = c(1L, 13L), rightIndex = c(12L, testn), value = c(0, -1))
  testy <- as.numeric(arima.sim(n = testn, list(ar = c(), ma = c(0.8, 0.5, 0.3)), sd = 1))

  testlengths <- c(1, 3:5, 10:14)
  test(testy, left = 1L, right = testn, value = 0, testfamily = "mDependentPS", testlengths = testlengths,
       finallengths = testlengths, indices = indicesAllLengths, testintervalSystem = "all",
       singleStat = singleStatmDependentPS, covariances = testCovariances, defaultPenalty = "sqrt")
  test(testy, testsignal = 1, left = 1L, right = testn, value = 1, testfamily = "mDependentPS",
       testlengths = testlengths, finallengths = testlengths, testintervalSystem = "all", defaultPenalty = "sqrt",
       indices = indicesAllLengths, singleStat = singleStatmDependentPS, covariances = testCovariances)
  test(testy, testsignal = testsignal, left = c(1L, 13L), right = c(12L, testn), value = c(0, -1),
       testfamily = "mDependentPS", testlengths = testlengths, finallengths = testlengths,
       testintervalSystem = "all", indices = indicesAllLengths, singleStat = singleStatmDependentPS,
       covariances = testCovariances, defaultPenalty = "sqrt")

  testlengths <- c(1, 4, 16)
  test(testy, left = 1L, right = testn, value = 0, testfamily = "mDependentPS", testlengths = testlengths,
       finallengths = testlengths, indices = indicesDyaLenLengths,
       singleStat = singleStatmDependentPS, covariances = testCovariances, defaultPenalty = "sqrt")
  test(testy, testsignal = 1, left = 1L, right = testn, value = 1, testfamily = "mDependentPS",
       testlengths = testlengths, finallengths = testlengths, defaultPenalty = "sqrt",
       indices = indicesDyaLenLengths, singleStat = singleStatmDependentPS, covariances = testCovariances)
  test(testy, testsignal = testsignal, left = c(1L, 13L), right = c(12L, testn), value = c(0, -1),
       testfamily = "mDependentPS", testlengths = testlengths, finallengths = testlengths, defaultPenalty = "sqrt",
       indices = indicesDyaLenLengths, singleStat = singleStatmDependentPS, covariances = testCovariances)

  testlengths <- c(2, 4, 16)
  test(testy, left = 1L, right = testn, value = 0, testfamily = "mDependentPS", testlengths = testlengths,
       finallengths = testlengths, defaultPenalty = "sqrt", indices = indicesDyaParLengths,
       singleStat = singleStatmDependentPS, covariances = testCovariances, testintervalSystem = "dyaPar")
  test(testy, testsignal = 1, left = 1L, right = testn, value = 1, testfamily = "mDependentPS",
       testlengths = testlengths, finallengths = testlengths, indices = indicesDyaParLengths,
       singleStat = singleStatmDependentPS, covariances = testCovariances,
       testintervalSystem = "dyaPar", defaultPenalty = "sqrt")
  test(testy, testsignal = testsignal, left = c(1L, 13L), right = c(12L, testn), value = c(0, -1),
       testfamily = "mDependentPS", testlengths = testlengths, finallengths = testlengths,
       indices = indicesDyaParLengths, singleStat = singleStatmDependentPS, covariances = testCovariances,
       testintervalSystem = "dyaPar", defaultPenalty = "sqrt")

  expect_identical(computeStat(testy, family = "mDependentPS", signal = 0, covariances = testCovariances),
                   computeStat(testy, family = "mDependentPS", signal = 0, covariances = testCovariances,
                               penalty = "sqrt"))

  test(testy, left = 1L, right = testn, value = 0, testfamily = "mDependentPS", finallengths = 2^(0:4),
       testnq = 100L, indices = indicesDyaLen, singleStat = singleStatmDependentPS,
       covariances = testCovariances, defaultPenalty = "sqrt")
})

test_that("arguments in ... are tested and work for family mDependentPS", {
  testn <- 17L
  testcovariances <- as.numeric(ARMAacf(ar = c(), ma = c(0.8, 0.5, 0.3), lag.max = 3)) * 2^2
  testcorrelations <- testcovariances / testcovariances[1]
  testsignal <- list(leftIndex = c(1L, 13L), rightIndex = c(12L, testn), value = c(0, -1))
  testfilter <- list(acf = testcorrelations)
  class(testfilter) <- c("lowpassFilter", class(testfilter))
  testy <- as.numeric(arima.sim(n = testn, list(ar = c(), ma = c(0.8, 0.5, 0.3)), sd = 2))

  expect_error(computeStat(testy, family = "mDependentPS", covariances = testCovariances, std = 1))
  expect_error(computeStat(testy, family = "mDependentPS", correlations = testcorrelations, std = 1))
  expect_error(computeStat(testy, family = "mDependentPS", filter = testfilter, std = 1))
  expect_error(computeStat(testy, family = "mDependentPS", intervalsystem = "all"))
  expect_error(computeStat(testy, family = "mDependentPS"))
  expect_error(computeStat(testy, family = "mDependentPS", sd = 1))
  expect_error(computeStat(testy, family = "mDependentPS", covariances = c(testcovariances, "s")))
  expect_error(computeStat(testy, family = "mDependentPS", covariances = c(testcovariances, NA)))
  expect_error(computeStat(testy, family = "mDependentPS", covariances = c(testcovariances, Inf)))
  expect_error(computeStat(testy, family = "mDependentPS", covariances = c(0.01, testcovariances)))
  expect_error(computeStat(testy, family = "mDependentPS", covariances = c(-10, testcovariances)))
  expect_error(computeStat(testy, family = "mDependentPS", covariances = c(1, testcovariances)))
  expect_error(computeStat(testy, family = "mDependentPS", covariances = c(testcovariances, 0)))

  expect_error(computeStat(testy, family = "mDependentPS", correlations = c(testcorrelations, "s")))
  expect_error(computeStat(testy, family = "mDependentPS", correlations = c(testcorrelations, NA)))
  expect_error(computeStat(testy, family = "mDependentPS", correlations = c(testcorrelations, Inf)))
  expect_error(computeStat(testy, family = "mDependentPS", correlations = c(testcorrelations, 1.1)))
  expect_error(computeStat(testy, family = "mDependentPS", correlations = c(0.99, testcorrelations[-1])))
  expect_error(computeStat(testy, family = "mDependentPS", correlations = c(testcorrelations, 0)))

  expect_error(computeStat(testy, sd = "s", family = "mDependentPS", correlations = testcorrelations))
  expect_error(computeStat(testy, sd = c(1, 2), family = "mDependentPS", correlations = testcorrelations))
  expect_error(computeStat(testy, sd = NA, family = "mDependentPS", correlations = testcorrelations))
  expect_error(computeStat(c(1, 2), family = "mDependentPS", correlations = testcorrelations))
  expect_error(computeStat(testy, sd = Inf, family = "mDependentPS", correlations = testcorrelations))
  expect_error(computeStat(testy, sd = 0, family = "mDependentPS", correlations = testcorrelations))
  expect_error(computeStat(testy, sd = -0.1, family = "mDependentPS", correlations = testcorrelations))

  expect_error(computeStat(testy, family = "mDependentPS",
                           filter = list(param = list(acf = testcorrelations))))

  expect_error(computeStat(testy, sd = "s", family = "mDependentPS", filter = testfilter))
  expect_error(computeStat(testy, sd = c(1, 2), family = "mDependentPS", filter = testfilter))
  expect_error(computeStat(testy, sd = NA, family = "mDependentPS", filter = testfilter))
  expect_error(computeStat(c(1, 2), family = "mDependentPS", filter = testfilter))
  expect_error(computeStat(testy, sd = Inf, family = "mDependentPS", filter = testfilter))
  expect_error(computeStat(testy, sd = 0, family = "mDependentPS", filter = testfilter))
  expect_error(computeStat(testy, sd = -0.1, family = "mDependentPS", filter = testfilter))

  expect_identical(computeStat(testy, family = "mDependentPS",
                               covariances = sdrobnorm(testy, lag = 4)^2 * testcorrelations),
                   computeStat(testy, family = "mDependentPS", correlations = testcorrelations))
  expect_identical(computeStat(testy, family = "mDependentPS",
                               covariances = 1.1^2 * testcorrelations),
                   computeStat(testy, family = "mDependentPS", correlations = testcorrelations, sd = 1.1))
  expect_identical(computeStat(testy, family = "mDependentPS",
                               covariances = sdrobnorm(testy, lag = 4)^2 * testcorrelations),
                   computeStat(testy, family = "mDependentPS", filter = testfilter))
  expect_identical(computeStat(testy, family = "mDependentPS",
                               covariances = 1.1^2 * testcorrelations),
                   computeStat(testy, family = "mDependentPS", filter = testfilter, sd = 1.1))
})

test_that("family jsmurf works", {
  testn <- 70L
  testy <- rnorm(testn)
  testfilter <- lowpassFilter::lowpassFilter(type = "bessel", param = list(pole = 4L, cutoff = 0.1), sr = 1e4)

  expect_error(computeStat(family = "jsmurf", filter = testfilter))
  expect_error(computeStat(numeric(0), family = "jsmurf", filter = testfilter))
  expect_identical(computeStat(testy, family = "jsmurf", filter = testfilter),
                   computeStat(testy, family = "jsmurf", intervalSystem = "dyaLen", lengths = 2^(4:6),
                               penalty = "sqrt", nq = testn, output = "list", filter = testfilter))
  expect_error(computeStat(as.integer(testy), family = "jsmurf", filter = testfilter))
  expect_error(computeStat(c(testy, "s"), family = "jsmurf", filter = testfilter))
  expect_error(computeStat(c(rnorm(10), NA), family = "jsmurf", filter = testfilter))

  test(testy, left = 1L, right = testn, value = 0, testfamily = "jsmurf", finallengths = 2^(4:6),
       indices = indicesDyaLen, singleStat = singleStatJsmurf, sd = 1,
       filter = testfilter, defaultPenalty = "sqrt", shift = - testfilter$len)
  test(testy, left = 1L, right = testn, value = 0, testfamily = "jsmurfPS", finallengths = 2^(4:6),
       indices = indicesDyaLen, singleStat = singleStatJsmurfPS, sd = 1,
       filter = testfilter, defaultPenalty = "sqrt", shift = - testfilter$len)
  test(testy, left = 1L, right = testn, value = 0, testfamily = "jsmurfLR", finallengths = 2^(4:6),
       indices = indicesDyaLen, singleStat = singleStatJsmurfLR, sd = 1,
       filter = testfilter, defaultPenalty = "sqrt", shift = - testfilter$len)

  testn <- 68L
  testy <- c(rnorm(testn / 2), rnorm(testn / 2, 1))
  testsignal <- list(leftIndex = c(1L, 22L), rightIndex = c(21L, 68L), value = c(0, 1))

  test(testy, testsignal = 1.234, left = c(1L), right = c(testn), value = c(1.234),
       testfamily = "jsmurf", finallengths = 2^(4:6), indices = indicesDyaLen, sd = 1,
       singleStat = singleStatJsmurf, filter = testfilter, defaultPenalty = "sqrt", shift = -testfilter$len)
  test(testy, testsignal = 1.234, left = c(1L), right = c(testn), value = c(1.234),
       testfamily = "jsmurfPS", finallengths = 2^(4:6), indices = indicesDyaLen, sd = 1,
       singleStat = singleStatJsmurfPS, filter = testfilter, defaultPenalty = "sqrt", shift = -testfilter$len)
  test(testy, testsignal = 1.234, left = c(1L), right = c(testn), value = c(1.234),
       testfamily = "jsmurfLR", finallengths = 2^(4:6), indices = indicesDyaLen, sd = 1,
       singleStat = singleStatJsmurfLR, filter = testfilter, defaultPenalty = "sqrt", shift = -testfilter$len)

  test(testy, testsignal = testsignal, left = c(1L, 22L), right = c(21L, testn), value = c(0, 1),
       testfamily = "jsmurf", finallengths = 2^(4:6), indices = indicesDyaLen, sd = 1,
       singleStat = singleStatJsmurf, filter = testfilter, defaultPenalty = "sqrt", shift = -testfilter$len)
  test(testy, testsignal = testsignal, left = c(1L, 22L), right = c(21L, testn), value = c(0, 1),
       testfamily = "jsmurfPS", finallengths = 2^(4:6), indices = indicesDyaLen, sd = 1,
       singleStat = singleStatJsmurfPS, filter = testfilter, defaultPenalty = "sqrt",
       shift = -testfilter$len)
  test(testy, testsignal = testsignal, left = c(1L, 22L), right = c(21L, testn), value = c(0, 1),
       testfamily = "jsmurfLR", finallengths = 2^(4:6), indices = indicesDyaLen, sd = 1,
       singleStat = singleStatJsmurfLR, filter = testfilter, defaultPenalty = "sqrt",
       shift = -testfilter$len)

  expect_identical(computeStat(testy, signal = 0, family = "jsmurf", filter = testfilter),
                   computeStat(testy, signal = 0, family = "jsmurf", intervalSystem = "dyaLen",
                               filter = testfilter))

  test(testy, left = 1L, right = testn, value = 0, testfamily = "jsmurf", finallengths = (testfilter$len + 1):testn,
       testintervalSystem = "all", indices = indicesAll, sd = 1, singleStat = singleStatJsmurf,
       filter = testfilter, defaultPenalty = "sqrt", shift = -testfilter$len)
  test(testy, testsignal = 1, left = 1L, right = testn, value = 1, testfamily = "jsmurfPS",
       finallengths = (testfilter$len + 1):testn, testintervalSystem = "all", indices = indicesAll, sd = 1,
       singleStat = singleStatJsmurfPS, filter = testfilter, defaultPenalty = "sqrt", shift = -testfilter$len)
  test(testy, testsignal = testsignal, left = c(1L, 22L), right = c(21L, testn), value = c(0, 1),
       testfamily = "jsmurfLR", finallengths = (testfilter$len + 1):testn,
       testintervalSystem = "all", indices = indicesAll, sd = 1, singleStat = singleStatJsmurfLR,
       filter = testfilter, defaultPenalty = "sqrt", shift = -testfilter$len)

  test(testy, left = 1L, right = testn, value = 0, testfamily = "jsmurf", finallengths = 2^(4:6),
       testintervalSystem = "dyaPar", indices = indicesDyaPar, sd = 1, singleStat = singleStatJsmurf,
       filter = testfilter, defaultPenalty = "sqrt", shift = -testfilter$len)
  test(testy, testsignal = 1, left = 1L, right = testn, value = 1, testfamily = "jsmurfPS",
       finallengths = 2^(4:6), testintervalSystem = "dyaPar", indices = indicesDyaPar, sd = 1,
       singleStat = singleStatJsmurfPS, filter = testfilter, defaultPenalty = "sqrt", shift = -testfilter$len)
  test(testy, testsignal = testsignal, left = c(1L, 22L), right = c(21L, testn), value = c(0, 1),
       testfamily = "jsmurfLR", finallengths = 2^(4:6),
       testintervalSystem = "dyaPar", indices = indicesDyaPar, sd = 1, singleStat = singleStatJsmurfLR,
       filter = testfilter, defaultPenalty = "sqrt", shift = -testfilter$len)

  expect_error(computeStat(y = testy, family = "jsmurf", filter = testfilter, lengths = c(8L, 16L)))
  expect_error(computeStat(y = testy, family = "jsmurfPS", filter = testfilter, lengths = c(8L, 16L)))
  expect_error(computeStat(y = testy, family = "jsmurfLR", filter = testfilter, lengths = c(8L, 16L)))

  testlengths <- c(16L, 64L)
  test(testy, left = 1L, right = testn, value = 0, testfamily = "jsmurfLR", finallengths = testlengths,
       indices = indicesDyaLen, testlengths = testlengths, singleStat = singleStatJsmurfLR, sd = 1,
       filter = testfilter, defaultPenalty = "sqrt", shift = - testfilter$len)
  test(testy, testsignal = 1.234, left = c(1L), right = c(testn), value = c(1.234), testlengths = testlengths,
       testfamily = "jsmurf", finallengths = testlengths, indices = indicesDyaLen, sd = 1,
       singleStat = singleStatJsmurf, filter = testfilter, defaultPenalty = "sqrt", shift = -testfilter$len)
  test(testy, testsignal = testsignal, left = c(1L, 22L), right = c(21L, testn), value = c(0, 1),
       testfamily = "jsmurfPS", finallengths = testlengths, indices = indicesDyaLen, sd = 1,
       singleStat = singleStatJsmurfPS, filter = testfilter, defaultPenalty = "sqrt",
       testlengths = testlengths, shift = -testfilter$len)

  testlengths <- c(13L, 15:18, 64L)
  test(testy, left = 1L, right = testn, value = 0, testfamily = "jsmurfLR", finallengths = testlengths,
       testintervalSystem = "all", indices = indicesAll, testlengths = testlengths, singleStat = singleStatJsmurfLR,
       sd = 1, filter = testfilter, defaultPenalty = "sqrt", shift = - testfilter$len)
  test(testy, testsignal = 1.234, left = c(1L), right = c(testn), value = c(1.234), testlengths = testlengths,
       testintervalSystem = "all", testfamily = "jsmurf", finallengths = testlengths, indices = indicesAll, sd = 1,
       singleStat = singleStatJsmurf, filter = testfilter, defaultPenalty = "sqrt", shift = -testfilter$len)
  test(testy, testsignal = testsignal, left = c(1L, 22L), right = c(21L, testn), value = c(0, 1),
       testfamily = "jsmurfPS", finallengths = testlengths, indices = indicesAll, sd = 1,
       testintervalSystem = "all", singleStat = singleStatJsmurfPS, filter = testfilter, defaultPenalty = "sqrt",
       testlengths = testlengths, shift = -testfilter$len)

  testlengths <- c(32L, 64L)
  test(testy, left = 1L, right = testn, value = 0, testfamily = "jsmurfLR", finallengths = testlengths,
       testintervalSystem = "dyaPar", indices = indicesDyaPar, testlengths = testlengths,
       singleStat = singleStatJsmurfLR, tolerance = 1e-10,
       sd = 1, filter = testfilter, defaultPenalty = "sqrt", shift = - testfilter$len)
  test(testy, testsignal = 1.234, left = c(1L), right = c(testn), value = c(1.234), testlengths = testlengths,
       testintervalSystem = "dyaPar", testfamily = "jsmurf", finallengths = testlengths, indices = indicesDyaPar,
       sd = 1, singleStat = singleStatJsmurf, filter = testfilter, defaultPenalty = "sqrt", shift = -testfilter$len)
  test(testy, testsignal = testsignal, left = c(1L, 22L), right = c(21L, testn), value = c(0, 1),
       testfamily = "jsmurfPS", finallengths = testlengths, indices = indicesDyaPar, sd = 1,
       testintervalSystem = "dyaPar", singleStat = singleStatJsmurfPS, filter = testfilter, defaultPenalty = "sqrt",
       testlengths = testlengths, shift = -testfilter$len)

  expect_identical(computeStat(testy, family = "jsmurf", signal = 0, filter = testfilter),
                   computeStat(testy, family = "jsmurf", signal = 0, filter = testfilter,
                               penalty = "sqrt"))

  test(testy, left = 1L, right = testn, value = 0, testfamily = "jsmurf", finallengths = 2^(4:6),
       testnq = 100L, indices = indicesDyaLen, singleStat = singleStatJsmurf, sd = 1,
       filter = testfilter, defaultPenalty = "sqrt", shift = -testfilter$len)
  test(testy, left = 1L, right = testn, value = 0, testfamily = "jsmurfPS", finallengths = 2^(4:6),
       testnq = 100L, indices = indicesDyaLen, singleStat = singleStatJsmurfPS, sd = 1,
       filter = testfilter, defaultPenalty = "sqrt", shift = -testfilter$len)
  test(testy, left = 1L, right = testn, value = 0, testfamily = "jsmurfLR", finallengths = 2^(4:6),
       testnq = 100L, indices = indicesDyaLen, singleStat = singleStatJsmurfLR, sd = 1,
       filter = testfilter, defaultPenalty = "sqrt", shift = -testfilter$len)

  suppressWarnings(ret <- computeStat(rnorm(10), family = "jsmurfLR", filter = testfilter, sd = 1))
  expect_identical(ret$maximum, -Inf)
  expect_identical(ret$stat, numeric(0))
  expect_identical(ret$lengths, integer(0))
})

test_that("arguments in ... are tested and work for family jsmurf", {
  testn <- 35L
  testy <- rnorm(testn)
  testfilter <- lowpassFilter::lowpassFilter(type = "bessel", param = list(pole = 4L, cutoff = 0.1), sr = 1e4)

  expect_error(computeStat(testy, family = "jsmurf", filter = testfilter, std = 1))
  expect_error(computeStat(testy, family = "jsmurf", filter = testfilter, intervalsystem = "all"))
  expect_error(computeStat(testy, family = "jsmurf"))
  expect_error(computeStat(testy, family = "jsmurf", sd = 1))

  expect_error(computeStat(testy, family = "jsmurf", filter = unclass(testfilter)))

  expect_error(computeStat(testy, sd = "s", family = "jsmurf", filter = testfilter))
  expect_error(computeStat(testy, sd = c(1, 2), family = "jsmurf", filter = testfilter))
  expect_error(computeStat(testy, sd = NA, family = "jsmurf", filter = testfilter))
  expect_error(computeStat(c(1, 2), family = "jsmurf", filter = testfilter))
  expect_error(computeStat(testy, sd = Inf, family = "jsmurf", filter = testfilter))
  expect_error(computeStat(testy, sd = 0, family = "jsmurf", filter = testfilter))
  expect_error(computeStat(testy, sd = -0.1, family = "jsmurf", filter = testfilter))

  expect_identical(computeStat(testy, family = "jsmurf", filter = testfilter),
                   computeStat(testy, family = "jsmurf", filter = testfilter,
                               sd = sdrobnorm(testy, lag = testfilter$len + 1)))


  expect_error(computeStat(testy, family = "jsmurfPS", filter = testfilter, std = 1))
  expect_error(computeStat(testy, family = "jsmurfPS", filter = testfilter, intervalsystem = "all"))
  expect_error(computeStat(testy, family = "jsmurfPS"))
  expect_error(computeStat(testy, family = "jsmurfPS", sd = 1))

  expect_error(computeStat(testy, family = "jsmurfPS", filter = unclass(testfilter)))

  expect_error(computeStat(testy, sd = "s", family = "jsmurfPS", filter = testfilter))
  expect_error(computeStat(testy, sd = c(1, 2), family = "jsmurfPS", filter = testfilter))
  expect_error(computeStat(testy, sd = NA, family = "jsmurfPS", filter = testfilter))
  expect_error(computeStat(c(1, 2), family = "jsmurfPS", filter = testfilter))
  expect_error(computeStat(testy, sd = Inf, family = "jsmurfPS", filter = testfilter))
  expect_error(computeStat(testy, sd = 0, family = "jsmurfPS", filter = testfilter))
  expect_error(computeStat(testy, sd = -0.1, family = "jsmurfPS", filter = testfilter))

  expect_identical(computeStat(testy, family = "jsmurfPS", filter = testfilter),
                   computeStat(testy, family = "jsmurfPS", filter = testfilter,
                               sd = sdrobnorm(testy, lag = testfilter$len + 1)))


  expect_error(computeStat(testy, family = "jsmurfLR", filter = testfilter, std = 1))
  expect_error(computeStat(testy, family = "jsmurfLR", filter = testfilter, intervalsystem = "all"))
  expect_error(computeStat(testy, family = "jsmurfLR"))
  expect_error(computeStat(testy, family = "jsmurfLR", sd = 1))

  expect_error(computeStat(testy, family = "jsmurfLR", filter = unclass(testfilter)))

  expect_error(computeStat(testy, sd = "s", family = "jsmurfLR", filter = testfilter))
  expect_error(computeStat(testy, sd = c(1, 2), family = "jsmurfLR", filter = testfilter))
  expect_error(computeStat(testy, sd = NA, family = "jsmurfLR", filter = testfilter))
  expect_error(computeStat(c(1, 2), family = "jsmurfLR", filter = testfilter))
  expect_error(computeStat(testy, sd = Inf, family = "jsmurfLR", filter = testfilter))
  expect_error(computeStat(testy, sd = 0, family = "jsmurfLR", filter = testfilter))
  expect_error(computeStat(testy, sd = -0.1, family = "jsmurfLR", filter = testfilter))

  expect_identical(computeStat(testy, family = "jsmurfLR", filter = testfilter),
                   computeStat(testy, family = "jsmurfLR", filter = testfilter,
                               sd = sdrobnorm(testy, lag = testfilter$len + 1)))


  testfilter <- lowpassFilter::lowpassFilter(type = "bessel", param = list(pole = 6L, cutoff = 0.25),
                                             sr = 5e4, shift = 0.2)
  test(testy, left = 1L, right = testn, value = 0, testfamily = "jsmurf", finallengths = 2^(2:5),
       testnq = 100L, indices = indicesDyaLen, singleStat = singleStatJsmurf, sd = 1,
       filter = testfilter, defaultPenalty = "sqrt", shift = -testfilter$len)
  test(testy, left = 1L, right = testn, value = 0, testfamily = "jsmurfPS", finallengths = 2^(2:5),
       testnq = 100L, indices = indicesDyaLen, singleStat = singleStatJsmurfPS, sd = 1,
       filter = testfilter, defaultPenalty = "sqrt", shift = -testfilter$len)
  test(testy, left = 1L, right = testn, value = 0, testfamily = "jsmurfLR", finallengths = 2^(2:5),
       testnq = 100L, indices = indicesDyaLen, singleStat = singleStatJsmurfLR, sd = 1,
       filter = testfilter, defaultPenalty = "sqrt", shift = -testfilter$len)
})

test_that("family hjsmurf works", {
  testn <- 70L
  testy <- rnorm(testn)
  testfilter <- lowpassFilter::lowpassFilter(type = "bessel", param = list(pole = 4L, cutoff = 0.1), sr = 1e4)

  expect_error(computeStat(family = "hjsmurf", filter = testfilter))
  expect_error(computeStat(numeric(0), family = "hjsmurf", filter = testfilter))
  expect_identical(computeStat(testy, family = "hjsmurf", filter = testfilter),
                   computeStat(testy, family = "hjsmurf", intervalSystem = "dyaLen", lengths = 2^(4:6),
                               penalty = "none", nq = testn, output = "list", filter = testfilter))
  expect_error(computeStat(as.integer(testy), family = "hjsmurf", filter = testfilter))
  expect_error(computeStat(c(testy, "s"), family = "hjsmurf", filter = testfilter))
  expect_error(computeStat(c(rnorm(10), NA), family = "hjsmurf", filter = testfilter))

  test(testy, left = 1L, right = testn, value = 0, testfamily = "hjsmurf", finallengths = 2^(4:6),
       indices = indicesDyaLen, singleStat = singleStatHjsmurf,
       filter = testfilter, defaultPenalty = "none", shift = -testfilter$len)
  test(testy, left = 1L, right = testn, value = 0, testfamily = "hjsmurfSPS", finallengths = 2^(4:6),
       indices = indicesDyaLen, singleStat = singleStatHjsmurfSPS,
       filter = testfilter, defaultPenalty = "none", shift = -testfilter$len)
  test(testy, left = 1L, right = testn, value = 0, testfamily = "hjsmurfLR", finallengths = 2^(4:6),
       indices = indicesDyaLen, singleStat = singleStatHjsmurfLR,
       filter = testfilter, defaultPenalty = "none", shift = -testfilter$len)

  testn <- 68L
  testy <- c(rnorm(testn / 2), rnorm(testn / 2, 1))
  testsignal <- list(leftIndex = c(1L, 22L), rightIndex = c(21L, 68L), value = c(0, 1))

  test(testy, testsignal = 1.234, left = c(1L), right = c(testn), value = c(1.234),
       testfamily = "hjsmurf", finallengths = 2^(4:6), indices = indicesDyaLen,
       singleStat = singleStatHjsmurf, filter = testfilter, defaultPenalty = "none", shift = -testfilter$len)
  test(testy, testsignal = 1.234, left = c(1L), right = c(testn), value = c(1.234),
       testfamily = "hjsmurfSPS", finallengths = 2^(4:6), indices = indicesDyaLen,
       singleStat = singleStatHjsmurfSPS, filter = testfilter, defaultPenalty = "none", shift = -testfilter$len)
  test(testy, testsignal = 1.234, left = c(1L), right = c(testn), value = c(1.234),
       testfamily = "hjsmurfLR", finallengths = 2^(4:6), indices = indicesDyaLen,
       singleStat = singleStatHjsmurfLR, filter = testfilter, defaultPenalty = "none", shift = -testfilter$len)

  test(testy, testsignal = testsignal, left = c(1L, 22L), right = c(21L, testn), value = c(0, 1),
       testfamily = "hjsmurf", finallengths = 2^(4:6), indices = indicesDyaLen,
       singleStat = singleStatHjsmurf, filter = testfilter, defaultPenalty = "none", shift = -testfilter$len)
  test(testy, testsignal = testsignal, left = c(1L, 22L), right = c(21L, testn), value = c(0, 1),
       testfamily = "hjsmurfSPS", finallengths = 2^(4:6), indices = indicesDyaLen,
       singleStat = singleStatHjsmurfSPS, filter = testfilter, defaultPenalty = "none",
       shift = -testfilter$len)
  test(testy, testsignal = testsignal, left = c(1L, 22L), right = c(21L, testn), value = c(0, 1),
       testfamily = "hjsmurfLR", finallengths = 2^(4:6), indices = indicesDyaLen,
       singleStat = singleStatHjsmurfLR, filter = testfilter, defaultPenalty = "none",
       shift = -testfilter$len)

  expect_identical(computeStat(testy, signal = 0, family = "hjsmurf", filter = testfilter),
                   computeStat(testy, signal = 0, family = "hjsmurf", intervalSystem = "dyaLen",
                               filter = testfilter))

  test(testy, left = 1L, right = testn, value = 0, testfamily = "hjsmurf", finallengths = (testfilter$len + 2):testn,
       testintervalSystem = "all", indices = indicesAll, singleStat = singleStatHjsmurf,
       filter = testfilter, defaultPenalty = "none", shift = -testfilter$len, tolerance = 1e-6)
  test(testy, testsignal = 1, left = 1L, right = testn, value = 1, testfamily = "hjsmurfSPS",
       finallengths = (testfilter$len + 2):testn, testintervalSystem = "all", indices = indicesAll,
       singleStat = singleStatHjsmurfSPS, filter = testfilter, defaultPenalty = "none", shift = -testfilter$len,
       tolerance = 1e-6)
  test(testy, testsignal = testsignal, left = c(1L, 22L), right = c(21L, testn), value = c(0, 1),
       testfamily = "hjsmurfLR", finallengths = (testfilter$len + 2):testn,
       testintervalSystem = "all", indices = indicesAll, singleStat = singleStatHjsmurfLR,
       filter = testfilter, defaultPenalty = "none", shift = -testfilter$len, tolerance = 1e-6)

  test(testy, left = 1L, right = testn, value = 0, testfamily = "hjsmurf", finallengths = 2^(4:6),
       testintervalSystem = "dyaPar", indices = indicesDyaPar, singleStat = singleStatHjsmurf,
       filter = testfilter, defaultPenalty = "none", shift = -testfilter$len)
  test(testy, testsignal = 1, left = 1L, right = testn, value = 1, testfamily = "hjsmurfSPS",
       finallengths = 2^(4:6), testintervalSystem = "dyaPar", indices = indicesDyaPar,
       singleStat = singleStatHjsmurfSPS, filter = testfilter, defaultPenalty = "none", shift = -testfilter$len)
  test(testy, testsignal = testsignal, left = c(1L, 22L), right = c(21L, testn), value = c(0, 1),
       testfamily = "hjsmurfLR", finallengths = 2^(4:6),
       testintervalSystem = "dyaPar", indices = indicesDyaPar, singleStat = singleStatHjsmurfLR,
       filter = testfilter, defaultPenalty = "none", shift = -testfilter$len)

  expect_error(computeStat(y = testy, family = "hjsmurf", filter = testfilter, lengths = c(8L, 16L)))
  expect_error(computeStat(y = testy, family = "hjsmurfPS", filter = testfilter, lengths = c(8L, 16L)))
  expect_error(computeStat(y = testy, family = "hjsmurfLR", filter = testfilter, lengths = c(8L, 16L)))

  testlengths <- c(16L, 64L)
  test(testy, left = 1L, right = testn, value = 0, testfamily = "hjsmurfLR", finallengths = testlengths,
       indices = indicesDyaLen, testlengths = testlengths, singleStat = singleStatHjsmurfLR,
       filter = testfilter, defaultPenalty = "none", shift = - testfilter$len)
  test(testy, testsignal = 1.234, left = c(1L), right = c(testn), value = c(1.234), testlengths = testlengths,
       testfamily = "hjsmurf", finallengths = testlengths, indices = indicesDyaLen,
       singleStat = singleStatHjsmurf, filter = testfilter, defaultPenalty = "none", shift = -testfilter$len)
  test(testy, testsignal = testsignal, left = c(1L, 22L), right = c(21L, testn), value = c(0, 1),
       testfamily = "hjsmurfSPS", finallengths = testlengths, indices = indicesDyaLen,
       singleStat = singleStatHjsmurfSPS, filter = testfilter, defaultPenalty = "none",
       testlengths = testlengths, shift = -testfilter$len)

  testlengths <- c(13L, 15:18, 64L)
  test(testy, left = 1L, right = testn, value = 0, testfamily = "hjsmurfLR", finallengths = testlengths,
       testintervalSystem = "all", indices = indicesAll, testlengths = testlengths, singleStat = singleStatHjsmurfLR,
       filter = testfilter, defaultPenalty = "none", shift = - testfilter$len, tolerance = 1e-6)
  test(testy, testsignal = 1.234, left = c(1L), right = c(testn), value = c(1.234), testlengths = testlengths,
       testintervalSystem = "all", testfamily = "hjsmurf", finallengths = testlengths, indices = indicesAll,
       singleStat = singleStatHjsmurf, filter = testfilter, defaultPenalty = "none", shift = -testfilter$len,
       tolerance = 1e-6)
  test(testy, testsignal = testsignal, left = c(1L, 22L), right = c(21L, testn), value = c(0, 1),
       testfamily = "hjsmurfSPS", finallengths = testlengths, indices = indicesAll,
       testintervalSystem = "all", singleStat = singleStatHjsmurfSPS, filter = testfilter, defaultPenalty = "none",
       testlengths = testlengths, shift = -testfilter$len, tolerance = 1e-6)

  testlengths <- c(32L, 64L)
  test(testy, left = 1L, right = testn, value = 0, testfamily = "hjsmurfLR", finallengths = testlengths,
       testintervalSystem = "dyaPar", indices = indicesDyaPar, testlengths = testlengths,
       singleStat = singleStatHjsmurfLR, filter = testfilter, defaultPenalty = "none", shift = - testfilter$len)
  test(testy, testsignal = 1.234, left = c(1L), right = c(testn), value = c(1.234), testlengths = testlengths,
       testintervalSystem = "dyaPar", testfamily = "hjsmurf", finallengths = testlengths, indices = indicesDyaPar,
       singleStat = singleStatHjsmurf, filter = testfilter, defaultPenalty = "none", shift = -testfilter$len)
  test(testy, testsignal = testsignal, left = c(1L, 22L), right = c(21L, testn), value = c(0, 1),
       testfamily = "hjsmurfSPS", finallengths = testlengths, indices = indicesDyaPar,
       testintervalSystem = "dyaPar", singleStat = singleStatHjsmurfSPS, filter = testfilter, defaultPenalty = "none",
       testlengths = testlengths, shift = -testfilter$len)

  expect_identical(computeStat(testy, family = "hjsmurf", signal = 0, filter = testfilter),
                   computeStat(testy, family = "hjsmurf", signal = 0, filter = testfilter,
                               penalty = "none"))

  test(testy, left = 1L, right = testn, value = 0, testfamily = "hjsmurf", finallengths = 2^(4:6),
       testnq = 100L, indices = indicesDyaLen, singleStat = singleStatHjsmurf,
       filter = testfilter, defaultPenalty = "none", shift = -testfilter$len)
  test(testy, left = 1L, right = testn, value = 0, testfamily = "hjsmurfSPS", finallengths = 2^(4:6),
       testnq = 100L, indices = indicesDyaLen, singleStat = singleStatHjsmurfSPS,
       filter = testfilter, defaultPenalty = "none", shift = -testfilter$len)
  test(testy, left = 1L, right = testn, value = 0, testfamily = "hjsmurfLR", finallengths = 2^(4:6),
       testnq = 100L, indices = indicesDyaLen, singleStat = singleStatHjsmurfLR,
       filter = testfilter, defaultPenalty = "none", shift = -testfilter$len)

  suppressWarnings(ret <- computeStat(rnorm(10), family = "hjsmurfLR", filter = testfilter))
  expect_identical(ret$maximum, -Inf)
  expect_identical(ret$stat, numeric(0))
  expect_identical(ret$lengths, integer(0))
})

test_that("arguments in ... are tested and work for family hjsmurf", {
  testn <- 35L
  testy <- rnorm(testn)
  testfilter <- lowpassFilter::lowpassFilter(type = "bessel", param = list(pole = 4L, cutoff = 0.1), sr = 1e4)

  expect_error(computeStat(testy, family = "hjsmurf", filter = testfilter, sd = 1))
  expect_error(computeStat(testy, family = "hjsmurf", filter = testfilter, intervalsystem = "all"))
  expect_error(computeStat(testy, family = "hjsmurf"))
  expect_error(computeStat(testy, family = "hjsmurf", filter = unclass(testfilter)))

  expect_error(computeStat(testy, family = "hjsmurfSPS", filter = testfilter, sd = 1))
  expect_error(computeStat(testy, family = "hjsmurfSPS", filter = testfilter, intervalsystem = "all"))
  expect_error(computeStat(testy, family = "hjsmurfSPS"))
  expect_error(computeStat(testy, family = "hjsmurfSPS", filter = unclass(testfilter)))


  expect_error(computeStat(testy, family = "hjsmurfLR", filter = testfilter, sd = 1))
  expect_error(computeStat(testy, family = "hjsmurfLR", filter = testfilter, intervalsystem = "all"))
  expect_error(computeStat(testy, family = "hjsmurfLR"))
  expect_error(computeStat(testy, family = "hjsmurfLR", filter = unclass(testfilter)))
})

context("local tests statistics in computeStat")

testJump <- function(grid, time, data, filter, correlations, leftValue, rightValue) {
  len <- length(data)
  m <- min(len, length(correlations) - 1L)

  if (len == 1) {
    A <- matrix(correlations[1], 1, 1)
  } else {
    A <- matrix(0, len, len)
    for (i in 1:(len - 1)) {
      A[i, i] <- correlations[1]
      A[i, i + 1:min(m, len - i)] <- correlations[2:min(m + 1, len - i + 1)]
      A[i + 1:min(m, len - i), i] <- correlations[2:min(m + 1, len - i + 1)]
    }
    A[len, len] <- correlations[1]
  }

  costs <- numeric(length(grid))

  for (i in seq(along = grid)) {
    mu <- lowpassFilter::getConvolutionJump(time, grid[i], leftValue, rightValue, filter)
    costs[i] <- sum((data - mu) * solve(A, data - mu))
  }

  grid[which.min(costs)]
}

statAll <- function(y, filter, fit, singleStat, lengths, add, ...) {
  stat <- rep(-Inf, length(lengths))
  y[add] <- NA

  for (inLen in seq(along = lengths)) {
    len <- lengths[inLen]

    start <- fit$leftEnd[1] + filter$len - 1L
    end <- fit$rightEnd[1] - len - filter$len + 1L
    if (start <= end) {
      for (li in start:end) {
        ri <- li + len
        obs <- y[(li + 1):(ri + filter$len - 1)]
        time <- (li + 1):(ri + filter$len - 1) / filter$sr
        newStat <- singleStat(obs = obs, time = time, filter = filter,
                              left = li / filter$sr, right = ri / filter$sr,
                              leftValue = fit$value[1], rightValue = fit$value[1],
                              leftVar = fit$var[1], rightVar = fit$var[1], cp = 0, ...)
        if (!is.na(newStat)) {
          stat[inLen] <- max(stat[inLen], newStat)
        }
      }
    }

    for (inSeg in seq(along = fit$leftEnd)[-1]) {
      start <- max(fit$rightEnd[inSeg - 1L] - len - filter$len + 2L, start)
      end <- min(fit$leftEnd[inSeg] + filter$len - 2L, fit$rightEnd[inSeg] - len - filter$len + 1L)

      if (start <= end) {
        for (li in start:end) {
          ri <- li + len
          obs <- y[(li + 1):(ri + filter$len - 1)]
          time <- (li + 1):(ri + filter$len - 1) / filter$sr
          newStat <- singleStat(obs = obs, time = time, filter = filter,
                                left = li / filter$sr, right = ri / filter$sr,
                                leftValue = fit$value[inSeg - 1L], rightValue = fit$value[inSeg],
                                leftVar = fit$var[inSeg - 1L], rightVar = fit$var[inSeg],
                                cp = (fit$leftEnd[inSeg] - 1L) / filter$sr, ...)

          if (!is.na(newStat)) {
            stat[inLen] <- max(stat[inLen], newStat)
          }
        }
      }

      start <- fit$leftEnd[inSeg] + filter$len - 1L
      end <- fit$rightEnd[inSeg] - len - filter$len + 1L

      if (start <= end) {
        for (li in start:end) {
          ri <- li + len
          obs <- y[(li + 1):(ri + filter$len - 1)]
          time <- (li + 1):(ri + filter$len - 1) / filter$sr
          newStat <- singleStat(obs = obs, time = time, filter = filter,
                                left = li / filter$sr, right = ri / filter$sr,
                                leftValue = fit$value[inSeg], rightValue = fit$value[inSeg],
                                leftVar = fit$var[inSeg], rightVar = fit$var[inSeg], cp = 0, ...)

          if (!is.na(newStat)) {
            stat[inLen] <- max(stat[inLen], newStat)
          }
        }
      }
    }
  }

  stat
}

test_that("family LR works", {
  testn <- 50L
  testfilter <- lowpassFilter::lowpassFilter(type = "bessel", param = list(pole = 4, cutoff = 0.1), sr = 1e4)
  testy <- lowpassFilter::randomGenerationMA(n = testn, filter = testfilter, seed = "no")
  testfit <-  stepblock(0, leftEnd = 0, rightEnd = testn / testfilter$sr, x0 = 0)

  expect_error(computeStat(family = "LR", filter = testfilter, fit = testfit))
  expect_error(computeStat(y = numeric(0), family = "LR", filter = testfilter, fit = testfit))
  expect_error(computeStat(y = as.integer(testy), family = "LR", filter = testfilter, fit = testfit))
  expect_error(computeStat(y = c(testy[-1], "s"), family = "LR", filter = testfilter, fit = testfit))
  expect_error(computeStat(y = c(testy[-1], NA), family = "LR", filter = testfilter, fit = testfit))

  expect_identical(computeStat(testy, family = "LR", filter = testfilter, fit = testfit),
                   computeStat(testy, family = "LR", intervalSystem = "all", lengths = 1:20,
                               penalty = "none", nq = testn, output = "list", filter = testfilter, fit = testfit))

  ret <- computeStat(y = testy, filter = testfilter, fit = testfit, family = "LR")
  testcor <- testfilter$acf
  testcor[1] <- 2
  comparefit <- stepblock(stats::median(testy[11:39]), leftEnd = 1, rightEnd = c(testn), x0 = 0)
  compare <- statAll(y = testy, filter = testfilter, fit = comparefit, singleStat = singleStatLR,
                     lengths = 1:20, add = integer(0), sd = sdrobnorm(testy, lag = testfilter$len + 1L), regu = 1)
  expect_equal(ret$stat, compare, tolerance = 1e-12)

  expect_identical(computeStat(testy, family = "LR", filter = testfilter, lengths = c(1, 2, 4, 8, 16)),
                   computeStat(testy, family = "LR", intervalSystem = "dyaLen", filter = testfilter))

  ret <- computeStat(y = testy, filter = testfilter, fit = testfit, family = "LR", lengths = c(1:3, 5, 25))
  expect_identical(ret, computeStat(y = testy, filter = testfilter, family = "LR", lengths = c(1:3, 5, 25)))
  testcor <- testfilter$acf
  testcor[1] <- 2
  comparefit <- stepblock(stats::median(testy[11:39]), leftEnd = c(1), rightEnd = c(testn), x0 = 0)
  compare <- statAll(y = testy, filter = testfilter, fit = comparefit, singleStat = singleStatLR,
                     lengths = c(1:3, 5, 25), add = integer(0), sd = sdrobnorm(testy, lag = testfilter$len + 1L),
                     regu = 1)
  expect_equal(ret$stat, compare, tolerance = 1e-12)

  expect_warning(ret2 <- computeStat(y = testy, filter = testfilter, signal = testfit,
                                     fit = testfit, family = "LR", lengths = c(1:3, 5, 25)))
  expect_identical(ret, ret2)

  expect_error(computeStat(y = testy, filter = testfilter, intervalSystem = "dyaPar",
                           fit = testfit, family = "LR", lengths = c(1, 4, 16)))

  expect_identical(computeStat(y = testy, filter = testfilter, family = "LR", lengths = c(1:3, 5, 25),
                               output = "maximum"), ret$maximum)
  expect_identical(computeStat(y = testy, filter = testfilter, family = "LR", lengths = c(1:3, 5, 25),
                               output = "vector"), ret$stat)

  ret <- computeStat(y = testy, filter = testfilter, fit = testfit, family = "LR", lengths = c(1:3, 5, 25),
                     penalty = "sqrt", nq = 100L, output = "maximum")
  expect_equal(ret, max(sqrt(2 * compare) - sqrt(2 * log(exp(1) * 100L / c(1:3, 5, 25)))))

  suppressWarnings(ret <- computeStat(y = rnorm(20), filter = testfilter, family = "LR"))
  expect_identical(ret$maximum, -Inf)
  expect_identical(ret$stat, rep(-Inf, 20))
  expect_identical(ret$lengths, 1:20)

  testn <- 100L
  testfit <-  stepblock(c(0, 1), leftEnd = c(0, 50 / testfilter$sr),
                        rightEnd = c(50 / testfilter$sr, testn / testfilter$sr), x0 = 0)
  testy <- lowpassFilter::randomGenerationMA(n = testn, filter = testfilter, seed = "no")

  ret <- computeStat(y = testy, filter = testfilter, fit = testfit, family = "LR")
  testcor <- testfilter$acf
  testcor[1] <- 2
  j <- as.integer(testJump(39:50 / testfilter$sr, 40:60 / testfilter$sr, testy[40:60], testfilter, testcor,
                           stats::median(testy[11:39]), stats::median(testy[61:89])) * testfilter$sr + 1e-6)
  comparefit <- stepblock(c(stats::median(testy[11:39]), stats::median(testy[61:89])),
                          leftEnd = c(0, j), rightEnd = c(j, 100), x0 = 0)
  comparefit <- stepblock(c(stats::median(testy[11:39]), stats::median(testy[61:89])),
                          leftEnd = c(1, j + 1), rightEnd = c(j, 100), x0 = 0)
  compare <- statAll(y = testy, filter = testfilter, fit = comparefit, singleStat = singleStatLR,
                     lengths = 1:20, add = integer(0), sd = sdrobnorm(testy, lag = testfilter$len + 1L), regu = 1)
  expect_equal(ret$stat, compare, tolerance = 1e-12)

  testfit <-  stepblock(c(0, 1, 0), leftEnd = c(0, 50 / testfilter$sr, 53 / testfilter$sr),
                        rightEnd = c(50 / testfilter$sr, 53 / testfilter$sr, testn / testfilter$sr), x0 = 0)
  ret <- computeStat(y = testy, filter = testfilter, fit = testfit, family = "LR")
  comparefit <- stepblock(c(stats::median(testy[11:39]), NA, stats::median(testy[64:89])),
                          leftEnd = c(1, 51, 54), rightEnd = c(50, 53, 100), x0 = 0)
  compare <- statAll(y = testy, filter = testfilter, fit = comparefit, singleStat = singleStatLR,
                     lengths = 1:20, add = 39:53, sd = sdrobnorm(testy, lag = testfilter$len + 1L), regu = 1)
  expect_equal(ret$stat, compare, tolerance = 1e-12)


  testfit <-  stepblock(c(0, 1, 0, 1, 0), leftEnd = c(0, 50:53 / testfilter$sr),
                        rightEnd = c(50:53 / testfilter$sr, testn / testfilter$sr), x0 = 0)
  expect_warning(ret <- computeStat(y = testy, filter = testfilter, fit = testfit, family = "LR"))
  comparefit <- stepblock(c(stats::median(testy[11:39]), NA, stats::median(testy[64:89])),
                          leftEnd = c(1, 51, 54), rightEnd = c(50, 53, 100), x0 = 0)
  compare <- statAll(y = testy, filter = testfilter, fit = comparefit, singleStat = singleStatLR,
                     lengths = 1:20, add = 39:53, sd = sdrobnorm(testy, lag = testfilter$len + 1L), regu = 1)
  expect_equal(ret$stat, compare, tolerance = 1e-12)


  testfit <-  stepblock(c(0, 1, 0, 1), leftEnd = c(0, 25, 50, 75) / testfilter$sr,
                        rightEnd = c(25, 50, 75, testn) / testfilter$sr, x0 = 0)
  expect_warning(ret <- computeStat(y = testy, filter = testfilter, fit = testfit, family = "LR"))
  expect_identical(ret$stat, rep(-Inf, 20))


  testfit <-  stepblock(c(0, 1), leftEnd = c(0, 20 / testfilter$sr),
                        rightEnd = c(20 / testfilter$sr, testn / testfilter$sr), x0 = 0)
  expect_warning(ret <- computeStat(y = testy, filter = testfilter, fit = testfit, family = "LR"))
  comparefit <- stepblock(c(NA, stats::median(testy[31:89])),
                          leftEnd = c(1, 21), rightEnd = c(20, 100), x0 = 0)
  compare <- statAll(y = testy, filter = testfilter, fit = comparefit, singleStat = singleStatLR,
                     lengths = 1:20, add = 1:20, sd = sdrobnorm(testy, lag = testfilter$len + 1L), regu = 1)
  expect_equal(ret$stat, compare, tolerance = 1e-12)


  testfit <-  stepblock(c(0, 1), leftEnd = c(0, 80 / testfilter$sr),
                        rightEnd = c(80 / testfilter$sr, testn / testfilter$sr), x0 = 0)
  expect_warning(ret <- computeStat(y = testy, filter = testfilter, fit = testfit, family = "LR"))
  comparefit <- stepblock(c(stats::median(testy[11:69]), NA),
                          leftEnd = c(1, 81), rightEnd = c(80, 100), x0 = 0)
  compare <- statAll(y = testy, filter = testfilter, fit = comparefit, singleStat = singleStatLR,
                     lengths = 1:20, add = 69:100, sd = sdrobnorm(testy, lag = testfilter$len + 1L), regu = 1)
  expect_equal(ret$stat, compare, tolerance = 1e-12)

  testfit <-  stepR::stepblock(0, leftEnd = 1, rightEnd = 1 + testn / testfilter$sr, x0 = 1)
  expect_error(computeStat(testy, family = "LR", filter = testfilter, fit = testfit))
})

test_that("arguments in ... are tested and work for family LR", {
  testn <- 50L
  testfilter <- lowpassFilter::lowpassFilter(type = "bessel", param = list(pole = 4, cutoff = 0.1), sr = 1e4)
  testy <- lowpassFilter::randomGenerationMA(n = testn, filter = testfilter, seed = "no")
  testfit <-  stepblock(0, leftEnd = 0, rightEnd = testn / testfilter$sr, x0 = 0)

  expect_error(computeStat(testy, family = "LR", filter = unclass(testfilter), fit = testfit))

  expect_error(computeStat(testy, family = "LR", filter = testfilter, fit = unclass(testfit)))
  expect_identical(computeStat(testy, family = "LR", filter = testfilter, fit = list(fit = testfit)),
                   computeStat(testy, family = "LR", filter = testfilter, fit = testfit))

  expect_identical(computeStat(testy, family = "LR", filter = testfilter, fit = testfit),
                   computeStat(testy, family = "LR", filter = testfilter, fit = testfit,
                               startTime = 0, thresholdLongSegment = 10L, localValue = stats::median,
                               regularization = 1, suppressWarningNoDeconvolution = FALSE, localList = NULL))

  expect_error(computeStat(testy, family = "LR", filter = testfilter, fit = testfit, startTime = NA))
  expect_error(computeStat(testy, family = "LR", filter = testfilter, fit = testfit, startTime = c(0, 1)))
  expect_error(computeStat(testy, family = "LR", filter = testfilter, fit = testfit, startTime = 1))
  expect_error(computeStat(testy, family = "LR", filter = testfilter, fit = testfit, startTime = -1))

  testfitStartTime <-  stepblock(0, leftEnd = -1, rightEnd = testn / testfilter$sr - 1, x0 = -1)
  expect_identical(computeStat(testy, family = "LR", filter = testfilter, fit = testfitStartTime, startTime = -1),
                   computeStat(testy, family = "LR", filter = testfilter, fit = testfit, startTime = 0))

  testfitStartTime <-  stepblock(0, leftEnd = 0.001, rightEnd = testn / testfilter$sr + 0.001, x0 = 0.001)
  expect_identical(computeStat(testy, family = "LR", filter = testfilter, fit = testfitStartTime, startTime = 0.001),
                   computeStat(testy, family = "LR", filter = testfilter, fit = testfit, startTime = 0))

  expect_error(computeStat(testy, family = "LR", filter = testfilter, fit = testfit, thresholdLongSegment = NA))
  expect_error(computeStat(testy, family = "LR", filter = testfilter, fit = testfit, thresholdLongSegment = -1L))
  expect_error(computeStat(testy, family = "LR", filter = testfilter, fit = testfit,
                           thresholdLongSegment = c(10L, 23L)))

  expect_identical(computeStat(testy, family = "LR", filter = testfilter, fit = testfit, thresholdLongSegment = 10.5),
                   computeStat(testy, family = "LR", filter = testfilter, fit = testfit))
  expect_warning(computeStat(testy, family = "LR", filter = testfilter, fit = testfit, thresholdLongSegment = 30))


  expect_error(computeStat(testy, family = "LR", filter = testfilter, fit = testfit, localValue = 1))
  expect_error(computeStat(testy, family = "LR", filter = testfilter, fit = testfit, localValue = function() 1))
  expect_error(computeStat(testy, family = "LR", filter = testfilter, fit = testfit, localValue = function(a) NA))
  expect_error(computeStat(testy, family = "LR", filter = testfilter, fit = testfit, localValue = function(a) c(1, 2)))

  ret <- computeStat(testy, family = "LR", filter = testfilter, fit = testfit, localValue = mean)
  testcor <- testfilter$acf
  testcor[1] <- 2
  comparefit <- stepblock(mean(testy[11:39]), leftEnd = c(1), rightEnd = c(testn), x0 = 0)
  compare <- statAll(y = testy, filter = testfilter, fit = comparefit, singleStat = singleStatLR,
                     lengths = 1:20, add = integer(0), sd = sdrobnorm(testy, lag = testfilter$len + 1L), regu = 1)
  expect_equal(ret$stat, compare, tolerance = 1e-12)


  expect_error(computeStat(testy, family = "LR", filter = testfilter, fit = testfit, regularization = NULL))
  expect_error(computeStat(testy, family = "LR", filter = testfilter, fit = testfit, regularization = c(1, NA)))

  ret <- computeStat(testy, family = "LR", filter = testfilter, fit = testfit, regularization = 2)
  testcor <- testfilter$acf
  testcor[1] <- 3
  comparefit <- stepblock(stats::median(testy[11:39]), leftEnd = c(1), rightEnd = c(testn), x0 = 0)
  compare <- statAll(y = testy, filter = testfilter, fit = comparefit, singleStat = singleStatLR,
                     lengths = 1:20, add = integer(0), sd = sdrobnorm(testy, lag = testfilter$len + 1L), regu = 1)
  expect_equal(ret$stat, compare, tolerance = 1e-12)

  ret <- computeStat(testy, family = "LR", filter = testfilter, fit = testfit, regularization = c(2, 1))
  testcor <- testfilter$acf
  testcor[1] <- 3
  testcor[2] <- testcor[2] + 1
  comparefit <- stepblock(stats::median(testy[11:39]), leftEnd = c(1), rightEnd = c(testn), x0 = 0)
  compare <- statAll(y = testy, filter = testfilter, fit = comparefit, singleStat = singleStatLR,
                     lengths = 1:20, add = integer(0), sd = sdrobnorm(testy, lag = testfilter$len + 1L), regu = 1)
  expect_equal(ret$stat, compare, tolerance = 1e-12)


  expect_error(computeStat(testy, family = "LR", filter = testfilter, fit = testfit, correlations = "s"))
  expect_error(computeStat(testy, family = "LR", filter = testfilter, fit = testfit, correlations = c(Inf, 3)))
  expect_error(computeStat(testy, family = "LR", filter = testfilter, fit = testfit, correlations = c(1, 0)))
  expect_identical(computeStat(testy, family = "LR", filter = testfilter, fit = testfit,
                               correlations = testcor), ret)

  expect_error(computeStat(testy, family = "LR", filter = testfilter, fit = testfit,
                           suppressWarningNoDeconvolution = 1))
  expect_error(computeStat(testy, family = "LR", filter = testfilter, fit = testfit,
                           suppressWarningNoDeconvolution = as.logical(NA)))
  expect_error(computeStat(testy, family = "LR", filter = testfilter, fit = testfit,
                           suppressWarningNoDeconvolution = c(TRUE, TRUE)))

  testfit <- stepblock(c(0, 1, 2, 3), leftEnd = c(0, 40, 41, 42) / testfilter$sr,
                       rightEnd = c(40, 41, 42, testn) / testfilter$sr, x0 = 0)
  expect_warning(ret <- computeStat(testy, family = "LR", filter = testfilter, fit = testfit))
  expect_identical(computeStat(testy, family = "LR", filter = testfilter, fit = testfit,
                               suppressWarningNoDeconvolution = TRUE), ret)
})

test_that("family 2Param works (short version)", {
  testn <- 70L
  testfilter <- lowpassFilter::lowpassFilter(type = "bessel", param = list(pole = 4, cutoff = 0.1), sr = 1e4)
  testy <- lowpassFilter::randomGenerationMA(n = testn, filter = testfilter, seed = "no")
  testfit <-  stepblock(0, leftEnd = 0, rightEnd = testn / testfilter$sr, x0 = 0)

  expect_error(computeStat(family = "2Param", filter = testfilter, fit = testfit))
  expect_error(computeStat(y = numeric(0), family = "2Param", filter = testfilter, fit = testfit))
  expect_error(computeStat(y = as.integer(testy), family = "2Param", filter = testfilter, fit = testfit))
  expect_error(computeStat(y = c(testy[-1], "s"), family = "2Param", filter = testfilter, fit = testfit))
  expect_error(computeStat(y = c(testy[-1], NA), family = "2Param", filter = testfilter, fit = testfit))

  # this takes very long
  # expect_identical(computeStat(testy, family = "2Param", filter = testfilter, fit = testfit),
  #                  computeStat(testy, family = "2Param", intervalSystem = "all", lengths = 1:65,
  #                              penalty = "none", nq = testn, output = "list", filter = testfilter, fit = testfit))

  testlengths = c(1:3, 15, 64:65)
  ret <- computeStat(y = testy, filter = testfilter, fit = testfit, family = "2Param", lengths = testlengths)
  comparefit <- stepblock(stats::median(testy[11:59]), leftEnd = c(1), rightEnd = c(testn), x0 = 0)
  comparefit$var <- sdrobnorm(testy[11:59], lag = testfilter$len + 1L)^2
  compare <- statAll(y = testy, filter = testfilter, fit = comparefit, singleStat = singleStat2Param,
                     lengths = testlengths, add = integer(0), regu = 1)
  expect_equal(ret$stat, compare, tolerance = 1e-12)
})

test_that("family 2Param works (long version)", {
  testthat::skip_on_cran()
  testn <- 70L
  testfilter <- lowpassFilter::lowpassFilter(type = "bessel", param = list(pole = 4, cutoff = 0.1), sr = 1e4)
  testy <- lowpassFilter::randomGenerationMA(n = testn, filter = testfilter, seed = "no")
  testfit <-  stepblock(0, leftEnd = 0, rightEnd = testn / testfilter$sr, x0 = 0)

  expect_identical(computeStat(testy, family = "2Param", filter = testfilter, lengths = c(1, 2, 4, 8, 16, 32, 64)),
                   computeStat(testy, family = "2Param", intervalSystem = "dyaLen", filter = testfilter))

  ret <- computeStat(y = testy, filter = testfilter, fit = testfit, family = "2Param", lengths = c(1:3, 5, 25))
  expect_identical(ret, computeStat(y = testy, filter = testfilter, family = "2Param", lengths = c(1:3, 5, 25)))
  testcor <- testfilter$acf
  testcor[1] <- 2
  comparefit <- stepblock(stats::median(testy[11:59]), leftEnd = c(1), rightEnd = c(testn), x0 = 0)
  comparefit$var <- sdrobnorm(testy[11:59], lag = testfilter$len + 1L)^2
  compare <- statAll(y = testy, filter = testfilter, fit = comparefit, singleStat = singleStat2Param,
                     lengths = c(1:3, 5, 25), add = integer(0), regu = 1)
  expect_equal(ret$stat, compare, tolerance = 1e-12)

  expect_warning(ret2 <- computeStat(y = testy, filter = testfilter, signal = testfit,
                                     fit = testfit, family = "2Param", lengths = c(1:3, 5, 25)))
  expect_identical(ret, ret2)

  expect_error(computeStat(y = testy, filter = testfilter, intervalSystem = "dyaPar",
                           fit = testfit, family = "2Param", lengths = c(1, 4, 16)))

  expect_identical(computeStat(y = testy, filter = testfilter, family = "2Param", lengths = c(1:3, 5, 25),
                               output = "maximum"), ret$maximum)
  expect_identical(computeStat(y = testy, filter = testfilter, family = "2Param", lengths = c(1:3, 5, 25),
                               output = "vector"), ret$stat)

  ret <- computeStat(y = testy, filter = testfilter, fit = testfit, family = "2Param", lengths = c(1:3, 5, 25),
                     penalty = "sqrt", nq = 100L, output = "maximum")
  expect_equal(ret, max(sqrt(2 * compare) - sqrt(2 * log(exp(1) * 100L / c(1:3, 5, 25)))))

  suppressWarnings(ret <- computeStat(y = rnorm(20), filter = testfilter, family = "2Param"))
  expect_identical(ret$maximum, -Inf)
  expect_identical(ret$stat, rep(-Inf, 20))
  expect_identical(ret$lengths, 1:20)

  testn <- 100L
  testfit <-  stepblock(c(0, 1), leftEnd = c(0, 50 / testfilter$sr),
                        rightEnd = c(50 / testfilter$sr, testn / testfilter$sr), x0 = 0)
  testy <- lowpassFilter::randomGenerationMA(n = testn, filter = testfilter, seed = "no")

  ret <- computeStat(y = testy, filter = testfilter, fit = testfit, family = "2Param", lengths = c(1:3, 5, 25))
  testcor <- testfilter$acf
  testcor[1] <- 2
  j <- as.integer(testJump(39:50 / testfilter$sr, 40:60 / testfilter$sr, testy[40:60], testfilter, testcor,
                           stats::median(testy[11:39]), stats::median(testy[61:89])) * testfilter$sr + 1e-6)
  comparefit <- stepblock(c(stats::median(testy[11:39]), stats::median(testy[61:89])),
                          leftEnd = c(1, j + 1), rightEnd = c(j, 100), x0 = 0)
  comparefit$var <- c(sdrobnorm(testy[11:39], lag = testfilter$len + 1L)^2,
                      sdrobnorm(testy[61:89], lag = testfilter$len + 1L)^2)
  compare <- statAll(y = testy, filter = testfilter, fit = comparefit, singleStat = singleStat2Param,
                     lengths = c(1:3, 5, 25), add = integer(0), regu = 1)
  expect_equal(ret$stat, compare, tolerance = 1e-4)


  testfit <-  stepblock(c(0, 1, 0), leftEnd = c(0, 50 / testfilter$sr, 53 / testfilter$sr),
                        rightEnd = c(50 / testfilter$sr, 53 / testfilter$sr, testn / testfilter$sr), x0 = 0)
  ret <- computeStat(y = testy, filter = testfilter, fit = testfit, family = "2Param", lengths = c(1:2, 5, 13))
  comparefit <- stepblock(c(stats::median(testy[11:39]), NA, stats::median(testy[64:89])),
                          leftEnd = c(1, 51, 54), rightEnd = c(50, 53, 100), x0 = 0)
  comparefit$var <- c(sdrobnorm(testy[11:39], lag = testfilter$len + 1L)^2, NA,
                      sdrobnorm(testy[64:89], lag = testfilter$len + 1L)^2)
  compare <- statAll(y = testy, filter = testfilter, fit = comparefit, singleStat = singleStat2Param,
                     lengths = c(1:2, 5, 13), add = 39:53, regu = 1)
  expect_equal(ret$stat, compare, tolerance = 1e-12)


  testfit <-  stepblock(c(0, 1, 0, 1, 0), leftEnd = c(0, 50:53 / testfilter$sr),
                        rightEnd = c(50:53 / testfilter$sr, testn / testfilter$sr), x0 = 0)
  expect_warning(ret <- computeStat(y = testy, filter = testfilter, fit = testfit,
                                    family = "2Param", lengths = c(1:3, 5, 15)))
  comparefit <- stepblock(c(stats::median(testy[11:39]), NA, stats::median(testy[64:89])),
                          leftEnd = c(1, 51, 54), rightEnd = c(50, 53, 100), x0 = 0)
  comparefit$var <- c(sdrobnorm(testy[11:39], lag = testfilter$len + 1L)^2, NA,
                      sdrobnorm(testy[64:89], lag = testfilter$len + 1L)^2)
  compare <- statAll(y = testy, filter = testfilter, fit = comparefit, singleStat = singleStat2Param,
                     lengths = c(1:3, 5, 15), add = 39:53, regu = 1)
  expect_equal(ret$stat, compare, tolerance = 1e-12)


  testfit <-  stepblock(c(0, 1, 0, 1), leftEnd = c(0, 25, 50, 75) / testfilter$sr,
                        rightEnd = c(25, 50, 75, testn) / testfilter$sr, x0 = 0)
  expect_warning(ret <- computeStat(y = testy, filter = testfilter, fit = testfit, family = "2Param",
                                    lengths = c(1:3, 5, 25)))
  expect_identical(ret$stat, rep(-Inf, 5))


  testfit <-  stepblock(c(0, 1), leftEnd = c(0, 20 / testfilter$sr),
                        rightEnd = c(20 / testfilter$sr, testn / testfilter$sr), x0 = 0)
  expect_warning(ret <- computeStat(y = testy, filter = testfilter, fit = testfit, family = "2Param",
                                    lengths = c(7:9, 12)))
  comparefit <- stepblock(c(NA, stats::median(testy[31:89])),
                          leftEnd = c(1, 21), rightEnd = c(20, 100), x0 = 0)
  comparefit$var <- c(NA, sdrobnorm(testy[31:89], lag = testfilter$len + 1L)^2)
  compare <- statAll(y = testy, filter = testfilter, fit = comparefit, singleStat = singleStat2Param,
                     lengths = c(7:9, 12), add = 1:20, regu = 1)
  expect_equal(ret$stat, compare, tolerance = 1e-12)


  testfit <-  stepblock(c(0, 1), leftEnd = c(0, 80 / testfilter$sr),
                        rightEnd = c(80 / testfilter$sr, testn / testfilter$sr), x0 = 0)
  expect_warning(ret <- computeStat(y = testy, filter = testfilter, fit = testfit, family = "2Param",
                                    lengths = c(11:13)))
  comparefit <- stepblock(c(stats::median(testy[11:69]), NA),
                          leftEnd = c(1, 81), rightEnd = c(80, 100), x0 = 0)
  comparefit$var <- c(sdrobnorm(testy[11:69], lag = testfilter$len + 1L)^2, NA)
  compare <- statAll(y = testy, filter = testfilter, fit = comparefit, singleStat = singleStat2Param,
                     lengths = c(11:13), add = 69:100, regu = 1)
  expect_equal(ret$stat, compare, tolerance = 1e-12)

  testfit <-  stepR::stepblock(0, leftEnd = 1, rightEnd = 1 + testn / testfilter$sr, x0 = 1)
  expect_error(computeStat(testy, family = "2Param", filter = testfilter, fit = testfit, lengths = testlengths))
})

test_that("arguments in ... are tested and work for family 2Param", {
  testthat::skip_on_cran()
  testn <- 50L
  testfilter <- lowpassFilter::lowpassFilter(type = "bessel", param = list(pole = 4, cutoff = 0.1), sr = 1e4)
  testy <- lowpassFilter::randomGenerationMA(n = testn, filter = testfilter, seed = "no")
  testfit <-  stepblock(0, leftEnd = 0, rightEnd = testn / testfilter$sr, x0 = 0)

  expect_error(computeStat(testy, family = "2Param", filter = unclass(testfilter), fit = testfit))

  expect_error(computeStat(testy, family = "2Param", filter = testfilter, fit = unclass(testfit)))
  expect_identical(computeStat(testy, family = "2Param", filter = testfilter, fit = list(fit = testfit),
                               lengths = c(7:9, 12)),
                   computeStat(testy, family = "2Param", filter = testfilter, fit = testfit, lengths = c(7:9, 12)))

  expect_identical(computeStat(testy, family = "2Param", filter = testfilter, fit = testfit, lengths = c(7:9, 12)),
                   computeStat(testy, family = "2Param", filter = testfilter, fit = testfit, lengths = c(7:9, 12),
                               startTime = 0, thresholdLongSegment = 10L, localValue = stats::median,
                               localVar = function(data) sdrobnorm(data, lag = testfilter$len + 1L)^2,
                               regularization = 1, suppressWarningNoDeconvolution = FALSE, localList = NULL))

  expect_error(computeStat(testy, family = "2Param", filter = testfilter, fit = testfit, startTime = NA))
  expect_error(computeStat(testy, family = "2Param", filter = testfilter, fit = testfit, startTime = c(0, 1)))
  expect_error(computeStat(testy, family = "2Param", filter = testfilter, fit = testfit,
                           startTime = (testn + 1) / testfilter$sr))
  expect_error(computeStat(testy, family = "2Param", filter = testfilter, fit = testfit, startTime = -1 / testfilter$sr))

  testfitStartTime <-  stepblock(0, leftEnd = -1, rightEnd = testn / testfilter$sr - 1, x0 = -1)
  expect_identical(computeStat(testy, family = "2Param", filter = testfilter, fit = testfitStartTime, startTime = -1),
                   computeStat(testy, family = "2Param", filter = testfilter, fit = testfit, startTime = 0))

  testfitStartTime <-  stepblock(0, leftEnd = 0.001, rightEnd = testn / testfilter$sr + 0.001, x0 = 0.001)
  expect_identical(computeStat(testy, family = "2Param", filter = testfilter, fit = testfitStartTime, startTime = 0.001),
                   computeStat(testy, family = "2Param", filter = testfilter, fit = testfit, startTime = 0))

  expect_error(computeStat(testy, family = "2Param", filter = testfilter, fit = testfit, thresholdLongSegment = NA))
  expect_error(computeStat(testy, family = "2Param", filter = testfilter, fit = testfit, thresholdLongSegment = -1L))
  expect_error(computeStat(testy, family = "2Param", filter = testfilter, fit = testfit,
                           thresholdLongSegment = c(10L, 23L)))

  expect_identical(computeStat(testy, family = "2Param", filter = testfilter, fit = testfit,
                               thresholdLongSegment = 25.5, lengths = c(7:9, 12, 37)),
                   computeStat(testy, family = "2Param", filter = testfilter, fit = testfit, lengths = c(7:9, 12, 37)))
  expect_warning(computeStat(testy, family = "2Param", filter = testfilter, fit = testfit, thresholdLongSegment = 30,
                             lengths = c(7:9, 12)))

  expect_error(computeStat(testy, family = "2Param", filter = testfilter, fit = testfit, localValue = 1))
  expect_error(computeStat(testy, family = "2Param", filter = testfilter, fit = testfit, localValue = function() 1))
  expect_error(computeStat(testy, family = "2Param", filter = testfilter, fit = testfit, localValue = function(a) NA))
  expect_error(computeStat(testy, family = "2Param", filter = testfilter, fit = testfit,
                           localValue = function(a) c(1, 2)))

  expect_error(computeStat(testy, family = "2Param", filter = testfilter, fit = testfit, localVar = 1))
  expect_error(computeStat(testy, family = "2Param", filter = testfilter, fit = testfit, localVar = function() 1))
  expect_error(computeStat(testy, family = "2Param", filter = testfilter, fit = testfit, localVar = function(a) NA))
  expect_error(computeStat(testy, family = "2Param", filter = testfilter, fit = testfit,
                           localVar = function(a) c(1, 2)))
  expect_error(computeStat(testy, family = "2Param", filter = testfilter, fit = testfit, localVar = function(a) -1))

  ret <- computeStat(testy, family = "2Param", filter = testfilter, fit = testfit, localValue = mean,
                     localVar = function(data) 1, lengths = c(7:9, 12))
  testcor <- testfilter$acf
  testcor[1] <- 2
  comparefit <- stepblock(mean(testy[11:39]), leftEnd = c(1), rightEnd = c(testn), x0 = 0)
  comparefit$var <- 1
  compare <- statAll(y = testy, filter = testfilter, fit = comparefit, singleStat = singleStat2Param,
                     lengths = c(7:9, 12), add = integer(0), regu = 1)
  expect_equal(ret$stat, compare, tolerance = 1e-12)

  expect_error(computeStat(testy, family = "2Param", filter = testfilter, fit = testfit, regularization = NULL))
  expect_error(computeStat(testy, family = "2Param", filter = testfilter, fit = testfit, regularization = c(1, NA)))

  ret <- computeStat(testy, family = "2Param", filter = testfilter, fit = testfit, regularization = 2,
                     lengths = c(7:9, 12))
  testcor <- testfilter$acf
  testcor[1] <- 3
  comparefit <- stepblock(stats::median(testy[11:39]), leftEnd = c(1), rightEnd = c(testn), x0 = 0)
  comparefit$var <- sdrobnorm(testy[11:39], lag = 12)^2
  compare <- statAll(y = testy, filter = testfilter, fit = comparefit, singleStat = singleStat2Param,
                     lengths = c(7:9, 12), add = integer(0), regu = 1)
  expect_equal(ret$stat, compare, tolerance = 1e-12)


  expect_error(computeStat(testy, family = "2Param", filter = testfilter, fit = testfit, correlations = "s"))
  expect_error(computeStat(testy, family = "2Param", filter = testfilter, fit = testfit, correlations = c(Inf, 3)))
  expect_error(computeStat(testy, family = "2Param", filter = testfilter, fit = testfit, correlations = c(1, 0)))
  expect_identical(computeStat(testy, family = "2Param", filter = testfilter, fit = testfit,
                               correlations = testcor, lengths = c(7:9, 12)), ret)

  expect_error(computeStat(testy, family = "2Param", filter = testfilter, fit = testfit,
                           suppressWarningNoDeconvolution = 1))
  expect_error(computeStat(testy, family = "2Param", filter = testfilter, fit = testfit,
                           suppressWarningNoDeconvolution = as.logical(NA)))
  expect_error(computeStat(testy, family = "2Param", filter = testfilter, fit = testfit,
                           suppressWarningNoDeconvolution = c(TRUE, TRUE)))

  testfit <- stepblock(c(0, 1, 2, 3), leftEnd = c(0, 47, 48, 49) / testfilter$sr,
                       rightEnd = c(47, 48, 49, testn) / testfilter$sr, x0 = 0)
  expect_warning(ret <- computeStat(testy, family = "2Param", filter = testfilter, fit = testfit, lengths = c(7:9, 12)))
  expect_identical(computeStat(testy, family = "2Param", filter = testfilter, fit = testfit, lengths = c(7:9, 12),
                               suppressWarningNoDeconvolution = TRUE), ret)
})
