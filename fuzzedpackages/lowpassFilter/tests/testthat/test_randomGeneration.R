
compareRandomGeneration <- function(signal, noise, filter, oversampling) {
  n <- length(signal)
  
  diffFilter <- diff(filter$truncatedStepfun(seq(0, filter$len, 1 / oversampling) / filter$sr))
  eps <- stats::rnorm((n + filter$len - 1) * oversampling, 0, noise / sqrt(sum(diffFilter^2)))
  
  y <- numeric(n)
  for (i in 1:n) {
    y[i] <- sum(rev(eps[(i - 1) * oversampling + seq(1, filter$len * oversampling, 1)]) * diffFilter)
  }
  signal + y
}

compareRandomGenerationMA <- function(signal, noise, filter) {
  testComputeMA <- function(cov) {
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
  
  n <- length(signal)
  kern <- c(1, testComputeMA(filter$acf))
  z <- stats::rnorm(n + length(kern) - 1, sd = 1)
  signal + .convolve(z, kern) / sqrt(sum(kern^2)) * noise
}

test_that("everything works for default arguments", {
  testn <- 100L
  testfilter <- lowpassFilter(type = "bessel", param = list(pole = 4, cutoff = 0.1))
  ret <- randomGeneration(n = testn, filter = testfilter)
  
  expect_error(randomGeneration(filter = testfilter))
  expect_error(randomGeneration(n = testn))
  
  expect_true(is.numeric(ret))
  expect_true(length(ret) == testn)
  expect_true(all(is.finite(ret)))
  
  expect_identical(ret, randomGeneration(n = testn, filter = testfilter, signal = 0, noise = 1,
                                         oversampling = 100L, seed = testn, startTime = 0, truncated = TRUE))
  
  set.seed(testn)
  compare <- compareRandomGeneration(signal = rep(0, 100L), noise = rep(1, 1001000),
                                     filter = testfilter, oversampling = 100L)
  expect_equal(ret, compare, tolerance = 1e-12)
  
  
  ret <- randomGenerationMA(n = testn, filter = testfilter)
  
  expect_error(randomGenerationMA(filter = testfilter))
  expect_error(randomGenerationMA(n = testn))
  
  expect_true(is.numeric(ret))
  expect_true(length(ret) == testn)
  expect_true(all(is.finite(ret)))
  
  expect_identical(ret, randomGenerationMA(n = testn, filter = testfilter, signal = 0, noise = 1,
                                           seed = testn, startTime = 0, truncated = TRUE))
  
  set.seed(testn)
  compare <- compareRandomGenerationMA(signal = rep(0, 100L), noise = 1, filter = testfilter)
  expect_equal(ret, compare, tolerance = 1e-12)
})

test_that("n is tested and works", {
  testfilter <- lowpassFilter(type = "bessel", param = list(pole = 4, cutoff = 0.1))

  expect_error(randomGeneration(n = "s", filter = testfilter))
  expect_error(randomGeneration(n = c(1L, 2L), filter = testfilter))
  expect_error(randomGeneration(n = as.integer(NA), filter = testfilter))
  expect_error(randomGeneration(n = NULL, filter = testfilter))
  expect_error(randomGeneration(n = Inf, filter = testfilter))
  expect_error(randomGeneration(n = 0L, filter = testfilter))
  expect_error(randomGeneration(n = -1.2, filter = testfilter))
  
  expect_identical(randomGeneration(n = 30, filter = testfilter),
                   randomGeneration(n = 30L, filter = testfilter))
  expect_identical(randomGeneration(n = 30.5, filter = testfilter),
                   randomGeneration(n = 30L, filter = testfilter))
  
  
  expect_error(randomGenerationMA(n = "s", filter = testfilter))
  expect_error(randomGenerationMA(n = c(1L, 2L), filter = testfilter))
  expect_error(randomGenerationMA(n = as.integer(NA), filter = testfilter))
  expect_error(randomGenerationMA(n = NULL, filter = testfilter))
  expect_error(randomGenerationMA(n = Inf, filter = testfilter))
  expect_error(randomGenerationMA(n = 0L, filter = testfilter))
  expect_error(randomGenerationMA(n = -1.2, filter = testfilter))
  
  expect_identical(randomGenerationMA(n = 30, filter = testfilter),
                   randomGenerationMA(n = 30L, filter = testfilter))
  expect_identical(randomGenerationMA(n = 30.5, filter = testfilter),
                   randomGenerationMA(n = 30L, filter = testfilter))
})

test_that("filter is tested and works", {
  testn <- 100L
  testfilter <- lowpassFilter(type = "bessel", param = list(pole = 4, cutoff = 0.1))
  ret <- randomGeneration(n = testn, filter = testfilter)
  
  expect_error(randomGeneration(n = testn, filter = unclass(testfilter)))
  expect_error(randomGenerationMA(n = testn, filter = unclass(testfilter)))

  expect_equal(ret, randomGeneration(n = testn, filter = 
                                       lowpassFilter(type = "bessel", param = list(pole = 4, cutoff = 0.1), sr = 1e4)))
  expect_equal(randomGenerationMA(n = testn, filter = testfilter),
               randomGenerationMA(n = testn, filter = 
                                    lowpassFilter(type = "bessel", param = list(pole = 4, cutoff = 0.1), sr = 1e4)))
  
  testfilter <- lowpassFilter(type = "bessel", param = list(pole = 4, cutoff = 0.1), len = 8L)
  ret <- randomGeneration(n = testn, filter = testfilter)
  expect_true(is.numeric(ret))
  expect_true(length(ret) == testn)
  expect_true(all(is.finite(ret)))
  set.seed(testn)
  compare <- compareRandomGeneration(signal = rep(0, 100L), noise = rep(1, 1000700),
                                     filter = testfilter, oversampling = 100L)
  expect_equal(ret, compare, tolerance = 1e-12)
  ret <- randomGenerationMA(n = testn, filter = testfilter)
  expect_true(is.numeric(ret))
  expect_true(length(ret) == testn)
  expect_true(all(is.finite(ret)))
  set.seed(testn)
  compare <- compareRandomGenerationMA(signal = rep(0, 100L), noise = 1, filter = testfilter)
  expect_equal(ret, compare, tolerance = 1e-12)
  
  testfilter <- lowpassFilter(param = list(pole = 4L, cutoff = 0.05), len = 14)
  ret <- randomGeneration(n = testn, filter = testfilter)
  expect_true(is.numeric(ret))
  expect_true(length(ret) == testn)
  expect_true(all(is.finite(ret)))
  set.seed(testn)
  compare <- compareRandomGeneration(signal = rep(0, 100L), noise = rep(1, 1000700),
                                     filter = testfilter, oversampling = 100L)
  expect_equal(ret, compare, tolerance = 1e-12)
  ret <- randomGenerationMA(n = testn, filter = testfilter)
  expect_true(is.numeric(ret))
  expect_true(length(ret) == testn)
  expect_true(all(is.finite(ret)))
  set.seed(testn)
  compare <- compareRandomGenerationMA(signal = rep(0, 100L), noise = 1, filter = testfilter)
  expect_equal(ret, compare, tolerance = 1e-12)
  
  testfilter <- lowpassFilter(param = list(pole = 6L, cutoff = 0.25), len = 5, shift = 0.1)
  ret <- randomGeneration(n = testn, filter = testfilter)
  expect_true(is.numeric(ret))
  expect_true(length(ret) == testn)
  expect_true(all(is.finite(ret)))
  set.seed(testn)
  compare <- compareRandomGeneration(signal = rep(0, 100L), noise = rep(1, 1000700),
                                     filter = testfilter, oversampling = 100L)
  expect_equal(ret, compare, tolerance = 1e-12)
  ret <- randomGenerationMA(n = testn, filter = testfilter)
  expect_true(is.numeric(ret))
  expect_true(length(ret) == testn)
  expect_true(all(is.finite(ret)))
  set.seed(testn)
  compare <- compareRandomGenerationMA(signal = rep(0, 100L), noise = 1, filter = testfilter)
  expect_equal(ret, compare, tolerance = 1e-12)
})

test_that("signal is tested and works", {
  testn <- 100L
  testfilter <- lowpassFilter(type = "bessel", param = list(pole = 4, cutoff = 0.1))
  testsignal <- data.frame(leftEnd = c(0, 25, 50, 75),
                      rightEnd = c(25, 50, 75, 100),
                      value = c(0, 1, 0, 1))
  attr(testsignal, "x0") <- 0
  class(testsignal) <- c("stepblock", class(testsignal))
  
  expect_error(randomGeneration(n = testn, filter = testfilter, signal = unclass(testsignal)))
  expect_error(randomGeneration(n = testn, filter = testfilter, signal = c(numeric(99), "s")))
  expect_error(randomGeneration(n = testn, filter = testfilter, signal = c(numeric(99), Inf)))
  expect_error(randomGeneration(n = testn, filter = testfilter, signal = numeric(90)))
  
  expect_error(randomGenerationMA(n = testn, filter = testfilter, signal = unclass(testsignal)))
  expect_error(randomGenerationMA(n = testn, filter = testfilter, signal = c(numeric(99), "s")))
  expect_error(randomGenerationMA(n = testn, filter = testfilter, signal = c(numeric(99), Inf)))
  expect_error(randomGenerationMA(n = testn, filter = testfilter, signal = numeric(90)))
  
  
  ret <- randomGeneration(n = testn, filter = testfilter, signal = testsignal)
  
  expect_true(is.numeric(ret))
  expect_true(length(ret) == testn)
  expect_true(all(is.finite(ret)))

  set.seed(testn)
  compare <- compareRandomGeneration(signal = getConvolution(t = 1:testn / testfilter$sr,
                                                             stepfun = testsignal, filter = testfilter,
                                                             truncated = TRUE),
                                     noise = rep(1, 1001000), filter = testfilter, oversampling = 100L)
  expect_equal(ret, compare, tolerance = 1e-12)
  
  ret <- randomGenerationMA(n = testn, filter = testfilter, signal = testsignal)
  expect_true(is.numeric(ret))
  expect_true(length(ret) == testn)
  expect_true(all(is.finite(ret)))

  set.seed(testn)
  compare <- compareRandomGenerationMA(signal = getConvolution(t = 1:testn / testfilter$sr,
                                                               stepfun = testsignal, filter = testfilter,
                                                               truncated = TRUE),
                                       noise = 1, filter = testfilter)
  expect_equal(ret, compare, tolerance = 1e-12)
  
  
  testsignal <- rep(1, testn)
  ret <- randomGeneration(n = testn, filter = testfilter, signal = testsignal)
  expect_identical(ret, randomGeneration(n = testn, filter = testfilter, signal = 1))
  
  expect_true(is.numeric(ret))
  expect_true(length(ret) == testn)
  expect_true(all(is.finite(ret)))
  
  set.seed(testn)
  compare <- compareRandomGeneration(signal = testsignal, noise = rep(1, 1001000),
                                     filter = testfilter, oversampling = 100L)
  expect_equal(ret, compare, tolerance = 1e-12)
  
  ret <- randomGenerationMA(n = testn, filter = testfilter, signal = testsignal)
  expect_identical(ret, randomGenerationMA(n = testn, filter = testfilter, signal = 1))
  expect_true(is.numeric(ret))
  expect_true(length(ret) == testn)
  expect_true(all(is.finite(ret)))
  
  set.seed(testn)
  compare <- compareRandomGenerationMA(signal = testsignal, noise = 1, filter = testfilter)
  expect_equal(ret, compare, tolerance = 1e-12)
  
  
  testsignal <- 1:testn / 10
  testfilter <- lowpassFilter(param = list(pole = 6L, cutoff = 0.25), len = 5, shift = 0.1)
  ret <- randomGeneration(n = testn, filter = testfilter, signal = testsignal)
  
  expect_true(is.numeric(ret))
  expect_true(length(ret) == testn)
  expect_true(all(is.finite(ret)))
  
  set.seed(testn)
  compare <- compareRandomGeneration(signal = testsignal, noise = rep(1, 1001000),
                                     filter = testfilter, oversampling = 100L)
  expect_equal(ret, compare, tolerance = 1e-12)
  
  ret <- randomGenerationMA(n = testn, filter = testfilter, signal = testsignal)
  expect_true(is.numeric(ret))
  expect_true(length(ret) == testn)
  expect_true(all(is.finite(ret)))
  
  set.seed(testn)
  compare <- compareRandomGenerationMA(signal = testsignal, noise = 1, filter = testfilter)
  expect_equal(ret, compare, tolerance = 1e-12)
})
