library(testthat)
library(gif)
library(MASS)
context("hgt function")
skip_on_cran()

test_that("Error if computation matrix for hgt is wrong!", {
  target_matrix <- matrix(c(1.0725, 0.4958, 0, 0.4958, 1.0196, 0.4319, 0, 0.4319, 0.9802), nrow = 3, ncol = 3)
  set.seed(1)
  num <- 100
  p <- 3
  sigma_inv <- diag(1, p, p)
  for(i in 1:(p - 1)) {
    sigma_inv[i, i + 1] <- 0.5
    sigma_inv[i + 1, i] <- 0.5
  }
  x <- ggm.generator(num, sigma_inv)
  non_zero_num <- sum(sigma_inv != 0) - p
  res <- hgt(x, size = non_zero_num / 2)
  omega <- as.matrix(res[["Omega"]])
  omega <- round(omega, 4)
  for(i in 1:3) {
    for(j in 1:3) {
      expect_equal(omega[i, j], target_matrix[i, j])
    }
  }
})

test_that("Error if the size of the output active.entry doesn't match the arguments size in hgt!", {
  set.seed(1)
  num <- 100
  p <- 3
  sigma_inv <- diag(1, p, p)
  for(i in 1:(p - 1)) {
    sigma_inv[i, i + 1] <- 0.5
    sigma_inv[i + 1, i] <- 0.5
  }
  non_zero_num <- sum(sigma_inv != 0) - p
  x <- ggm.generator(num, sigma_inv)
  res <- hgt(x, size = non_zero_num / 2)
  active.entry <- res[["active.entry"]]
  expect_equal(nrow(active.entry), non_zero_num / 2)
})

test_that("Error if active.entry is changed when it is given directly in hgt!", {
  set.seed(1)
  num <- 100
  p <- 3
  sigma_inv <- diag(1, p, p)
  for(i in 1:(p - 1)) {
    sigma_inv[i, i + 1] <- 0.5
    sigma_inv[i + 1, i] <- 0.5
  }
  x <- ggm.generator(num, sigma_inv)
  non_zero_index <- which(as.matrix(sigma_inv) != 0, arr.ind = TRUE)
  active.entry <- non_zero_index[which(non_zero_index[,1] < non_zero_index[,2]),]
  res <- hgt(x, active.entry = active.entry)
  expect_equal(res[["active.entry"]], active.entry)
})

test_that("Error if computation matrix for sgt is wrong!", {
  target_matrix <- matrix(c(1.1362, 0.4846, -0.2459, 0.4846, 1.0048, 0.4224, -0.2459, 0.4224, 1.0430), nrow = 3, ncol = 3)
  set.seed(1)
  num <- 100
  p <- 3
  sigma_inv <- diag(1, p, p)
  for(i in 1:(p - 1)) {
    sigma_inv[i, i + 1] <- 0.5
    sigma_inv[i + 1, i] <- 0.5
  }
  x <- ggm.generator(num, sigma_inv)
  res <- sgt(x, lambda = 0.01)
  omega <- as.matrix(res[["Omega"]])
  omega <- round(omega, 4)
  for(i in 1:3) {
    for(j in 1:3) {
      expect_equal(omega[i, j], target_matrix[i, j])
    }
  }
})

test_that("Error if the size of the output active.entry doesn't match the arguments size in sgt!", {
  set.seed(1)
  num <- 100
  p <- 3
  sigma_inv <- diag(1, p, p)
  for(i in 1:(p - 1)) {
    sigma_inv[i, i + 1] <- 0.5
    sigma_inv[i + 1, i] <- 0.5
  }
  non_zero_num <- sum(sigma_inv != 0) - p
  x <- ggm.generator(num, sigma_inv)
  res <- sgt(x, size = non_zero_num / 2)
  active.entry <- res[["active.entry"]]
  expect_equal(nrow(active.entry), non_zero_num / 2)
})

test_that("Error if the acyclic diagnostic is wrong in sgt!", {
  set.seed(1)
  num <- 100
  p <- 3
  sigma_inv <- diag(1, p, p)
  for(i in 1:(p - 1)) {
    sigma_inv[i, i + 1] <- 0.5
    sigma_inv[i + 1, i] <- 0.5
  }
  non_zero_num <- sum(sigma_inv != 0) - p
  x <- ggm.generator(num, sigma_inv)
  res <- sgt(x, size = non_zero_num / 2)
  expect_equal(res[["is.acyclic"]], TRUE)
})
