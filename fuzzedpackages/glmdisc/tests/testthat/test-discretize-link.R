context("test-discretize-link")

test_that("discretize_link works with one-dimensional continuous data", {
  x <- matrix(runif(40), nrow = 40, ncol = 1)
  cuts <- seq(0, 1, length.out = 4)
  xd <- as.numeric(cut(x, cuts))
  link <- list(nnet::multinom(e ~ x, data = data.frame(e = factor(xd), x = x), maxit = 20))
  res <- discretize_link(link, data.frame(x), 3)
  expect_type(res, "character")
  expect_length(res, 40)
  expect_equal(nlevels(factor(res)), 3)
})

test_that("discretize_link works with multi-dimensional continuous data", {
  x <- matrix(runif(120), nrow = 40, ncol = 3)
  cuts <- seq(0, 1, length.out = 4)
  xd <- apply(x, 2, function(col) as.numeric(cut(col, cuts)))
  link <- list()
  for (j in 1:3) {
    link[[j]] <- nnet::multinom(e ~ x, data = data.frame(e = factor(xd[, j]), x = x[, j]), maxit = 20)
  }
  res <- discretize_link(link, data.frame(x), 3)
  expect_type(res, "character")
  expect_length(res, 120)
  for (j in 1:3) {
    expect_equal(nlevels(factor(res[, j])), 3)
  }
})

test_that("discretize_link works with one-dimensional categorical data", {
  x <- matrix(runif(40), nrow = 40, ncol = 1)
  cuts <- seq(0, 1, length.out = 4)
  xd <- as.numeric(cut(x, cuts))
  e <- ifelse(xd %in% c(1, 2), 1, 2)
  link <- list(table(factor(e), factor(xd)))
  res <- discretize_link(link, data.frame(x = factor(xd)), 2)
  expect_type(res, "character")
  expect_length(res, 40)
  expect_equal(nlevels(factor(res)), 2)
})

test_that("discretize_link works with multi-dimensional categorical data", {
  x <- matrix(runif(120), nrow = 40, ncol = 3)
  cuts <- seq(0, 1, length.out = 4)
  xd <- apply(x, 2, function(col) as.numeric(cut(col, cuts)))
  e <- apply(xd, 2, function(col) ifelse(col %in% c(1, 2), 1, 2))
  link <- list()
  for (j in 1:3) {
    link[[j]] <- table(factor(e[, j]), factor(xd[, j]))
  }
  res <- discretize_link(link, data.frame(apply(xd, 2, factor), stringsAsFactors = TRUE), 2)
  expect_type(res, "character")
  expect_length(res, 120)
  for (j in 1:3) {
    expect_equal(nlevels(factor(res[, j])), 2)
  }
})

test_that("discretize_link works with multi-dimensional mixed data", {
  x <- matrix(runif(120), nrow = 40, ncol = 3)
  cuts <- seq(0, 1, length.out = 4)
  xd <- apply(x, 2, function(col) as.numeric(cut(col, cuts)))
  e <- apply(xd[, 1:2], 2, function(col) ifelse(col %in% c(1, 2), 1, 2))
  link <- list()
  for (j in 1:2) {
    link[[j]] <- table(factor(e[, j]), factor(xd[, j]))
  }
  link[[3]] <- nnet::multinom(e ~ x, data = data.frame(e = factor(xd[, 3]), x = x[, 3]), maxit = 50)
  res <- discretize_link(link, cbind(data.frame(apply(xd[, 1:2], 2, factor), stringsAsFactors = TRUE), x[, 3]), 2)
  expect_type(res, "character")
  expect_length(res, 120)
  for (j in 1:2) {
    expect_equal(nlevels(factor(res[, j])), 2)
  }
  expect_equal(nlevels(factor(res[, 3])), 3)
})

test_that("discretize_link errors", {
  x <- matrix(runif(120), nrow = 40, ncol = 3)
  cuts <- seq(0, 1, length.out = 4)
  xd <- apply(x, 2, function(col) as.numeric(cut(col, cuts)))
  e <- apply(xd[, 1:2], 2, function(col) ifelse(col %in% c(1, 2), 1, 2))
  link <- list()
  for (j in 1:2) {
    link[[j]] <- table(factor(e[, j]), factor(xd[, j]))
  }
  link[[3]] <- nnet::multinom(e ~ x, data = data.frame(e = factor(xd[, 3]), x = x[, 3]), maxit = 50)
  test <- data.frame(xd[, 1:2]) %>% dplyr::mutate_at(dplyr::vars("X1", "X2"), dplyr::funs(as.integer))
  test <- cbind(test, x[, 3])
  expect_error(discretize_link(link, test, 2))
})
