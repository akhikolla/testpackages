context("test-methods-disc")

test_that("show works", {
  x <- matrix(runif(120), nrow = 40, ncol = 3)
  cuts <- seq(0, 1, length.out = 4)
  xd <- apply(x, 2, function(col) as.numeric(cut(col, cuts)))
  theta <- t(matrix(c(0, 0, 0, 2, 2, 2, -2, -2, -2), ncol = 3, nrow = 3))
  log_odd <- rowSums(t(sapply(seq_along(xd[, 1]), function(row_id) sapply(seq_along(xd[row_id, ]), function(element) theta[xd[row_id, element], element]))))
  y <- rbinom(40, 1, 1 / (1 + exp(-log_odd)))
  sem_disc <- glmdisc(x, y, iter = 15, m_start = 4, test = FALSE, validation = FALSE, criterion = "aic")
  expect_output(show(sem_disc))
})

test_that("print works", {
  x <- matrix(runif(120), nrow = 40, ncol = 3)
  cuts <- seq(0, 1, length.out = 4)
  xd <- apply(x, 2, function(col) as.numeric(cut(col, cuts)))
  theta <- t(matrix(c(0, 0, 0, 2, 2, 2, -2, -2, -2), ncol = 3, nrow = 3))
  log_odd <- rowSums(t(sapply(seq_along(xd[, 1]), function(row_id) sapply(seq_along(xd[row_id, ]), function(element) theta[xd[row_id, element], element]))))
  y <- rbinom(40, 1, 1 / (1 + exp(-log_odd)))
  sem_disc <- glmdisc(x, y, iter = 15, m_start = 4, test = FALSE, validation = FALSE, criterion = "aic")
  expect_output(print(sem_disc))
})

test_that("summary works", {
  x <- matrix(runif(120), nrow = 40, ncol = 3)
  cuts <- seq(0, 1, length.out = 4)
  xd <- apply(x, 2, function(col) as.numeric(cut(col, cuts)))
  theta <- t(matrix(c(0, 0, 0, 2, 2, 2, -2, -2, -2), ncol = 3, nrow = 3))
  log_odd <- rowSums(t(sapply(seq_along(xd[, 1]), function(row_id) sapply(seq_along(xd[row_id, ]), function(element) theta[xd[row_id, element], element]))))
  y <- rbinom(40, 1, 1 / (1 + exp(-log_odd)))
  sem_disc <- glmdisc(x, y, iter = 15, m_start = 4, test = FALSE, validation = FALSE, criterion = "aic")
  expect_type(summary(sem_disc), "character")
})
