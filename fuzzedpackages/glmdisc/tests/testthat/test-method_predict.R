context("test-method_predict")

test_that("predict works with continuous data", {
  x <- matrix(runif(150), nrow = 50, ncol = 3)
  cuts <- seq(0, 1, length.out = 4)
  xd <- apply(x, 2, function(col) as.numeric(cut(col, cuts)))
  theta <- t(matrix(c(0, 0, 0, 2, 2, 2, -2, -2, -2), ncol = 3, nrow = 3))
  log_odd <- rowSums(t(sapply(seq_along(xd[, 1]), function(row_id) sapply(seq_along(xd[row_id, ]), function(element) theta[xd[row_id, element], element]))))
  y <- rbinom(50, 1, 1 / (1 + exp(-log_odd)))
  sem_disc <- glmdisc(x, y, iter = 20, m_start = 4, test = FALSE, validation = FALSE, criterion = "aic")
  pred_sem <- predict(sem_disc, x)

  expect_length(pred_sem, 50)
  expect_true(all(pred_sem <= 1), TRUE)
  expect_true(all(pred_sem >= 0), TRUE)
})

test_that("predict works with categorical data", {
  x <- matrix(runif(150), nrow = 50, ncol = 3)
  cuts <- seq(0, 1, length.out = 4)
  xd <- apply(x, 2, function(col) as.numeric(cut(col, cuts)))
  theta <- t(matrix(c(0, 0, 0, 2, 2, 2, -2, -2, -2), ncol = 3, nrow = 3))
  log_odd <- rowSums(t(sapply(seq_along(xd[, 1]), function(row_id) sapply(seq_along(xd[row_id, ]), function(element) theta[xd[row_id, element], element]))))
  y <- rbinom(50, 1, 1 / (1 + exp(-log_odd)))
  sem_disc <- glmdisc(data.frame(x.1 = factor(xd[, 1]), x.2 = factor(xd[, 2]), x.3 = factor(xd[, 3])), y, iter = 20, m_start = 4, test = FALSE, validation = FALSE, criterion = "aic")
  pred_sem <- predict(sem_disc, rbind(data.frame(x.1 = factor(xd[, 1]), x.2 = factor(xd[, 2]), x.3 = factor(xd[, 3]))))

  expect_length(pred_sem, 50)
  expect_true(all(pred_sem <= 1), TRUE)
  expect_true(all(pred_sem >= 0), TRUE)
})

test_that("predict errors with other data columns", {
  x <- matrix(runif(150), nrow = 50, ncol = 3)
  cuts <- seq(0, 1, length.out = 4)
  xd <- apply(x, 2, function(col) as.numeric(cut(col, cuts)))
  theta <- t(matrix(c(0, 0, 0, 2, 2, 2, -2, -2, -2), ncol = 3, nrow = 3))
  log_odd <- rowSums(t(sapply(seq_along(xd[, 1]), function(row_id) sapply(seq_along(xd[row_id, ]), function(element) theta[xd[row_id, element], element]))))
  y <- rbinom(50, 1, 1 / (1 + exp(-log_odd)))
  sem_disc <- glmdisc(data.frame(x.1 = factor(xd[, 1]), x.2 = factor(xd[, 2]), x.3 = factor(xd[, 3])), y, iter = 20, m_start = 4, test = FALSE, validation = FALSE, criterion = "aic")
  expect_error(predict(sem_disc, x[, 1:2]))
})

test_that("predict errors with other data types", {
  x <- matrix(runif(150), nrow = 50, ncol = 3)
  cuts <- seq(0, 1, length.out = 4)
  xd <- apply(x, 2, function(col) as.numeric(cut(col, cuts)))
  theta <- t(matrix(c(0, 0, 0, 2, 2, 2, -2, -2, -2), ncol = 3, nrow = 3))
  log_odd <- rowSums(t(sapply(seq_along(xd[, 1]), function(row_id) sapply(seq_along(xd[row_id, ]), function(element) theta[xd[row_id, element], element]))))
  y <- rbinom(50, 1, 1 / (1 + exp(-log_odd)))
  sem_disc <- glmdisc(data.frame(x.1 = factor(xd[, 1]), x.2 = factor(xd[, 2]), x.3 = factor(xd[, 3])), y, iter = 20, m_start = 4, test = FALSE, validation = FALSE, criterion = "aic")
  expect_error(predict(sem_disc, x))
})

test_that("predict errors with unseen categorical data", {
  x <- matrix(runif(150), nrow = 50, ncol = 3)
  cuts <- seq(0, 1, length.out = 4)
  xd <- apply(x, 2, function(col) as.numeric(cut(col, cuts)))
  theta <- t(matrix(c(0, 0, 0, 2, 2, 2, -2, -2, -2), ncol = 3, nrow = 3))
  log_odd <- rowSums(t(sapply(seq_along(xd[, 1]), function(row_id) sapply(seq_along(xd[row_id, ]), function(element) theta[xd[row_id, element], element]))))
  y <- rbinom(50, 1, 1 / (1 + exp(-log_odd)))
  sem_disc <- glmdisc(data.frame(x.1 = factor(xd[, 1]), x.2 = factor(xd[, 2]), x.3 = factor(xd[, 3])), y, iter = 20, m_start = 4, test = FALSE, validation = FALSE, criterion = "aic")
  expect_error(predict(sem_disc, rbind(
    data.frame(x.1 = factor(xd[, 1]), x.2 = factor(xd[, 2]), x.3 = factor(xd[, 3])),
    data.frame(x.1 = factor(10), x.2 = factor(10), x.3 = factor(10))
  )))
})

# test_that("predict warns when removing levels", {
#   x <- matrix(runif(150), nrow = 50, ncol = 3)
#   cuts <- seq(0, 1, length.out = 4)
#   xd <- apply(x, 2, function(col) as.numeric(cut(col, cuts)))
#   theta <- t(matrix(c(0, 0, 0, 2, 2, 2, -2, -2, -2), ncol = 3, nrow = 3))
#   log_odd <- rowSums(t(sapply(seq_along(xd[, 1]), function(row_id) sapply(seq_along(xd[row_id, ]), function(element) theta[xd[row_id, element], element]))))
#   y <- rbinom(50, 1, 1 / (1 + exp(-log_odd)))
#   sem_disc <- glmdisc(data.frame(x.1 = factor(xd[, 1]), x.2 = factor(xd[, 2]), x.3 = factor(xd[, 3])), y, iter = 20, m_start = 4, test = FALSE, validation = FALSE, criterion = "aic")
#   sem_disc@parameters$encoder$lvls$V1 <- c("1", "2")
#   sem_disc@best.disc[[1]]$coefficients <- sem_disc@best.disc[[1]]$coefficients[1:(length(sem_disc@best.disc[[1]]$coefficients) - 1)]
#   expect_warning(predict(sem_disc, data.frame(x.1 = factor(xd[, 1]), x.2 = factor(xd[, 2]), x.3 = factor(xd[, 3]))))
# })
