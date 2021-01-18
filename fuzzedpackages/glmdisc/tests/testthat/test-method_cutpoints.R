context("test-method_cutpoints")

test_that("cutpoints works with one-dimensional continuous data", {
  x <- matrix(runif(40), nrow = 40, ncol = 1)
  cuts <- seq(0, 1, length.out = 4)
  xd <- as.numeric(cut(x, cuts))
  xd <- t(t(xd))
  theta <- matrix(c(0, 2, -2), ncol = 1, nrow = 3)
  log_odd <- sapply(
    seq_along(xd[, 1]),
    function(row_id) {
      sapply(
        seq_along(xd[row_id, ]),
        function(element) theta[xd[row_id, element], element]
      )
    }
  )
  y <- rbinom(40, 1, 1 / (1 + exp(-log_odd)))
  sem_disc <- glmdisc(x,
    y,
    iter = 15,
    m_start = 4,
    test = FALSE,
    validation = FALSE,
    criterion = "aic",
    interact = FALSE
  )
  les_cuts <- cutpoints(sem_disc)
  expect_type(les_cuts, "list")
  expect_length(les_cuts, 1)
  expect_type(les_cuts[[1]], "double")
  expect_equal(names(les_cuts), "V1")
  x <- data.frame(x)
  sem_disc <- glmdisc(x,
    y,
    iter = 15,
    m_start = 4,
    test = FALSE,
    validation = FALSE,
    criterion = "aic",
    interact = FALSE
  )
  les_cuts <- cutpoints(sem_disc)
  expect_equal(names(les_cuts), "x")
})

test_that("cutpoints works with multi-dimensional continuous data", {
  x <- matrix(runif(120), nrow = 40, ncol = 3)
  cuts <- seq(0, 1, length.out = 4)
  xd <- apply(x, 2, function(col) as.numeric(cut(col, cuts)))
  theta <- t(matrix(c(0, 0, 0, 2, 2, 2, -2, -2, -2), ncol = 3, nrow = 3))
  log_odd <- rowSums(t(sapply(seq_along(xd[, 1]), function(row_id) sapply(seq_along(xd[row_id, ]), function(element) theta[xd[row_id, element], element]))))
  y <- rbinom(40, 1, 1 / (1 + exp(-log_odd)))
  sem_disc <- glmdisc(x, y, iter = 15, m_start = 4, test = FALSE, validation = FALSE, criterion = "aic")
  les_cuts <- cutpoints(sem_disc)
  expect_type(les_cuts, "list")
  expect_length(les_cuts, 3)
  expect_type(les_cuts[[1]], "double")
  expect_equal(names(les_cuts), c("V1", "V2", "V3"))
})

test_that("cutpoints works with multi-dimensional categorical data", {
  x <- matrix(runif(120), nrow = 40, ncol = 3)
  cuts <- seq(0, 1, length.out = 4)
  xd <- apply(x, 2, function(col) as.numeric(cut(col, cuts)))
  theta <- t(matrix(c(0, 0, 0, 2, 2, 2, -2, -2, -2), ncol = 3, nrow = 3))
  log_odd <- rowSums(t(sapply(seq_along(xd[, 1]), function(row_id) sapply(seq_along(xd[row_id, ]), function(element) theta[xd[row_id, element], element]))))
  y <- rbinom(40, 1, 1 / (1 + exp(-log_odd)))
  sem_disc <- glmdisc(data.frame(apply(xd, 2, factor), stringsAsFactors = TRUE), y, iter = 15, m_start = 4, test = FALSE, validation = FALSE, criterion = "aic")
  les_cuts <- cutpoints(sem_disc)
  expect_type(les_cuts, "list")
  expect_length(les_cuts, 3)
  expect_type(les_cuts[[1]], "integer")
  expect_equal(names(les_cuts), c("X1", "X2", "X3"))
})
