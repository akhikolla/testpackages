test_that("discretize works", {
  set.seed(1)
  x <- matrix(runif(80), nrow = 80, ncol = 1)
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
  y <- rbinom(80, 1, 1 / (1 + exp(-log_odd)))
  sem_disc <- glmdisc(x,
    y,
    iter = 25,
    m_start = 4,
    test = FALSE,
    validation = FALSE,
    criterion = "aic",
    interact = FALSE
  )
  expect_length(discretize(sem_disc, x), 80)
})
