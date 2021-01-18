test_that("plot works", {
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
  sem_disc <- glmdisc(data.frame(x),
    y,
    iter = 15,
    m_start = 4,
    test = FALSE,
    validation = FALSE,
    criterion = "aic",
    interact = FALSE
  )
  plot(sem_disc)
  sem_disc <- glmdisc(data.frame(x = factor(xd)),
    y,
    iter = 15,
    m_start = 4,
    test = FALSE,
    validation = FALSE,
    criterion = "aic",
    interact = FALSE
  )
  plot(sem_disc)
})
