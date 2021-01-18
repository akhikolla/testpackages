context('arms_gibbs')

log_dnorm <- function(...) sum(dnorm(..., log = TRUE))

test_that('arms_gibbs returns the same result for the same seed', {
  run_arms_gibbs <- function() {
    arms_gibbs(10, c(0, 0), log_dnorm, -10, 10, metropolis = FALSE)
  }
  set.seed(1)
  output1 <- run_arms_gibbs()
  set.seed(1)
  output2 <- run_arms_gibbs()

  expect_equal(output1, output2)
})

test_that('arms_gibbs returns output of the expected length', {
  n_samples <- 10
  n_dim <- 2
  output <- arms_gibbs(
    n_samples,
    rep(0, n_dim),
    log_dnorm,
    -10,
    10,
    metropolis = FALSE
  )
  expect_equal(dim(output), c(n_samples, n_dim))
})

test_that('arms_gibbs with include_evaluations has the right structure', {
  n_samples <- 10
  n_dim <- 2
  output <- arms_gibbs(
    n_samples,
    rep(0, n_dim),
    log_dnorm,
    -10,
    10,
    metropolis = FALSE,
    include_n_evaluations = TRUE
  )
  expect_equal(names(output), c('n_evaluations', 'samples'))
  expect_true(output$n_evaluations > 0)
  expect_equal(dim(output$samples), c(n_samples, n_dim))
})

test_that('arms_gibbs accepts show_progress', {
  n_samples <- 10
  n_dim <- 2
  captured_output <- capture.output({
    output <- arms_gibbs(
      n_samples,
      rep(0, n_dim),
      log_dnorm,
      -10,
      10,
      metropolis = FALSE,
      show_progress = TRUE
    )
  })
  expect_true(nchar(captured_output) > 0)
  expect_equal(dim(output), c(n_samples, n_dim))
})
