context('arms')

log_dnorm <- function(...) dnorm(..., log = TRUE)

log_dnormmixture <- function(x) {
  parts <- log(c(0.4, 0.6)) + dnorm(x, mean = c(-1, 4), log = TRUE)
  log(sum(exp(parts - max(parts)))) + max(parts)
}

test_that('arms returns the same result for the same seed', {
  run_arms <- function() {
    arms(10, log_dnorm, -10, 10, metropolis = FALSE)
  }
  set.seed(1)
  output1 <- run_arms()
  set.seed(1)
  output2 <- run_arms()

  expect_equal(output1, output2)
})

test_that('arms returns output of the expected length', {
  n_samples <- 10
  output <- arms(n_samples, log_dnorm, -10, 10, metropolis = FALSE)
  expect_equal(length(output), n_samples)
})

test_that('arms with include_evaluations has the right structure', {
  n_samples <- 10
  output <- arms(
    n_samples,
    log_dnorm,
    -10,
    10,
    metropolis = FALSE,
    include_n_evaluations = TRUE
  )
  expect_equal(names(output), c('n_evaluations', 'samples'))
  expect_true(output$n_evaluations > 0)
  expect_equal(length(output$samples), n_samples)
})

test_that('arms with metropolis gives output with the expected length', {
  n_samples <- 10
  output <- arms(
    n_samples,
    log_dnormmixture,
    -100,
    100,
    previous = 0,
    metropolis = TRUE
  )
  expect_equal(length(output), n_samples)
})

test_that('arms accepts arguments passed to log pdf', {
  run_arms <- function(arguments) {
    n_samples <- 10
    output <- arms(
      n_samples,
      log_dnorm,
      -100,
      100,
      arguments = arguments,
      metropolis = FALSE
    )
    expect_equal(length(output), n_samples)
  }

  # Named argument
  run_arms(list(sd = 1))

  # Unnamed argument
  run_arms(list(1))

  # Non-vector argument
  output <- arms(
    1,
    function(x, y) {
      stopifnot(is.matrix(y))
      dnorm(x, log = TRUE)
    },
    -100,
    100,
    arguments = list(y = matrix(0, nrow = 2, ncol = 2)),
    metropolis = FALSE
  )
  expect_equal(length(output), 1)
})

test_that('arms accepts initial argument', {
  n_samples <- 10
  output <- arms(
    n_samples,
    log_dnorm,
    -100,
    100,
    initial = c(-5, 0, 5),
    metropolis = FALSE
  )
  expect_equal(length(output), n_samples)
})

test_that('arms recycles arguments', {
  samples <- matrix(
    arms(200, log_dnorm, -1000, 1000, metropolis = FALSE, arguments = list(
      mean = c(0, 10),
      sd = 1
    )),
    ncol = 2,
    byrow = TRUE
  )
  sample_means <- colMeans(samples)

  # Distribution of the sample mean has sd = 1 / 10, so this should be
  # astronomically unlikely to fail by chance
  expect_true(all(sample_means > c(-1, 9) & sample_means < c(1, 11)))
})
