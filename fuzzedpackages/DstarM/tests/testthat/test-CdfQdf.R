context("estCdf & estQdf")

tol <- 1e-4
xx <- seq(-5, 5, length.out = 1e5)

test_that("estCdf works.", {
  approxCdf <- estCdf(dnorm(xx))
  trueCdf<- pnorm(xx)
  expect_equal(
    object = approxCdf,
    expected = trueCdf,
    tolerance = tol
  )
})

test_that("estQdf works", {
  d <- dnorm(xx)
  p <- seq(0, 1, length.out = 21)
  cEst <- estCdf(d)
  approxQdf <- estQdf(p = p, x = xx, cdf = cEst)
  trueQdf <- qnorm(p)
  expect_equal(
    object = approxQdf[-c(1, length(approxQdf))],
    expected = trueQdf[-c(1, length(trueQdf))],
    tolerance = tol
  )
})

test_that("Normalize works", {

  d <- dnorm(xx)
  d0 <- normalize(d, xx)
  expect_equal(
    object = DstarM:::simpson(xx, d),
    expected = 0.9999994,
    label = "simpson integrates to approximately 1.",
    tolerance = tol
  )
  expect_equal(
    object = DstarM:::simpson(xx, d0),
    expected = 1,
    label = "After normalization integrates to 1.",
    tolerance = tol
  )
  d1 <- normalize(d, xx, props = .5)
  expect_equal(
    object = DstarM:::simpson(xx, d1),
    expected = .5,
    label = "After normalization integrates to 0.5.",
    tolerance = tol
  )
})
