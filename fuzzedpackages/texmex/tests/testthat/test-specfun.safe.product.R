context("specfun.safe.product")

test_that("specfun.safe.product behaves as it should", {
  skip_on_cran()
  skip_on_travis()

  prod <- .specfun.safe.product

  ## simple values
  x <- runif(10, -5, 5)
  y <- runif(length(x), -5, 5)

  expect_equal(prod(x, y), pmax(x * y, -1),
               label="safe product: simple values")

  ## complicated values
  x <- c(0, 0, 1, 1)
  y <- c(Inf, -Inf, Inf, -Inf)

  res <- c(0, 0, Inf, -1)
  expect_equal(prod(x, y), res,
               label="safe product: complicated values")
}
)

test_that("specfun.safe.product correctly recycles", {
  skip_on_cran()
    set.seed(123456)
    x <- runif(10, -5, 5)
    y <- c(0.5)

    expect_equal(.specfun.safe.product(x, y),
                 pmax(x * y, -1))
})

test_that("specfun.safe.product correctly recycles again", {
  skip_on_cran()
    x <- 0.5
    y <- c(1,2,3,4)

    expect_equal(.specfun.safe.product(x, y),
                 pmax(x * y, -1))

})
