context("multivariate normal samplers")

test_that("rmvn_arma", {
    ## check that the rmvn_arma function generates draws from the
    ## same distribution -- verify Monte Carlo mean and covariance
    ## are the same up to Monte Carlo error
    set.seed(11)
    N <- 10000
    d <- 6
    b <- rnorm(6)
    Z <- matrix(rnorm((d + 2) * d), d + 2, d)
    A <- t(Z) %*% Z

    X <- mvnfast::rmvn(N, solve(A) %*% b, solve(A))
    Y <- t(sapply(1:N, function(i) rmvn_arma(A, b)))
    expect_equal(
        colMeans(X), colMeans(Y), tol = 0.05
    )
    expect_equal(
        var(X), var(Y), tol = 0.05
    )

    set.seed(11)
    N <- 10000
    d <- 3
    b <- rnorm(3)
    ## matrix with third eigenvalue of 0
    A <- matrix(c(1, 0, 2, 0, 1, 2, 2, 2, 8), 3, 3)

    X <- mvnfast::rmvn(N, solve(A + 1e-6 * diag(d)) %*% b, solve(A + 1e-6 * diag(d)))
    Y <- t(sapply(1:N, function(i) rmvn_arma(A, b)))
    expect_equal(
        colMeans(X), colMeans(Y), tol = 0.05
    )
    expect_equal(
        var(X), var(Y), tol = 0.05
    )
})


test_that("rmvn_arma_chol", {
    ## check that the rmvn_arma function generates draws from the
    ## same distribution -- verify Monte Carlo mean and covariance
    ## are the same up to Monte Carlo error
    set.seed(11)
    N <- 10000
    d <- 6
    b <- rnorm(6)
    Z <- matrix(rnorm((d + 2) * d), d + 2, d)
    A <- t(Z) %*% Z
    A_chol <- chol(A)

    X <- mvnfast::rmvn(N, solve(A) %*% b, solve(A))
    Y <- t(sapply(1:N, function(i) rmvn_arma_chol(A_chol, b)))
    expect_equal(
        colMeans(X), colMeans(Y), tol = 0.05
    )
    expect_equal(
        var(X), var(Y), tol = 0.05
    )
})


test_that("rmvn_arma_scalar", {
    ## check that the rmvn_arma_scalar function generates draws from the
    ## same distribution -- verify Monte Carlo mean and covariance
    ## are the same up to Monte Carlo error
    set.seed(11)
    N <- 10000
    b <- rnorm(1)
    a <- rgamma(1, 1, 1)

    X <- rnorm(N, b / a, sqrt(1 / a))
    Y <- sapply(1:N, function(i) rmvn_arma_scalar(a, b))
    expect_equal(
        mean(X), mean(Y), tol = 0.05
    )
    expect_equal(
        var(X), var(Y), tol = 0.05
    )
})
