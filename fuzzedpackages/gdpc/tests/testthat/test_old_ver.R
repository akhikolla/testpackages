library(gdpc)
context("Compare gdpc output with stored values")

T <- 50 #length of series
m <- 10 #number of series
set.seed(1234)
f <- rnorm(T + 1)
x_small <- matrix(0, T, m)
u <- matrix(rnorm(T * m), T, m)
for (i in 1:m) {
  x_small[, i] <- 10 * sin(2 * pi * (i/m)) * f[1:T] + 10 * cos(2 * pi * (i/m)) * f[2:(T + 1)] + u[, i]
}
# from gdpc 1.0.3
old_small_mse_gdpc <- 0.9171641
old_small_mse_auto.gdpc <- 0.9171641
old_small_k_auto.gdpc <- 1
old_small_expvar_auto.gdpc <- 0.9888407

fit <- gdpc(x_small, k = 1)

test_that('Equality for gdpc with stored mse, small data', {
  expect_equal(fit$mse, old_small_mse_gdpc, tolerance=1e-2)
})

aut <- auto.gdpc(x_small)

test_that('Equality for auto.gdpc with stored mse, small data', {
  expect_equal(aut[[1]]$mse, old_small_mse_auto.gdpc, tolerance=1e-2)
})

test_that('Equality for auto.gdpc with stored k, small data', {
  expect_equal(aut[[1]]$k, old_small_k_auto.gdpc, tolerance=1e-2)
})

test_that('Equality for auto.gdpc with stored explained variance, small data', {
  expect_equal(aut[[1]]$expart, old_small_expvar_auto.gdpc, tolerance=1e-2)
})

T <- 50 #length of series
m <- 100 #number of series
set.seed(1234)
f <- rnorm(T + 1)
x_large <- matrix(0, T, m)
u <- matrix(rnorm(T * m), T, m)
for (i in 1:m) {
  x_large[, i] <- 10 * sin(2 * pi * (i/m)) * f[1:T] + u[, i]
}
# from gdpc 1.0.3
old_large_mse_gdpc <- 0.9086348
old_large_mse_auto.gdpc <- 0.9333897
old_large_k_auto.gdpc <- 0
old_large_expvar_auto.gdpc <- 0.9765623

fit <- gdpc(x_large, k = 1)

test_that('Equality for gdpc with stored mse, large data', {
  expect_equal(fit$mse, old_large_mse_gdpc, tolerance=1e-2)
})

aut <- auto.gdpc(x_large)

test_that('Equality for auto.gdpc with stored mse, large data', {
  expect_equal(aut[[1]]$mse, old_large_mse_auto.gdpc, tolerance=1e-2)
})

test_that('Equality for auto.gdpc with stored k, large data', {
  expect_equal(aut[[1]]$k, old_large_k_auto.gdpc, tolerance=1e-2)
})

test_that('Equality for auto.gdpc with stored explained variance, large data', {
  expect_equal(aut[[1]]$expart, old_large_expvar_auto.gdpc, tolerance=1e-2)
})
