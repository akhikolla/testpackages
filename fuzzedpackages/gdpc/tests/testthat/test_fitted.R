library(gdpc)
context("Test that fitted method gets mses right")

set.seed(1234)
N <- 50 #length of series
m <- 10 #number of series
set.seed(1234)
f <- rnorm(N + 1)
x_small <- matrix(0, N, m)
u <- matrix(rnorm(N * m), N, m)
for (i in 1:m) {
  x_small[, i] <- 10 * sin(2 * pi * (i/m)) * f[1:N] + 10 * cos(2 * pi * (i/m)) * f[2:(N + 1)] + u[, i]
}

N <- 70 #length of series
m <- 100 #number of series
set.seed(1234)
f <- rnorm(N + 1)
x_large <- matrix(0, N, m)
u <- matrix(rnorm(N * m), N, m)
for (i in 1:m) {
  x_large[, i] <- 10 * sin(2 * pi * (i/m)) * f[1:N] + 10 * cos(2 * pi * (i/m)) * f[2:(N + 1)] + u[, i]
}

fit1_small <- gdpc(x_small, k=1)
pred1_small <- fitted(fit1_small)
resid1_small <- x_small - pred1_small
fit2_small <- gdpc(resid1_small, k=2)
pred2_small <- fitted(fit2_small)
resid2_small <- resid1_small - pred2_small

test_that(paste0("Correct mse, first component for m<T"), {
  mse1 <- mean(resid1_small^2)
  expect_equal(mse1, fit1_small$mse, tolerance=1e-2)
})
test_that(paste0("Correct mse, second component for m<T"), {
  mse2 <- mean(resid2_small^2)
  expect_equal(mse2, fit2_small$mse, tolerance=1e-2)
})

fit1_large <- gdpc(x_large, k=1)
pred1_large <- fitted(fit1_large)
resid1_large <- x_large - pred1_large
fit2_large <- gdpc(resid1_large, k=2)
pred2_large <- fitted(fit2_large)
resid2_large <- resid1_large - pred2_large

test_that(paste0("Correct mse, first component for m<T"), {
  mse1 <- mean(resid1_large^2)
  expect_equal(mse1, fit1_large$mse, tolerance=1e-2)
})
test_that(paste0("Correct mse, second component for m<T"), {
  mse2 <- mean(resid2_large^2)
  expect_equal(mse2, fit2_large$mse, tolerance=1e-2)
})