context("mvrsquared tests")

### Test univariate rsquared ----
f <- stats::lm(mpg ~ cyl + disp + hp + wt, data = datasets::mtcars)

y <- f$model$mpg

yhat <- f$fitted.values

s <- summary(f)

test_that("We get the expected value for correct inputs to univariate rsquared",{

  r2 <- calc_rsquared(y = y, yhat = yhat)

  expect_equal(round(r2, 3), round(s$r.squared, 3))

  ss <- calc_rsquared(y = y, yhat = yhat, ybar = mean(y), return_ss_only = TRUE)

  expect_equal(length(ss), 2)

  expect_equal(r2, 1 - ss[[1]] / ss[[2]])

})


test_that("Get the right r-squared for single column matrix inputs", {

  y_mat <- matrix(y, ncol = 1)

  x <- cbind(1, as.matrix(f$model[, -1]))

  w <- matrix(s$coefficients[, 1], ncol = 1)

  expect_equal(calc_rsquared(y = y_mat, yhat = list(x, w)),
                         calc_rsquared(y = y, yhat = yhat))

})

test_that("Multithreading works as expected",{

  r2 <- calc_rsquared(y = y, yhat = yhat, threads = 2)

  expect_equal(length(r2), 1)

  # expect_equal(round(r2, 3), round(s$r.squared, 3))

  ss <- calc_rsquared(y = y, yhat = yhat, ybar = mean(y), return_ss_only = TRUE)

  expect_equal(length(ss), 2)

  # expect_equal(r2, 1 - ss[[1]] / ss[[2]])

})


### fancier stuff ----


test_that("can pass named 'w' and 'x' in list for 'yhat' out-of order and still get the same calculation", {

  x <- cbind(1, as.matrix(f$model[, -1]))

  w <- matrix(s$coefficients[, 1], ncol = 1)

  r2_1 <- calc_rsquared(y = y, yhat = list(w = w, x = x))

  expect_equal(r2_1, s$r.squared)

  # name only one and you should get a warning
  expect_warning(calc_rsquared(y = y, yhat = list(x = x, w)))

  # repeated names of 'x' or 'w' produce a warning
  expect_warning(calc_rsquared(y = y, yhat = list(x = x, w = w, x = x)))

  # confirm that naming nothing produces no warning
  calc_rsquared(y = y, yhat = list(x, w))

})

test_that("Get the right value for dgCMatrix inputs", {

  expect_type(calc_rsquared(y = Matrix::Matrix(y, ncol = 1, sparse = TRUE), yhat = yhat),
              "double")

})

test_that("get errors for incompatible dimensions",{

  # not enough columns for yhat
  expect_error(calc_rsquared(y = cbind(y, y), yhat = yhat))

  # too many columns for yhat
  expect_error(calc_rsquared(y = y, yhat = cbind(yhat, yhat)))

  # number of rows does not match
  expect_error(calc_rsquared(y = y[1:10], yhat = yhat))

  # pass a vector and matrix in list
  expect_error(
    calc_rsquared(y = y, yhat = list(yhat, matrix(1, nrow = 1, ncol = 3)))
    )

  # dimensions of matrices in list do not match
  expect_error(
    calc_rsquared(y = y, yhat = list(matrix(yhat, ncol = 1), matrix(1, nrow = 1, ncol = 3)))
    )

})

test_that("batch (for parallel) computation behaves nicely", {

  # define some batches
  batches <- list(list(y = cbind(y[1:16], y[1:16]), yhat = cbind(yhat[1:16], yhat[1:16])),
                  list(y = cbind(y[17:32], y[17:32]), yhat = cbind(yhat[17:32], yhat[17:32])))

  ybar <- c(mean(y), mean(y))

  # calc sum of squares by batch with correct inputs
  # Note: this uses lapply, but one could easily do this in parallel with mclapply or similar
  ss <- lapply(X = batches,
               FUN = function(ybatch){
                     calc_rsquared(y = ybatch$y, yhat = ybatch$yhat, ybar = ybar, return_ss_only = TRUE)
                   })


  sse <- sum(sapply(ss, function(x) x["sse"]))

  sst <- sum(sapply(ss, function(x) x["sst"]))

  r2_batch <- 1 - sse / sst # final r-squared value here

  expect_equal(round(r2_batch, 4), round(s$r.squared, 4))

  # leave out ybar and get a warning
  expect_warning(lapply(X = batches,
                                  FUN = function(ybatch){
                                    calc_rsquared(y = ybatch$y, yhat = ybatch$yhat, return_ss_only = TRUE)
                                  }))

})

test_that("Errors are triggered for malformed inputs",{

  expect_error(calc_rsquared(y = y, yhat = yhat, return_ss_only = NA))

  expect_error(calc_rsquared(y = y, yhat = yhat, ybar = c(1,1)))

  expect_error(calc_rsquared(y = list(y), yhat = yhat))

  expect_error(calc_rsquared(y = as.character(y), yhat = yhat))

  expect_error(calc_rsquared(y = as.character(y), yhat = yhat))

  expect_error(calc_rsquared(y = y, yhat = data.frame(yhat)))

  expect_error(calc_rsquared(y = y, yhat = 5))


})
