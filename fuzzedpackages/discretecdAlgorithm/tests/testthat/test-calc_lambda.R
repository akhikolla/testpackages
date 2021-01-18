context("calc_lambda")

# set up input variable
data_matrix <- matrix(c(1, 1, 0, 0, 1, 1,
                        1, 1, 0, 1, 1, 1,
                        0, 0, 1, 0, 0, 1,
                        0, 0, 1, 0, 0, 1,
                        0, 0, 0, 1, 1, 0,
                        0, 0, 0, 1, 1, 1,
                        1, 1, 1, 1, 0, 0,
                        1, 0, 1, 1, 0, 1,
                        0, 0, 0, 0, 1, 0,
                        1, 1, 1, 1, 0, 1,
                        1, 1, 0, 1, 1, 1,
                        0, 0, 1, 0, 0, 1,
                        1, 1, 0, 1, 0, 0,
                        1, 0, 1, 1, 0, 1,
                        1, 1, 1, 1, 1, 0,
                        1, 0, 1, 1, 1, 1,
                        0, 0, 1, 0, 0, 0,
                        1, 1, 0, 1, 1, 1,
                        1, 1, 1, 0, 0, 0,
                        1, 1, 1, 1, 0, 0,
                        0, 0, 0, 0, 1, 0,
                        0, 0, 1, 0, 0, 0,
                        1, 0, 0, 0, 1, 1,
                        0, 0, 1, 0, 0, 0,
                        1, 0, 1, 1, 0, 1,
                        0, 0, 0, 1, 1, 0,
                        0, 0, 0, 0, 0, 0,
                        0, 0, 1, 0, 1, 0,
                        0, 0, 1, 0, 0, 0,
                        0, 0, 1, 0, 0, 0), byrow = TRUE, ncol = 6)
data_matrix <- matrix(as.integer(data_matrix), ncol = 6)
node <- ncol(data_matrix); node <- as.integer(node)
dataSize <- nrow(data_matrix); dataSize <- as.integer(dataSize)
n_levels <- rep(2, node)
n_levels <- as.integer(n_levels)
obs <- 1:30
obsIndex_R <- vector("list", length = 6)
obsIndex_R <- lapply(obsIndex_R, function(x, obs){as.integer(obs-1)}, obs)
weights <- matrix(rep(1, 6*6), ncol = 6)
gamma <- 1; gamma <- as.numeric(gamma)
upperbound <- 100; upperbound <- as.numeric(upperbound)

# test
test_that("CD_path runs as expected", {
  ### throw error if parameter and initial values not explicitly specified
  expect_error(calc_lambda(node, dataSize, data_matrix, n_levels))

  ### no error
  expect_error(calc_lambda(node, dataSize, data_matrix, n_levels, obsIndex_R, weights, gamma, upperbound), NA)
})

test_that("Check input: node", {
  ### Throw an error if node is not an integer
  expect_error(calc_lambda(node = 6, dataSize, data_matrix, n_levels, obsIndex_R, weights, gamma, upperbound))

  ### Throw an error if node <= 0
  expect_error(calc_lambda(node = as.integer(-1), dataSize, data_matrix, n_levels, obsIndex_R, weights, gamma, upperbound))
})

test_that("Check input: dataSize", {
  ### throw an error if dataSize is not an integer
  expect_error(calc_lambda(node, dataSize = 30, data_matrix, n_levels, obsIndex_R, weights, gamma, upperbound))

  ### Throw an error if dataSize <= 0
  expect_error(calc_lambda(node, dataSize = as.integer(-10), data_matrix, n_levels, obsIndex_R, weights, gamma, upperbound))
})

test_that("Ckeck input: data_matrix", {
  ### throw an error if the input data_matrix is not an integer matrix.
  data_num <- data_matrix+1
  data_num <- data_num-1
  expect_error(calc_lambda(node, dataSize, data_matrix = data_num, n_levels, obsIndex_R, weights, gamma, upperbound))

  ### throw an error if the input data_matrix has a wrong number of rows.
  data_wrong_dataSize <- data_matrix[1:10, ]
  expect_error(calc_lambda(node, dataSize, data_matrix = data_wrong_dataSize, n_levels, obsIndex_R, weights, gamma, upperbound))

  ### throw an error if the input data_matrix has the wrong number of columns.
  data_wrong_node <- data_matrix[, 1:3]
  expect_error(calc_lambda(node, dataSize, data_matrix = data_wrong_node, n_levels, obsIndex_R, weights, gamma, upperbound))
})

test_that("Check input: n_levels", {
  ### throw an error if the input n_levels is not an integer vector
  expect_error(calc_lambda(node, dataSize, data_matrix, as.numeric(n_levels), obsIndex_R, weights, gamma, upperbound))

  ### throw an error if the input n_levels has the wrong dimension
  n_level_wrong_dim <- rep(2, 10)
  expect_error(calc_lambda(node, dataSize, data_matrix, n_levels = n_level_wrong_dim, obsIndex_R, weights, gamma, upperbound))

  ### throw an error if n_level is less than 2
  n_level_one <- n_levels
  n_level_one[1] <- 1
  expect_error(CD_path(node, dataSize, data_matrix, n_levels = n_level_one, obsIndex_R, eor_nr, eor, lambda_seq, nlam, eps, convLb, qtol, weights, gamma, upperbound, threshold))

  ### throw an error if n_level has less level than data set
  data_matrix_more_level <- data_matrix
  data_matrix_more_level[1, 1] <- 2L
  expect_error(CD_path(node, dataSize, data_matrix = data_matrix_more_level, n_levels = n_levels, obsIndex_R, eor_nr, eor, lambda_seq, nlam, eps, convLb, qtol, weights, gamma, upperbound, threshold))

})

test_that("Check input: obsIndex_R", {
  obsIndex_noList <- 1:100
  obsIndex_noInteger <- lapply(obsIndex_R, function(x){as.numeric(x)})
  obsIndex_negative <- obsIndex_R
  obsIndex_negative[[1]][1] <- -1
  obsIndex_wrong_dim <- obsIndex_R[1:4]

  ### throw an error if obsIndex_R is not a list
  expect_error(calc_lambda(node, dataSize, data_matrix, n_levels, obsIndex_R = obsIndex_noList, weights, gamma, upperbound))

  ### throw an error if element in obsIndex_R is not a integer vector
  expect_error(calc_lambda(node, dataSize, data_matrix, n_levels, obsIndex_R = obsIndex_noInteger, weights, gamma, upperbound))

  ### throw an error if element in obsIndex_R is negative
  expect_error(calc_lambda(node, dataSize, data_matrix, n_levels, obsIndex_R = obsIndex_negative, weights, gamma, upperbound))

  ### throw an error if obsIndex_R has wrong dimension
  expect_error(calc_lambda(node, dataSize, data_matrix, n_levels, obsIndex_R = obsIndex_wrong_dim, weights, gamma, upperbound))
})

test_that("Check input: weights", {
  weights_noNumeric <- matrix(as.character(weights), nrow = node)
  weights_integer <- matrix(as.integer(weights), nrow = node)
  weights_wrong_dim <- matrix(1, nrow=node-1, ncol = node+1)

  ### throw an error if weights is not a numeric matrix
  expect_error(calc_lambda(node, dataSize, data_matrix, n_levels, obsIndex_R, weights = weights_noNumeric, gamma, upperbound))

  ### throw an error if weights is an integer matrix
  expect_error(calc_lambda(node, dataSize, data_matrix, n_levels, obsIndex_R, weights = weights_integer, gamma, upperbound))

  ### throw an error if weights has a wrong dimension
  expect_error(calc_lambda(node, dataSize, data_matrix, n_levels, obsIndex_R, weights = weights_wrong_dim, gamma, upperbound))
})

test_that("Check input: gamma", {
  ### throw an error if gamma is not a numeric number
  expect_error(calc_lambda(node, dataSize, data_matrix, n_levels, obsIndex_R, weights, gamma = as.character(1), upperbound))

  ### throw an error is gamma is a non-positive number
  expect_error(calc_lambda(node, dataSize, data_matrix, n_levels, obsIndex_R, weights, gamma = -1, upperbound))
})

test_that("Check input: upperbound", {
  ### throw an error if upperbound is not a numeric number
  expect_error(calc_lambda(node, dataSize, data_matrix, n_levels, obsIndex_R, weights, gamma, upperbound = as.character(100)))

  ### throw an error if upperbound is negative and it is not -1
  expect_error(calc_lambda(node, dataSize, data_matrix, n_levels, obsIndex_R, weights, gamma, upperbound = -100))

  ### if upper bound is -1, run as usual
  expect_error(calc_lambda(node, dataSize, data_matrix, n_levels, obsIndex_R, weights, gamma, upperbound = -1), NA)
})
