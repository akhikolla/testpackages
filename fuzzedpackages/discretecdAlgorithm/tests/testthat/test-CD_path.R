context("CD_path")

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
eor_nr <- node*(node-1)/2
eor <- matrix(0, nrow=eor_nr, ncol=2)
cnt1 <- 1
for (i in 1:(node-1)) {
  for (j in (i+1):node) {
    eor[cnt1, 1] <- i;
    eor[cnt1, 2] <- j;
    cnt1 = cnt1+1;
  }
}
eor_nr <- as.integer(eor_nr)
eor <- matrix(as.integer(eor), ncol = 2)
nlam <- 30; nlam <- as.integer(nlam)
eps <- 0.0001; eps <- as.numeric(eps)
convLb <- 0.01; convLb <- as.numeric(convLb)
qtol <- 0.0001; qtol <- as.numeric(qtol)
weights <- matrix(rep(1, 6*6), ncol = 6)
gamma <- 1; gamma <- as.numeric(gamma)
fmlam <- 0.1
upperbound <- 100; upperbound <- as.numeric(upperbound)
threshold <- 5; threshold <- as.integer(threshold)
ivn <- vector("list", length = dataSize)
ivn <- lapply(ivn, function(x){return(as.integer(0))})
databn <- sparsebnUtils::sparsebnData(as.data.frame(data_matrix), ivn = ivn, type = "discrete")
lambda_m <- max_lambda(databn, weights, gamma, upperbound)
lambda_seq <- sparsebnUtils::generate.lambdas(lambda.max = lambda_m, lambdas.ratio = fmlam, lambdas.length = nlam)
lambda_seq <- as.numeric(lambda_seq)

# test
test_that("CD_path runs as expected", {
  ### throw error if parameter and initial values not explicitly specified
  expect_error(CD_path(node, dataSize, data_matrix, n_levels))

  ### no error
  expect_error(CD_path(node, dataSize, data_matrix, n_levels, obsIndex_R, eor_nr, eor, lambda_seq, nlam, eps, convLb, qtol, weights, gamma, upperbound, threshold), NA)
})

test_that("Check input: node", {
  ### Throw an error if node is not an integer
  expect_error(CD_path(node = 6, dataSize, data_matrix, n_levels, obsIndex_R, eor_nr, eor, lambda_seq, nlam, eps, convLb, qtol, weights, gamma, upperbound, threshold))

  ### Throw an error if node <= 0
  expect_error(CD_path(node = as.integer(-1), dataSize, data_matrix, n_levels, obsIndex_R, eor_nr, eor, lambda_seq, nlam, eps, convLb, qtol, weights, gamma, upperbound, threshold))
})

test_that("Check input: dataSize", {
  ### throw an error if dataSize is not an integer
  expect_error(CD_path(node, dataSize = 30, data_matrix, n_levels, obsIndex_R, eor_nr, eor, lambda_seq, nlam, eps, convLb, qtol, weights, gamma, upperbound, threshold))

  ### Throw an error if dataSize <= 0
  expect_error(CD_path(node, dataSize = as.integer(-10), data_matrix, n_levels, obsIndex_R, eor_nr, eor, lambda_seq, nlam, eps, convLb, qtol, weights, gamma, upperbound, threshold))
})

test_that("Ckeck input: data_matrix", {
  ### throw an error if the input data_matrix is not an integer matrix.
  data_num <- data_matrix+1
  data_num <- data_num-1
  expect_error(CD_path(node, dataSize, data_matrix = data_num, n_levels, obsIndex_R, eor_nr, eor, lambda_seq, nlam, eps, convLb, qtol, weights, gamma, upperbound, threshold))

  ### throw an error if the input data_matrix has a wrong number of rows.
  data_wrong_dataSize <- data_matrix[1:10, ]
  expect_error(CD_path(node, dataSize, data_matrix = data_wrong_dataSize, n_levels, obsIndex_R, eor_nr, eor, lambda_seq, nlam, eps, convLb, qtol, weights, gamma, upperbound, threshold))

  ### throw an error if the input data_matrix has the wrong number of columns.
  data_wrong_node <- data_matrix[, 1:3]
  expect_error(CD_path(node, dataSize, data_matrix = data_wrong_node, n_levels, obsIndex_R, eor_nr, eor, lambda_seq, nlam, eps, convLb, qtol, weights, gamma, upperbound, threshold))

  data_oneLevel <- data_matrix
  data_oneLevel[, 1] <- rep(0, nrow(data_matrix))
  expect_error(CD_call(indata = data_oneLevel, eor = NULL, permute = TRUE, weights = NULL, lambda_seq = NULL, fmlam = 0.1, nlam = 30, eps = 0.0001, convLb = 0.01, qtol = 0.0001, gamma = 1, upperbound = 100, threshold = 3))

})

test_that("Check input: n_levels", {
  ### throw an error if the input n_levels is not an integer vector
  expect_error(CD_path(node, dataSize, data_matrix, as.numeric(n_levels), obsIndex_R, eor_nr, eor, lambda_seq, nlam, eps, convLb, qtol, weights, gamma, upperbound, threshold))

  ### throw an error if the input n_levels has the wrong dimension
  n_level_wrong_dim <- rep(2, 10)
  expect_error(CD_path(node, dataSize, data_matrix, n_levels = n_level_wrong_dim, obsIndex_R, eor_nr, eor, lambda_seq, nlam, eps, convLb, qtol, weights, gamma, upperbound, threshold))

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
  expect_error(CD_path(node, dataSize, data_matrix, n_levels, obsIndex_R = obsIndex_noList, eor_nr, eor, lambda_seq, nlam, eps, convLb, qtol, weights, gamma, upperbound, threshold))

  ### throw an error if element in obsIndex_R is not a integer vector
  expect_error(CD_path(node, dataSize, data_matrix, n_levels, obsIndex_R = obsIndex_noInteger, eor_nr, eor, lambda_seq, nlam, eps, convLb, qtol, weights, gamma, upperbound, threshold))

  ### throw an error if element in obsIndex_R is negative
  expect_error(CD_path(node, dataSize, data_matrix, n_levels, obsIndex_R = obsIndex_negative, eor_nr, eor, lambda_seq, nlam, eps, convLb, qtol, weights, gamma, upperbound, threshold))

  ### throw an error if obsIndex_R has wrong dimension
  expect_error(CD_path(node, dataSize, data_matrix, n_levels, obsIndex_R = obsIndex_wrong_dim, eor_nr, eor, lambda_seq, nlam, eps, convLb, qtol, weights, gamma, upperbound, threshold))
})

test_that("Check input: eor_nr and eor", {
  eor_noInteger <- as.numeric(eor)
  eor_empty <- rep()
  eor_nr_noInteger <- as.numeric(eor_nr)
  eor_nr_wrong_dim <- 4

  ### throw an error if element in eor is not integer
  expect_error(CD_path(node, dataSize, data_matrix, n_levels, obsIndex_R, eor_nr, eor = eor_noInteger, lambda_seq, nlam, eps, convLb, qtol, weights, gamma, upperbound, threshold))

  ### throw an error if eor is empty
  expect_error(CD_path(node, dataSize, data_matrix, n_levels, obsIndex_R, eor_nr, eor = eor_empty, lambda_seq, nlam, eps, convLb, qtol, weights, gamma, upperbound, threshold))

  ### throw an error if eor_nr is not integer
  expect_error(CD_path(node, dataSize, data_matrix, n_levels, obsIndex_R, eor_nr = eor_nr_noInteger, eor, lambda_seq, nlam, eps, convLb, qtol, weights, gamma, upperbound, threshold))

  ### throw an error if eor_nr does not equal to the number of row of eor
  expect_error(CD_path(node, dataSize, data_matrix, n_levels, obsIndex_R, eor_nr = eor_nr_wrong_dim, eor, lambda_seq, nlam, eps, convLb, qtol, weights, gamma, upperbound, threshold))
})

test_that("Check input: lambda_seq", {
  lambda_seq_neg <- lambda_seq
  lambda_seq_neg[10] <- -lambda_seq[10]

  ### throw an error if lambda_seq is not a numeric vector
  expect_error(CD_path(node, dataSize, data_matrix, n_levels, obsIndex_R, eor_nr, eor, as.character(lambda_seq), nlam, eps, convLb, qtol, weights, gamma, upperbound, threshold))

  ### throw an error if lambda_seq has a negative element
  expect_error(CD_path(node, dataSize, data_matrix, n_levels, obsIndex_R, eor_nr, eor, lambda_seq = lambda_seq_neg, nlam, eps, convLb, qtol, weights, gamma, upperbound, threshold))
})

test_that("Check input: nlam", {
  ### throw an error if nlam is not an integer
  expect_error(CD_path(node, dataSize, data_matrix, n_levels, obsIndex_R, eor_nr, eor, lambda_seq, nlam = 10, eps, convLb, qtol, weights, gamma, upperbound, threshold))

  ### throw an error if nlam is non-positive
  expect_error(CD_path(node, dataSize, data_matrix, n_levels, obsIndex_R, eor_nr, eor, lambda_seq, nlam = 0L, eps, convLb, qtol, weights, gamma, upperbound, threshold))
  expect_error(CD_path(node, dataSize, data_matrix, n_levels, obsIndex_R, eor_nr, eor, lambda_seq, nlam = -10L, eps, convLb, qtol, weights, gamma, upperbound, threshold))

  ### throw an error if nlam is not the same with length of lambda_seq
  expect_error(CD_path(node, dataSize, data_matrix, n_levels, obsIndex_R, eor_nr, eor, lambda_seq, nlam = 5L, eps, convLb, qtol, weights, gamma, upperbound, threshold))
})

test_that("Check input: eps", {
  ### throw an error if eps is not an numeric number
  expect_error(CD_path(node, dataSize, data_matrix, n_levels, obsIndex_R, eor_nr, eor, lambda_seq, nlam, eps = as.character(0.00001), convLb, qtol, weights, gamma, upperbound, threshold))

  ### throw an error if eps is non-positive
  expect_error(CD_path(node, dataSize, data_matrix, n_levels, obsIndex_R, eor_nr, eor, lambda_seq, nlam, eps = -0.00001, convLb, qtol, weights, gamma, upperbound, threshold))
})

test_that("Check input: convLb", {
  ### throw an error if convLb is not an numeric number
  expect_error(CD_path(node, dataSize, data_matrix, n_levels, obsIndex_R, eor_nr, eor, lambda_seq, nlam, eps, convLb = as.character(0.01), qtol, weights, gamma, upperbound, threshold))

  ### throw an error if convLb is non-positive
  expect_error(CD_path(node, dataSize, data_matrix, n_levels, obsIndex_R, eor_nr, eor, lambda_seq, nlam, eps, convLb = -0.01, qtol, weights, gamma, upperbound, threshold))
})

test_that("Check input: qtol", {
  ### throw an error if qtol is not an numeric number
  expect_error(CD_path(node, dataSize, data_matrix, n_levels, obsIndex_R, eor_nr, eor, lambda_seq, nlam, eps, convLb, qtol = as.character(0.00001), weights, gamma, upperbound, threshold))

  ### throw an error if qtol is non-positive
  expect_error(CD_path(node, dataSize, data_matrix, n_levels, obsIndex_R, eor_nr, eor, lambda_seq, nlam, eps, convLb, qtol = -0.00001, weights, gamma, upperbound, threshold))
})

test_that("Check input: weights", {
  weights_noNumeric <- matrix(as.character(weights), nrow = node)
  weights_integer <- matrix(as.integer(weights), nrow = node)
  weights_wrong_dim <- matrix(1, nrow=node-1, ncol = node+1)

  ### throw an error if weights is not a numeric matrix
  expect_error(CD_path(node, dataSize, data_matrix, n_levels, obsIndex_R, eor_nr, eor, lambda_seq, nlam, eps, convLb, qtol, weights = weights_noNumeric, gamma, upperbound, threshold))

  ### throw an error if weights is an integer matrix
  expect_error(CD_path(node, dataSize, data_matrix, n_levels, obsIndex_R, eor_nr, eor, lambda_seq, nlam, eps, convLb, qtol, weights = weights_integer, gamma, upperbound, threshold))

  ### throw an error if weights has a wrong dimension
  expect_error(CD_path(node, dataSize, data_matrix, n_levels, obsIndex_R, eor_nr, eor, lambda_seq, nlam, eps, convLb, qtol, weights = weights_wrong_dim, gamma, upperbound, threshold))
})

test_that("Check input: gamma", {
  ### throw an error if gamma is not a numeric number
  expect_error(CD_path(node, dataSize, data_matrix, n_levels, obsIndex_R, eor_nr, eor, lambda_seq, nlam, eps, convLb, qtol, weights, gamma = as.character(1), upperbound, threshold))

  ### throw an error is gamma is a non-positive number
  expect_error(CD_path(node, dataSize, data_matrix, n_levels, obsIndex_R, eor_nr, eor, lambda_seq, nlam, eps, convLb, qtol, weights, gamma = -1, upperbound, threshold))
})

test_that("Check input: upperbound", {
  ### throw an error if upperbound is not a numeric number
  expect_error(CD_path(node, dataSize, data_matrix, n_levels, obsIndex_R, eor_nr, eor, lambda_seq, nlam, eps, convLb, qtol, weights, gamma, upperbound = as.character(100), threshold))

  ### throw an error if upperbound is negative and it is not -1
  expect_error(CD_path(node, dataSize, data_matrix, n_levels, obsIndex_R, eor_nr, eor, lambda_seq, nlam, eps, convLb, qtol, weights, gamma, upperbound = -100, threshold))

  ### if upper bound is -1, run as usual
  expect_error(CD_path(node, dataSize, data_matrix, n_levels, obsIndex_R, eor_nr, eor, lambda_seq, nlam, eps, convLb, qtol, weights, gamma, upperbound = -1, threshold), NA)
})

test_that("Check input: threshold", {
  ### throw an error if threshold is not an integer
  expect_error(CD_path(node, dataSize, data_matrix, n_levels, obsIndex_R, eor_nr, eor, lambda_seq, nlam, eps, convLb, qtol, weights, gamma, upperbound, threshold = 5))

  ### throw an error if threshold is negative
  expect_error(CD_path(node, dataSize, data_matrix, n_levels, obsIndex_R, eor_nr, eor, lambda_seq, nlam, eps, convLb, qtol, weights, gamma, upperbound, threshold = as.integer(-5)))
})

test_that("Check output", {
  ### output is a list of length 2
  fit <- CD_path(node, dataSize, data_matrix, n_levels, obsIndex_R, eor_nr, eor, lambda_seq, nlam, eps, convLb, qtol, weights, gamma, upperbound, threshold)
  expect_equal(length(fit), 3)
  expect_true("matrix" %in% class(fit[[1]]))
  expect_equal(class(fit[[2]]), "numeric")
  expect_true("matrix" %in% class(fit[[3]]))
})
