context("CD_call")

# set up input variable
data <- matrix(c(1, 1, 0, 0, 1, 1,
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
colnames(data) <- c("a", "b", "c", "d", "e", "f")
dataSize <- dim(data)[1]
node <- dim(data)[2]
ivn <- vector("list", length = dataSize)
ivn <- lapply(ivn, function(x){return(as.integer(0))})
databn <- sparsebnUtils::sparsebnData(data, ivn = ivn, type = "discrete")
black_list = matrix(c("c", "a", "a", "b", "c", "d"), byrow = TRUE, ncol = 2)
white_list = matrix(c("a", "e", "a", "d"), byrow = TRUE, ncol = 2)

# test
test_that("CD_call runs as expected", {
  ### no error
  expect_error(CD_call(databn, eor = NULL, permute = FALSE, weights = NULL, lambda_seq = NULL, fmlam = 0.1, nlam = 30, blacklist = NULL, whitelist = NULL, eps = 0.0001, convLb = 0.01, qtol = 0.0001, gamma = 1, upperbound = 100, threshold = 3), NA)
})

test_that("Check input for indata", {
  ### input with data.frame object
  data_frame_obj <- as.data.frame(data)
  expect_warning(CD_call(data_frame_obj, eor = NULL, permute = FALSE, weights = NULL, lambda_seq = NULL, fmlam = 0.1, nlam = 30, blacklist = NULL, whitelist = NULL, eps = 0.0001, convLb = 0.01, qtol = 0.0001, gamma = 1, upperbound = 100, threshold = 3))

  ### input wrong indata should throw error
  data_list <- as.list(data)
  expect_error(CD_call(indata = data_list, eor = NULL, permute = FALSE, weights = NULL, lambda_seq = NULL, fmlam = 0.1, nlam = 30, blacklist = NULL, whitelist = NULL, eps = 0.0001, convLb = 0.01, qtol = 0.0001, gamma = 1, upperbound = 100, threshold = 3))

  ### input wrong levels should throw an error
  data_oneLevel <- databn
  data_oneLevel$levels[[1]] <-0
  expect_error(CD_call(indata = data_oneLevel, eor = NULL, permute = FALSE, weights = NULL, lambda_seq = NULL, fmlam = 0.1, nlam = 30, blacklist = NULL, whitelist = NULL, eps = 0.0001, convLb = 0.01, qtol = 0.0001, gamma = 1, upperbound = 100, threshold = 3))
})

test_that("Check input ivn", {
  ### Input ivn
  ivn_list <- vector("list", length = dataSize)
  ivn_list <- lapply(1:length(ivn_list), function(x){return(as.integer(x%%7))})
  databn_ivn <- sparsebnUtils::sparsebnData(data, ivn = ivn_list, type = "discrete")
  expect_error(CD_call(databn_ivn, eor = NULL, permute = FALSE, weights = NULL, lambda_seq = NULL, fmlam = 0.1, nlam = 30, blacklist = NULL, whitelist = NULL, eps = 0.0001, convLb = 0.01, qtol = 0.0001, gamma = 1, upperbound = 100, threshold = 3), NA)

  ### Throw error when length of ivn list is wrong
  ivn_list <- vector("list", length = (dataSize+10))
  ivn_list <- lapply(1:length(ivn_list), function(x){return(as.integer(x%%7))})
  databn_ivn <- sparsebnUtils::sparsebnData(data, ivn = ivn_list, type = "discrete")
  expect_error(CD_call(databn_ivn, eor = NULL, permute = FALSE, weights = NULL, lambda_seq = NULL, fmlam = 0.1, nlam = 30, blacklist = NULL, whitelist = NULL, eps = 0.0001, convLb = 0.01, qtol = 0.0001, gamma = 1, upperbound = 100, threshold = 3))

  ### Works if ivn list contains NA
  ivn_list <- vector("list", length = dataSize)
  ivn_list <- lapply(1:length(ivn_list), function(x){return(ifelse(!x%%7, NA, as.integer(x%%7)))})
  databn_ivn <- sparsebnUtils::sparsebnData(data, ivn = ivn_list, type = "discrete")
  expect_error(CD_call(databn_ivn, eor = NULL, permute = FALSE, weights = NULL, lambda_seq = NULL, fmlam = 0.1, nlam = 30, blacklist = NULL, whitelist = NULL, eps = 0.0001, convLb = 0.01, qtol = 0.0001, gamma = 1, upperbound = 100, threshold = 3), NA)

  ### Throw error if ivn list indicate intervention on a node for every observation
  ivn_list <- vector("list", length = dataSize)
  ivn_list <- lapply(1:length(ivn_list), function(x){return(as.integer(c(x%%6, 6)))})
  databn_ivn <- sparsebnUtils::sparsebnData(data, ivn = ivn_list, type = "discrete")
  expect_error(CD_call(databn_ivn, eor = NULL, permute = FALSE, weights = NULL, lambda_seq = NULL, fmlam = 0.1, nlam = 30, blacklist = NULL, whitelist = NULL, eps = 0.0001, convLb = 0.01, qtol = 0.0001, gamma = 1, upperbound = 100, threshold = 3))

})

test_that("Check input weights", {
  ### Throw error if weight matrix is not a square matrix
  weights_notSquare <- matrix(1, node, node+1)
  expect_error(CD_call(databn, eor = NULL, permute = FALSE, weights = weights_notSquare, lambda_seq = NULL, fmlam = 0.1, nlam = 30, blacklist = NULL, whitelist = NULL, eps = 0.0001, convLb = 0.01, qtol = 0.0001, gamma = 1, upperbound = 100, threshold = 3))

  ### Throw error if weight matrix has wrong dimension
  weight_wrongDim <- matrix(1, node+1, node+1)
  expect_error(CD_call(databn, eor = NULL, permute = FALSE, weights = weights_wrongDim, lambda_seq = NULL, fmlam = 0.1, nlam = 30, blacklist = NULL, whitelist = NULL, eps = 0.0001, convLb = 0.01, qtol = 0.0001, gamma = 1, upperbound = 100, threshold = 3))

})

test_that("Check black white list", {
  ### throw error if black white list have common edge
  black_mixed_white <- matrix(c("a", "e", "c", "a", "a", "b", "c", "d"), byrow = TRUE, ncol = 2)
  expect_error(CD_call(databn, eor = NULL, permute = FALSE, weights = NULL, lambda_seq = NULL, fmlam = 0.1, nlam = 30, blacklist = black_mixed_white, whitelist = white_list, eps = 0.0001, convLb = 0.01, qtol = 0.0001, gamma = 1, upperbound = 100, threshold = 3))

  ### throw error if black white list is not a matrix
  expect_error(CD_call(databn, eor = NULL, permute = FALSE, weights = NULL, lambda_seq = NULL, fmlam = 0.1, nlam = 30, blacklist = as.list(black_list), whitelist = white_list, eps = 0.0001, convLb = 0.01, qtol = 0.0001, gamma = 1, upperbound = 100, threshold = 3))
  expect_error(CD_call(databn, eor = NULL, permute = FALSE, weights = NULL, lambda_seq = NULL, fmlam = 0.1, nlam = 30, blacklist = black_list, whitelist = as.list(white_list), eps = 0.0001, convLb = 0.01, qtol = 0.0001, gamma = 1, upperbound = 100, threshold = 3))

  ### throw error if black white list has wrong dimension
  expect_error(CD_call(databn, eor = NULL, permute = FALSE, weights = NULL, lambda_seq = NULL, fmlam = 0.1, nlam = 30, blacklist = matrix(black_list, ncol = 3), whitelist = white_list, eps = 0.0001, convLb = 0.01, qtol = 0.0001, gamma = 1, upperbound = 100, threshold = 3))
  expect_error(CD_call(databn, eor = NULL, permute = FALSE, weights = NULL, lambda_seq = NULL, fmlam = 0.1, nlam = 30, blacklist = black_list, whitelist = matrix(white_list, ncol = 4), eps = 0.0001, convLb = 0.01, qtol = 0.0001, gamma = 1, upperbound = 100, threshold = 3))

  ### throw error if black white list has na
  black_list_na <- black_list
  black_list_na[1, 1] <- NA
  white_list_na <- white_list
  white_list_na[1, 1] <- NA
  expect_error(CD_call(databn, eor = NULL, permute = FALSE, weights = NULL, lambda_seq = NULL, fmlam = 0.1, nlam = 30, blacklist = black_list_na, whitelist = white_list, eps = 0.0001, convLb = 0.01, qtol = 0.0001, gamma = 1, upperbound = 100, threshold = 3))
  expect_error(CD_call(databn, eor = NULL, permute = FALSE, weights = NULL, lambda_seq = NULL, fmlam = 0.1, nlam = 30, blacklist = black_list, whitelist = white_list_na, eps = 0.0001, convLb = 0.01, qtol = 0.0001, gamma = 1, upperbound = 100, threshold = 3))
})

