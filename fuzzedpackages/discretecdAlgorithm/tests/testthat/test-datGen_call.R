context("datGen_call")

# set up input variables
edge_list <- vector("list", length = 5)
edge_list[[1]] <- integer(0)
edge_list[[2]] <- integer(0)
edge_list[[3]] <- 1
edge_list[[4]] <- c(1, 3)
edge_list[[5]] <- c(1, 2)
edge_list <- sparsebnUtils::as.edgeList(edge_list)
names(edge_list) <- c("V1", "V2", "V3", "V4", "V5")
nlevels <- c(3, 5, 2, 2, 3)
data_params <- coef_gen(edge_list = edge_list, n_levels = nlevels)
ivn_list <- rep(1:5, rep(10, 5))
ivn_list <- as.list(ivn_list)
ivn_list <- lapply(ivn_list, function(x){paste0("V", x)})
n_obs <- length(ivn_list)

# test
test_that("datGen_call runs as expected", {
  expect_error(datGen_call(edge_list = edge_list, dataSize = n_obs, ivn = ivn_list, ivn_rand = TRUE, nlevels = nlevels, coef = data_params), NA)

})

test_that("Check input edge_list", {
  ### Throw an error if edge_list is of the wrong type
  edge_list_list <- as.list(edge_list)
  expect_error(datGen_call(edge_list = edge_list_list, dataSize = n_obs, ivn = ivn_list, ivn_rand = TRUE, nlevels = nlevels, coef = data_params))

  ### Throw an error if the edge_list is not named
  edge_list_noName <- edge_list
  names(edge_list_noName) <- NULL
  expect_error(datGen_call(edge_list = edge_list_noName, dataSize = n_obs, ivn = ivn_list, ivn_rand = TRUE, nlevels = nlevels, coef = data_params))

})

test_that("Check input dataSize", {
  ### Throw error if dataSize is not a scaler
  n_obs_vector <- c(n_obs, 1)
  expect_error(datGen_call(edge_list = edge_list, dataSize = n_obs_vector, ivn = ivn_list, ivn_rand = TRUE, nlevels = nlevels, coef = data_params))

  ### Throw error if dataSize is not numeric
  n_obs_noNumeric <- character(50)
  expect_error(datGen_call(edge_list = edge_list, dataSize = n_obs_noNumeric, ivn = ivn_list, ivn_rand = TRUE, nlevels = nlevels, coef = data_params))

  ### Throw error if dataSize is not positive
  n_obs_neg <- -1
  expect_error(datGen_call(edge_list = edge_list, dataSize = n_obs_neg, ivn = ivn_list, ivn_rand = TRUE, nlevels = nlevels, coef = data_params))

})

test_that("Check input ivn", {
  ### Run as expected if ivn is NULL
  expect_error(datGen_call(edge_list = edge_list, dataSize = n_obs, ivn = NULL, ivn_rand = TRUE, nlevels = nlevels, coef = data_params), NA)

  ### Throw error if ivn is not a list
  ivn_vector <- rep(1:5, rep(10, 5))
  expect_error(datGen_call(edge_list = edge_list, dataSize = n_obs, ivn = ivn_vector, ivn_rand = TRUE, nlevels = nlevels, coef = data_params))

  ### Throw error if length of ivn is not compatiable with dataSize
  ivn_wrongSize <- c(rep(1:5, rep(10, 5)), rep(0, 10))
  ivn_wrongSize <- as.list(ivn_wrongSize)
  names(ivn_wrongSize) <- lapply(ivn_wrongSize, function(x){paste0("V", x)})
  expect_error(datGen_call(edge_list = edge_list, dataSize = n_obs, ivn = ivn_wrongSize, ivn_rand = TRUE, nlevels = nlevels, coef = data_params))

})

test_that("Check input nlevels", {
  ### Throw error if nlevels is not a vector
  nlevels_matrix <- matrix(0, 50, 50)
  expect_error(datGen_call(edge_list = edge_list, dataSize = n_obs, ivn = ivn_list, ivn_rand = TRUE, nlevels = nlevels_list, coef = data_params))

  ### Throw error if element of nlevels is smaller than 2
  nlevels_fewLevels <- c(1, 5, 2, 2, 3)
  expect_error(datGen_call(edge_list = edge_list, dataSize = n_obs, ivn = ivn_list, ivn_rand = TRUE, nlevels = nlevels_fewLevels, coef = data_params))

  ### Throw error if length of nlevels is not number of nodes
  nlevels_wrongLength <- c(3, 5, 2, 2)
  expect_error(datGen_call(edge_list = edge_list, dataSize = n_obs, ivn = ivn_list, ivn_rand = TRUE, nlevels = nlevels_wrongLength, coef = data_params))

})

test_that("Check input coef", {
  ### Throw error if coef is not a list
  coef_vector <- rnorm(50)
  expect_error(datGen_call(edge_list = edge_list, dataSize = n_obs, ivn = ivn_list, ivn_rand = TRUE, nlevels = nlevels, coef = coef_vector))

  ### Throw error if coef has wrong length
  coef_wrongDim <- data_params[1:4]
  expect_error(datGen_call(edge_list = edge_list, dataSize = n_obs, ivn = ivn_list, ivn_rand = TRUE, nlevels = nlevels, coef = coef_wrongDim))

  ### Throw error if coef is not a list of matrix or NA
  coef_wrongCont <- data_params
  coef_wrongCont[[4]] <- as.vector(coef_wrongCont[[4]])
  expect_error(datGen_call(edge_list = edge_list, dataSize = n_obs, ivn = ivn_list, ivn_rand = TRUE, nlevels = nlevels, coef = coef_wrongCont))

  ### Throw error if the coefficient not compatable with nlevels
  coef_wrongLevel <- data_params
  coef_wrongLevel[[3]] <- cbind(coef_wrongLevel[[3]], rnorm(2))
  expect_error(datGen_call(edge_list = edge_list, dataSize = n_obs, ivn = ivn_list, ivn_rand = TRUE, nlevels = nlevels, coef = coef_wrongLevel))

})
