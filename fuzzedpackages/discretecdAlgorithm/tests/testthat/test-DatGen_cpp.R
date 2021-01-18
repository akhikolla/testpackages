context("DatGen_cpp")

# set up input variables
edge_list <- vector("list", length = 5)
edge_list[[1]] <- integer(0)
edge_list[[2]] <- integer(0)
edge_list[[3]] <- 1
edge_list[[4]] <- c(1, 3)
edge_list[[5]] <- c(1, 2)
names(edge_list) <- c("V1", "V2", "V3", "V4", "V5")
nlevels <- c(3, 5, 2, 2, 3)
coef <- coef_gen(edge_list = sparsebnUtils::as.edgeList(edge_list), n_levels = nlevels)
ivn_list <- rep(1:5, rep(10, 5))
ivn_list <- as.list(ivn_list)
ivn_list <- lapply(ivn_list, function(x){paste0("V", x)})
dataSize <- length(ivn_list)

maxdeg <- max(sapply(edge_list, length))
maxdeg <- as.integer(maxdeg)
node <- length(edge_list)
node<- as.integer(node)
ordex <- matrix(c(0, 0, 1, 1, 1, 0, 0, 0, 3, 2), ncol = 5, byrow = TRUE)
ordex <- as.integer(ordex)
ordex <- matrix(ordex, nrow = 2)
ts <- as.integer(1:5)
dataSize <- as.integer(dataSize)
ivn <- as.list(rep(1:5, rep(10, 5)))
ivn_vals <- lapply(ivn, as.integer)
ivn <- lapply(ivn, function(x){
  as.integer(x-1)})
nlevels <- as.integer(nlevels)
coef_list <- lapply(1:node, function(x, coef, edge_list, nlevels){
  out <- NULL
  if (length(edge_list[[x]])) {
    index <- rep(1:(length(edge_list[[x]])+1), c(1, nlevels[edge_list[[x]]]-1))
    m_to_list <- vector("list", length = length(edge_list[[x]])+1)
    for (i in 1:(length(edge_list[[x]])+1)) {
      m_to_list[[i]] <- coef[[x]][, index==i, drop = FALSE]
    }
    out = m_to_list
  }
  out
}, coef, edge_list, nlevels)
coef_length <- sapply(coef_list, length)
coef_length <- as.integer(coef_length)

test_that("DatGen_cpp run as expected", {
  expect_error(DatGen_cpp(maxdeg = maxdeg, node = node, ordex = ordex, ts = ts, dataSize = dataSize, ivn = ivn, ivn_vals = ivn_vals, ivn_rand = TRUE, nlevels = nlevels, coef_list = coef_list, coef_length = coef_length), NA)
})

test_that("Check input maxdeg", {
  ### Throw error if maxdeg is not an integer
  maxdeg_num <- 2.0
  expect_error(DatGen_cpp(maxdeg = maxdeg_num, node = node, ordex = ordex, ts = ts, dataSize = dataSize, ivn = ivn, ivn_vals = ivn_vals, ivn_rand = TRUE, nlevels = nlevels, coef_list = coef_list, coef_length = coef_length))

  ### Throw error if maxdeg is not positive
  maxdeg_neg <- as.integer(-1)
  expect_error(DatGen_cpp(maxdeg = maxdeg_neg, node = node, ordex = ordex, ts = ts, dataSize = dataSize, ivn = ivn, ivn_vals = ivn_vals, ivn_rand = TRUE, nlevels = nlevels, coef_list = coef_list, coef_length = coef_length))

})

test_that("Check input node", {
  ### Throw error if node is not a integer
  node_num <- 5.0
  expect_error(DatGen_cpp(maxdeg = maxdeg, node = node_num, ordex = ordex, ts = ts, dataSize = dataSize, ivn = ivn, ivn_vals = ivn_vals, ivn_rand = TRUE, nlevels = nlevels, coef_list = coef_list, coef_length = coef_length))

  ### Throw error if node is not positive
  node_neg <- as.integer(-1)
  expect_error(DatGen_cpp(maxdeg = maxdeg, node = node_neg, ordex = ordex, ts = ts, dataSize = dataSize, ivn = ivn, ivn_vals = ivn_vals, ivn_rand = TRUE, nlevels = nlevels, coef_list = coef_list, coef_length = coef_length))

})

test_that("Check input ordex", {
  ### Throw error if ordex is not a matrix
  ordex_vec <- ordex[1, ]
  expect_error(DatGen_cpp(maxdeg = maxdeg, node = node, ordex = ordex_vec, ts = ts, dataSize = dataSize, ivn = ivn, ivn_vals = ivn_vals, ivn_rand = TRUE, nlevels = nlevels, coef_list = coef_list, coef_length = coef_length))

  ### Throw error if ordex is not a matrix with integer entries
  ordex_num <- as.numeric(ordex)
  ordex_num <- matrix(ordex_num, nrow = maxdeg)
  expect_error(DatGen_cpp(maxdeg = maxdeg, node = node, ordex = ordex_num, ts = ts, dataSize = dataSize, ivn = ivn, ivn_vals = ivn_vals, ivn_rand = TRUE, nlevels = nlevels, coef_list = coef_list, coef_length = coef_length))

  ### Throw error if dimension of ordex is not right
  ordex_wrongRow <- rbind(ordex, rep(0, 5))
  ordex_wrongCol <- ordex[, 1:4]
  expect_error(DatGen_cpp(maxdeg = maxdeg, node = node, ordex = ordex_wrongRow, ts = ts, dataSize = dataSize, ivn = ivn, ivn_vals = ivn_vals, ivn_rand = TRUE, nlevels = nlevels, coef_list = coef_list, coef_length = coef_length))
  expect_error(DatGen_cpp(maxdeg = maxdeg, node = node, ordex = ordex_wrongCol, ts = ts, dataSize = dataSize, ivn = ivn, ivn_vals = ivn_vals, ivn_rand = TRUE, nlevels = nlevels, coef_list = coef_list, coef_length = coef_length))

})

test_that("Check input ts", {
  ### Throw error if ts is not a vector
  ts_matrix <- matrix(1, 5, 5)
  expect_error(DatGen_cpp(maxdeg = maxdeg, node = node, ordex = ordex, ts = ts_matrix, dataSize = dataSize, ivn = ivn, ivn_vals = ivn_vals, ivn_rand = TRUE, nlevels = nlevels, coef_list = coef_list, coef_length = coef_length))

  ### Throw error if ts has wrong length
  ts_short <- as.integer(1:4)
  expect_error(DatGen_cpp(maxdeg = maxdeg, node = node, ordex = ordex, ts = ts_short, dataSize = dataSize, ivn = ivn, ivn_vals = ivn_vals, ivn_rand = TRUE, nlevels = nlevels, coef_list = coef_list, coef_length = coef_length))

  ### Throw error if ts is not an integer vector
  ts_num <- seq(1, 5, 1)
  expect_error(DatGen_cpp(maxdeg = maxdeg, node = node, ordex = ordex, ts = ts_num, dataSize = dataSize, ivn = ivn, ivn_vals = ivn_vals, ivn_rand = TRUE, nlevels = nlevels, coef_list = coef_list, coef_length = coef_length))

  ### Throw error if ts is not a valid ordering
  ts_notOrder <- as.integer(1, 2, 3, 3, 5)
  expect_error(DatGen_cpp(maxdeg = maxdeg, node = node, ordex = ordex, ts = ts_notOrder, dataSize = dataSize, ivn = ivn, ivn_vals = ivn_vals, ivn_rand = TRUE, nlevels = nlevels, coef_list = coef_list, coef_length = coef_length))

})

test_that("Check input dataSize", {
  ### Throw error if dataSize is not an integer
  dataSize_num <- as.numeric(dataSize)
  expect_error(DatGen_cpp(maxdeg = maxdeg, node = node, ordex = ordex, ts = ts, dataSize = dataSize_num, ivn = ivn, ivn_vals = ivn_vals, ivn_rand = TRUE, nlevels = nlevels, coef_list = coef_list, coef_length = coef_length))

  ### Throw error if dataSize is not an integer
  dataSize_neg <- 0L
  expect_error(DatGen_cpp(maxdeg = maxdeg, node = node, ordex = ordex, ts = ts, dataSize = dataSize_neg, ivn = ivn, ivn_vals = ivn_vals, ivn_rand = TRUE, nlevels = nlevels, coef_list = coef_list, coef_length = coef_length))

})

test_that("Check input ivn", {
  ### Throw error if ivn is not a list
  ivn_vec <- unlist(ivn)
  expect_error(DatGen_cpp(maxdeg = maxdeg, node = node, ordex = ordex, ts = ts, dataSize = dataSize, ivn = ivn_vec, ivn_vals = ivn_vals, ivn_rand = TRUE, nlevels = nlevels, coef_list = coef_list, coef_length = coef_length))

  ### Throw error if ivn is not a list of integer vectors
  ivn_num <- ivn
  ivn_num[[1]] <- as.numeric(ivn_num)
  expect_error(DatGen_cpp(maxdeg = maxdeg, node = node, ordex = ordex, ts = ts, dataSize = dataSize, ivn = ivn_num, ivn_vals = ivn_vals, ivn_rand = TRUE, nlevels = nlevels, coef_list = coef_list, coef_length = coef_length))

  ### Throw error if ivn has wrong length
  ivn_short <- ivn[1:10]
  expect_error(DatGen_cpp(maxdeg = maxdeg, node = node, ordex = ordex, ts = ts, dataSize = dataSize, ivn = ivn_short, ivn_vals = ivn_vals, ivn_rand = TRUE, nlevels = nlevels, coef_list = coef_list, coef_length = coef_length))

})

test_that("Check input nlevels", {
  ### Throw error if nlevels is not a vector
  nlevels_matrix <- matrix(1, 1, 5)
  expect_error(DatGen_cpp(maxdeg = maxdeg, node = node, ordex = ordex, ts = ts, dataSize = dataSize, ivn = ivn, ivn_vals = ivn_vals, ivn_rand = TRUE, nlevels = nlevels_matrix, coef_list = coef_list, coef_length = coef_length))

  ### Throw error if nlevels is not a vector of integers
  nlevels_num <- as.numeric(nlevels)
  expect_error(DatGen_cpp(maxdeg = maxdeg, node = node, ordex = ordex, ts = ts, dataSize = dataSize, ivn = ivn, ivn_vals = ivn_vals, ivn_rand = TRUE, nlevels = nlevels_num, coef_list = coef_list, coef_length = coef_length))

  ### Throw error if length of nlevels is wrong
  nlevels_short <- nlevels[1:4]
  expect_error(DatGen_cpp(maxdeg = maxdeg, node = node, ordex = ordex, ts = ts, dataSize = dataSize, ivn = ivn, ivn_vals = ivn_vals, ivn_rand = TRUE, nlevels = nlevels_short, coef_list = coef_list, coef_length = coef_length))

})

test_that("Check input coef", {
  ### Throw error if coef is not a list
  coef_vec <- unlist(coef_list)
  expect_error(DatGen_cpp(maxdeg = maxdeg, node = node, ordex = ordex, ts = ts, dataSize = dataSize, ivn = ivn, ivn_vals = ivn_vals, ivn_rand = TRUE, nlevels = nlevels, coef_list = coef_vec, coef_length = coef_length))

  ### Throw error if coef is not a list of list
  expect_error(DatGen_cpp(maxdeg = maxdeg, node = node, ordex = ordex, ts = ts, dataSize = dataSize, ivn = ivn, ivn_vals = ivn_vals, ivn_rand = TRUE, nlevels = nlevels, coef_list = coef, coef_length = coef_length))

  ### Throw error if element of list of coef_list is not a matrix or NA
  coef_list_vec <- coef_list
  coef_list_vec[[4]][[3]] <- as.vector(coef_list_vec[[4]][[3]])
  expect_error(DatGen_cpp(maxdeg = maxdeg, node = node, ordex = ordex, ts = ts, dataSize = dataSize, ivn = ivn, ivn_vals = ivn_vals, ivn_rand = TRUE, nlevels = nlevels, coef_list = coef_list_vec, coef_length = coef_length))

  ### Throw error if length of coef_list is not right
  coef_list_wrongSize <- coef_list[1:4]
  expect_error(DatGen_cpp(maxdeg = maxdeg, node = node, ordex = ordex, ts = ts, dataSize = dataSize, ivn = ivn, ivn_vals = ivn_vals, ivn_rand = TRUE, nlevels = nlevels, coef_list = coef_list_wrongSize, coef_length = coef_length))

  ### Throw error if length of sublist of coef_list is not right
  coef_list_sub <- coef_list
  coef_list_sub[[4]] <- coef_list_sub[[4]][1:2]
  expect_error(DatGen_cpp(maxdeg = maxdeg, node = node, ordex = ordex, ts = ts, dataSize = dataSize, ivn = ivn, ivn_vals = ivn_vals, ivn_rand = TRUE, nlevels = nlevels, coef_list = coef_list_sub, coef_length = coef_length))

  ### Throw error if dim of intercept is not right
  coef_list_inter_row <- coef_list
  coef_list_inter_col <- coef_list
  coef_list_inter_row[[3]][[1]] <- matrix(1, nrow = 3, ncol = 1)
  coef_list_inter_col[[3]][[1]] <- matrix(1, nrow = 2, ncol = 2)
  expect_error(DatGen_cpp(maxdeg = maxdeg, node = node, ordex = ordex, ts = ts, dataSize = dataSize, ivn = ivn, ivn_vals = ivn_vals, ivn_rand = TRUE, nlevels = nlevels, coef_list = coef_list_inter_row, coef_length = coef_length))
  expect_error(DatGen_cpp(maxdeg = maxdeg, node = node, ordex = ordex, ts = ts, dataSize = dataSize, ivn = ivn, ivn_vals = ivn_vals, ivn_rand = TRUE, nlevels = nlevels, coef_list = coef_list_inter_col, coef_length = coef_length))

  ###Throw error if dim of coefficient matrix is not right
  coef_list_matrix_row <- coef_list
  coef_list_matrix_col <- coef_list
  coef_list_matrix_row[[5]][[3]] <- coef_list_matrix_row[[5]][[3]][1:2, ]
  coef_list_matrix_col[[5]][[3]] <- coef_list_matrix_col[[5]][[3]][, 1:3]
  expect_error(DatGen_cpp(maxdeg = maxdeg, node = node, ordex = ordex, ts = ts, dataSize = dataSize, ivn = ivn, ivn_vals = ivn_vals, ivn_rand = TRUE, nlevels = nlevels, coef_list = coef_list_matrix_row, coef_length = coef_length))
  expect_error(DatGen_cpp(maxdeg = maxdeg, node = node, ordex = ordex, ts = ts, dataSize = dataSize, ivn = ivn, ivn_vals = ivn_vals, ivn_rand = TRUE, nlevels = nlevels, coef_list = coef_list_matrix_col, coef_length = coef_length))

})

test_that("Check input coef_length", {
  ### Throw error if coeff_length is not integer vector
  coef_length_num <- as.numeric(coef_length)
  expect_error(DatGen_cpp(maxdeg = maxdeg, node = node, ordex = ordex, ts = ts, dataSize = dataSize, ivn = ivn, ivn_vals = ivn_vals, ivn_rand = TRUE, nlevels = nlevels, coef_list = coef_list, coef_length = coef_length_num))

})
