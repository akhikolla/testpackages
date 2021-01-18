context("data_gen")

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
test_that("data_gen run as expected", {
  ### data_gen run with default settings
  expect_error(data_gen(graph = edge_list, n = n_obs), NA)

  ### data_gen run with manual settings
  expect_error(data_gen(graph = edge_list, n = n_obs, ivn = ivn_list, n_levels = nlevels, params = data_params, FUN = unif, flip = TRUE), NA)
})

test_that("check input graph", {
  graph_list <- as.list(edge_list)
  expect_error(data_gen(graph = graph_list, n = n_obs))

})

test_that("check input n", {
  n_neg <- -1
  expect_error(data_gen(graph = edge_list, n = n_neg))

})

test_that("check input n_levels", {
  n_levels_wrong <- nlevels
  n_levels_wrong[1] <- 1
  expect_error(data_gen(graph = edge_list, n = n_obs, n_levels = n_levels_wrong))

})

test_that("Check output", {
  generated_data <- data_gen(graph = edge_list, n = n_obs)
  expect_true("matrix" %in% class(generated_data))
  expect_equal(nrow(generated_data), n_obs)
  expect_equal(ncol(generated_data), length(edge_list))
})
