context("cd.run")

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
ivn_obs <- lapply(ivn, function(x){return(as.integer(0))})
ivn_int <- lapply(1:dataSize, function(x){return(as.integer(x/6))})
databn_obs <- sparsebnUtils::sparsebnData(data, ivn = ivn_obs, type = "discrete")
databn_int <- sparsebnUtils::sparsebnData(data, ivn = ivn_int, type = "discrete")
black_list = matrix(c("c", "a", "a", "b", "c", "d"), byrow = TRUE, ncol = 2)
white_list = matrix(c("a", "e", "a", "d"), byrow = TRUE, ncol = 2)
node_name <- colnames(data)
blocks <- lapply(node_name, function(x){
  # Allow all off-diagonal entries since we are no longer using the block decomposition
  row <- (node_name)[node_name != x]
  col <- rep(x, length(col))
  cbind(row, col)
})
blocks <- do.call("rbind", blocks)

# test
test_that("Testing default behaviour of cd.run", {
  final <- cd.run(databn_obs)
  # n_length <- length(final)
  # print(n_length)
  # print(final[[n_length]]$edges)
  # print(final[[n_length]]$lambda)
  # print(final[[n_length]]$nedge)

  ### check output type
  expect_is(final, "list")
  expect_is(final, "sparsebnPath")

  ### check element type of final
  for(i in seq_along(final)){
    expect_is(final[[i]], "sparsebnFit")
  }

  ### Check consistency of nedge
  for(i in seq_along(final)){
    matrix.nedge <- sum(sparsebnUtils::get.adjacency.matrix(final[[i]]$edges))
    edgeL.nedge <- sparsebnUtils::num.edges(final[[i]]$edges)
    expect_equal(final[[i]]$nedge, edgeL.nedge, matrix.nedge)
  }
})

test_that("Testing cd.run with manual settings", {
  weights <- matrix(1.5, nrow=node, ncol=node)
  final <- cd.run(databn_obs, weights=weights, lambdas=NULL, lambdas.length=10, error.tol=0.0003, convLb=0.02, weight.scale=1.5, upperbound=300.0, alpha = 5, permute = TRUE, adaptive = TRUE)

  ### check output type
  expect_is(final, "list")
  expect_is(final, "sparsebnPath")

  ### check element type of final
  for(i in seq_along(final)){
    expect_is(final[[i]], "sparsebnFit")
  }

  ### Check consistency of nedge
  for(i in seq_along(final)){
    matrix.nedge <- sum(sparsebnUtils::get.adjacency.matrix(final[[i]]$edges))
    edgeL.nedge <- sparsebnUtils::num.edges(final[[i]]$edges)
    expect_equal(final[[i]]$nedge, edgeL.nedge, matrix.nedge)
  }
})

test_that("Testing cd.run with intervention", {
  final <- cd.run(databn_int)
  # n_length <- length(final)
  # print(n_length)
  # print(final[[n_length]]$edges)
  # print(final[[n_length]]$lambda)
  # print(final[[n_length]]$nedge)

  ### check output type
  expect_is(final, "list")
  expect_is(final, "sparsebnPath")

  ### check element type of final
  for(i in seq_along(final)){
    expect_is(final[[i]], "sparsebnFit")
  }

  ### Check consistency of nedge
  for(i in seq_along(final)){
    matrix.nedge <- sum(sparsebnUtils::get.adjacency.matrix(final[[i]]$edges))
    edgeL.nedge <- sparsebnUtils::num.edges(final[[i]]$edges)
    expect_equal(final[[i]]$nedge, edgeL.nedge, matrix.nedge)
  }

  ### check result for black white list
  final <- cd.run(databn_obs, blacklist = black_list, whitelist = white_list)
  white_list_int <- matrix(match(white_list, colnames(data)), ncol = 2)
  black_list_int <- matrix(match(black_list, colnames(data)), ncol = 2)
  for(i in 1:length(final)){
    for(j in 1:nrow(white_list_int)){
      expect_true(white_list_int[j, 1] %in% final[[i]]$edges[[white_list_int[j, 2]]])
    }
    for(j in 1:nrow(black_list_int)){
      expect_true(!(black_list_int[j, 1] %in% final[[i]]$edges[[black_list_int[j, 2]]]))
    }
  }

  ### check result when all edges are in black list
  final_all_black <- cd.run(databn_obs, blacklist = blocks)
  expect_equal(length(final_all_black), 1)
  expect_equal(final_all_black[[1]]$nedge, 0)

  ### check result when all edges are in white list
  final_all_white <- cd.run(databn_obs, whitelist = blocks)
  n_edges <- lapply(final_all_white, function(x){x$nedge})
  expect_nedge <- node*(node+1)/2 - node
  expect_true(all(n_edges==expect_nedge))
})
