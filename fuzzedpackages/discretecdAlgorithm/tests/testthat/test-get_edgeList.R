context("get.edgeList")

# set up input of function get.edgeList
adj_matrix <- matrix(c(0, 1, 1, 0,
                       0, 0, 1, 0,
                       0, 0, 0, 1,
                       0, 0, 0, 0,
                       0, 0, 1, 1,
                       0, 0, 1, 1,
                       0, 0, 0, 1,
                       0, 0, 0, 0), byrow = TRUE, ncol = 4)
dataSize <- 50
lambda <- 3.5
time <- 1

# test get.edgeList
test_that("get.edgeList function well.", {
  ### Trivial case
  adj_matrix_trivial <- matrix(rep(0, 32), ncol = 4)

  expect_equal(length(get.edgeList(adj_matrix_trivial, dataSize, lambda, time)), 2)

  expect_equal(get.edgeList(adj_matrix_trivial, dataSize, lambda, time)[[1]]$nedge, 0)

  # no parent case
  no_parent <- vector("list", length = 4)
  no_parent <- lapply(no_parent, function(x){return(integer(0))})
  expect_equal(get.edgeList(adj_matrix_trivial, dataSize, lambda, time)[[1]]$edges, sparsebnUtils::edgeList(no_parent))

  expect_equal(get.edgeList(adj_matrix_trivial, dataSize, lambda, time)[[2]]$nedge, 0)
  expect_equal(get.edgeList(adj_matrix_trivial, dataSize, lambda, time)[[2]]$edges, sparsebnUtils::edgeList(no_parent))

  ### Non-trivial case
  expect_equal(length(get.edgeList(adj_matrix, dataSize, lambda, time)), 2)

  expect_equal(get.edgeList(adj_matrix, dataSize, lambda, time)[[1]]$nedge, 4)
  true_edgeList <- vector("list", length = 4)
  true_edgeList[[1]] <- which(adj_matrix == 2)
  true_edgeList[[2]] <- as.integer(1)
  true_edgeList[[3]] <- as.integer(c(1, 2))
  true_edgeList[[4]] <- as.integer(3)
  true_edgeList <- sparsebnUtils::edgeList(true_edgeList)
  expect_equal(get.edgeList(adj_matrix, dataSize, lambda, time)[[1]]$edges, true_edgeList)

  expect_equal(get.edgeList(adj_matrix, dataSize, lambda, time)[[2]]$nedge, 5)
  true_edgeList <- vector("list", length = 4)
  true_edgeList[[1]] <- which(adj_matrix == 2)
  true_edgeList[[2]] <- which(adj_matrix == 2)
  true_edgeList[[3]] <- as.integer(c(1, 2))
  true_edgeList[[4]] <- as.integer(c(1, 2, 3))
  true_edgeList <- sparsebnUtils::edgeList(true_edgeList)
  expect_equal(get.edgeList(adj_matrix, dataSize, lambda, time)[[2]]$edges, true_edgeList)
})
