context("get.summary")

# set up input of function get.summary
adj_matrix <- matrix(c(0, 1, 1, 0,
                       0, 0, 1, 0,
                       0, 0, 0, 1,
                       0, 0, 0, 0), byrow = TRUE, ncol = 4)
dataSize <- 50
lambda <- 3.5
time <- 1

# test get.summary
test_that("input adjacency matrix can be converted to the right edge list.", {
  ### Trivial case
  adj_matrix_trivial <- matrix(rep(0, 16), nrow = 4)
  no_parent <- vector("list", length = 4)
  no_parent <- lapply(no_parent, function(x){return(integer(0))})
  expect_equal(get.summary(adj_matrix_trivial, dataSize, lambda, time)$nedge, 0)
  expect_equal(get.summary(adj_matrix_trivial, dataSize, lambda, time)$edges, sparsebnUtils::edgeList(no_parent))

  ### when all nodes has just one parent
  adj_matrix_oneparent <- matrix(c(0, 0, 0, 0,
                         0, 0, 0, 0,
                         0, 0, 0, 0,
                         0, 1, 1, 1), byrow = TRUE, ncol = 4)
  expect_equal(get.summary(adj_matrix_oneparent, dataSize, lambda, time)$nedge, 3)
  true_edgeList_one <- vector("list", length = 4)
  true_edgeList_one[[1]] <- integer(0)
  true_edgeList_one[[2]] <- 4L
  true_edgeList_one[[3]] <- 4L
  true_edgeList_one[[4]] <- 4L
  true_edgeList_one <- sparsebnUtils::edgeList(true_edgeList_one)
  expect_equal(get.summary(adj_matrix_oneparent, dataSize, lambda, time)$edges, true_edgeList_one)

  ### Non-trivil case
  expect_equal(get.summary(adj_matrix, dataSize, lambda, time)$nedge, 4)
  true_edgeList <- vector("list", length = 4)
  true_edgeList[[1]] <- integer(0)
  true_edgeList[[2]] <- as.integer(1)
  true_edgeList[[3]] <- as.integer(c(1, 2))
  true_edgeList[[4]] <- as.integer(3)
  true_edgeList <- sparsebnUtils::edgeList(true_edgeList)
  expect_equal(get.summary(adj_matrix, dataSize, lambda, time)$edges, true_edgeList)
})

test_that("output the right lambda.", {
  expect_equal(get.summary(adj_matrix, dataSize, lambda, time)$lambda, 3.5)
})

test_that("output the right time.", {
  ### Trivial case
  expect_equal(get.summary(adj_matrix, dataSize, lambda, time=NA)$time, NA)

  ### Non-trivial case
  expect_equal(get.summary(adj_matrix, dataSize, lambda, time)$time, 1)
})

test_that("output the rigth nn.", {
  expect_equal(get.summary(adj_matrix, dataSize, lambda, time)$nn, 50)
})

test_that("output the right pp.", {
  expect_equal(get.summary(adj_matrix, dataSize, lambda, time)$pp, 4)
})
