context("get_obsIndex")

# set up input of function get_obsIndex
node <- 5
# ivn for observational data
ivn_obs <- list(0, 0, 0, 0, 0, 0)
# ivn, when some nodes are under intervention in all observations
ivn_node <- list(c(1, 2, 3), c(1, 2), c(1, 3), c(1, 4), c(1, 6), 1)
# ivn, when in some of the observation, all nodes are under intervention.
ivn_all <- list(c(1, 2, 3, 4, 5), c(2, 4), c(1, 3), c(1, 2, 3, 4, 5), 0, 0)
# normal ivn
ivn <- list(1, 2, 3, 4, 5, 2)

# test get_obsIndex
test_that("for observational data, get_obsIndex output a right outcome.", {
  true_ivn <- list(c(1, 2, 3, 4, 5, 6), c(1, 2, 3, 4, 5, 6), c(1, 2, 3, 4, 5, 6),
                   c(1, 2, 3, 4, 5, 6), c(1, 2, 3, 4, 5, 6))
  expect_equal(get_obsIndex(ivn_obs, node), true_ivn)
})

test_that("when some nodes are under intervention in all observations, get_obsIndex output a right outcome", {
  true_ivn <- list(0L, c(3, 4, 5, 6), c(2, 4, 5, 6), c(1, 2, 3, 5, 6), c(1, 2, 3, 4, 5, 6))
  expect_equal(get_obsIndex(ivn_node, node), true_ivn)
})

test_that("when in some observation, all nodes are under intervention, get_obsIndex output a right outcome", {
  true_ivn <- list(c(2, 5, 6), c(3, 5, 6), c(2, 5, 6), c(3, 5, 6), c(2, 3, 5, 6))
  expect_equal(get_obsIndex(ivn_all, node), true_ivn)
})

test_that("check in non-trivial case, if get_obsIndex output a right outcome.", {
  true_ivn <- list(c(2, 3, 4, 5, 6), c(1, 3, 4, 5), c(1, 2, 4, 5, 6), c(1, 2, 3, 5, 6), c(1, 2, 3, 4, 6))
  expect_equal(get_obsIndex(ivn, node), true_ivn)
})
