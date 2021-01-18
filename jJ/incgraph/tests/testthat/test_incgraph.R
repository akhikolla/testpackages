context("Testing IncGraph")

test_that("Simple test", {
  links <- tibble::tribble(
    ~from, ~to,
    1, 2,
    1, 3,
    1, 4,
    1, 5,
    1, 6,
    1, 7,
    2, 7,
    2, 8,
    2, 9,
    2, 10
  ) %>% as.matrix

  net <- new.incgraph.network(links = links)

  # Calculate the initial orbit counts using orca
  orb.counts <- calculate.orbit.counts(net)
  expect_equal(nrow(orb.counts), 10)
  expect_equal(ncol(orb.counts), 73)


  # add (5,10) and test
  flip(net, 5, 10)
  delta1 <- calculate.delta(net, 5, 10)

  expect_equal(nrow(delta1$add), 10)
  expect_equal(ncol(delta1$add), 73)
  expect_equal(nrow(delta1$rem), 10)
  expect_equal(ncol(delta1$rem), 73)

  # Verify that the new orbit counts equals the old orbit counts plus the delta counts
  new.inc.counts <- orb.counts + delta1$add - delta1$rem
  new.orb.counts <- calculate.orbit.counts(net)

  expect_equal(new.inc.counts, new.orb.counts)



  # Modify another edge
  flip(net, 6, 10) # add (6, 10)
  delta2 <- calculate.delta(net, 6, 10)

  # Verify that the new orbit counts equals the old orbit counts plus the delta counts
  new.inc.counts <- orb.counts +
    delta1$add - delta1$rem +
    delta2$add - delta2$rem
  new.orb.counts <- calculate.orbit.counts(net)

  expect_equal(new.inc.counts, new.orb.counts)



  # And another
  flip(net, 1, 5)  # remove (1, 5)
  delta3 <- calculate.delta(net, 1, 5)

  # Verify that the new orbit counts equals the old orbit counts plus the delta counts
  new.inc.counts <- orb.counts +
    delta1$add - delta1$rem +
    delta2$add - delta2$rem +
    delta3$add - delta3$rem
  new.orb.counts <- calculate.orbit.counts(net)

  expect_equal(new.inc.counts, new.orb.counts)

  ## Test Additional helper functions
  # Transform the network to a matrix
  newmat <- network.as.matrix(net)
  expect_equal(nrow(newmat), 11)
  expect_equal(ncol(newmat), 2)

  should_be_links <- tibble::tribble(
    ~from, ~to,
    1, 2,
    1, 3,
    1, 4,
    1, 6,
    1, 7,
    2, 7,
    2, 8,
    2, 9,
    2, 10,
    5, 10,
    6, 10
  ) %>% as.matrix
  pa.new <- c(
    paste0(newmat[,1], "|", newmat[,2]),
    paste0(newmat[,2], "|", newmat[,1])
  ) %>% sort
  pa.orig <- c(
    paste0(should_be_links[,1], "|", should_be_links[,2]),
    paste0(should_be_links[,2], "|", should_be_links[,1])
  ) %>% sort

  expect_equal(pa.new, pa.orig)

  # Get all neighbours of a node
  for (i in seq_len(10)) {
    neighs <- get.neighbours(net, i) %>% sort
    should_be_neighs <- c(
      should_be_links[should_be_links[,1] == i, 2],
      should_be_links[should_be_links[,2] == i, 1]
    ) %>% sort %>% setNames(NULL)
    expect_equal(neighs, should_be_neighs)
  }
  # Does the network contain a specific interaction?
  expect_true(contains(net, 5, 10))
  expect_false(contains(net, 7, 10))

  # Reinitialise to an empty network
  reset(net)

  new_mat <- network.as.matrix(net)
  expect_equal(nrow(new_mat), 0)
})

