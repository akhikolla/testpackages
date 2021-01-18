# test-GridOnClusters.R
#
# tests discretize.jointly function
# Created by: Jiandong Wang, Sajal Kumar and Dr. Mingzhou (Joe) Song
# Date Created: 9th March, 2020

library(testthat)
library(FunChisq)
library(GridOnClusters)
library(cluster)
library(dqrng)

context("Testing GridOnClusters")

test_that("Testing discretize.jointly", {

  # test 1
  # y = f(x)
  # z = f(x)
  # k = constant

  dqset.seed(123)
  x = dqrnorm(100, mean=5, sd=1)
  y = sin(x)
  z = cos(x)
  data = cbind(x, y, z)
  discr = discretize.jointly(data, k=3)

  # test marginal levels
  expect_equivalent(length(unique(discr$D[,1])), 3)
  expect_equivalent(length(unique(discr$D[,2])), 2)
  expect_equivalent(length(unique(discr$D[,3])), 2)

  # test marginal distribution
  expect_equivalent(table(discr$D[,1]), as.table(c(39, 43, 18)))
  expect_equivalent(table(discr$D[,2]), as.table(c(79, 21)))
  expect_equivalent(table(discr$D[,3]), as.table(c(39, 61)))

  # test 2d joint distributions
  dim12 = table(discr$D[,1], discr$D[,2])
  expect_equivalent(dim12, as.table(matrix(c(36, 3,
                                             43, 0,
                                             0, 18),
                                           nrow=3, ncol=2, byrow = T)))

  dim13 = table(discr$D[,1], discr$D[,3])
  expect_equivalent(dim13, as.table(matrix(c(39, 0,
                                             0, 43,
                                             0, 18),
                                           nrow=3, ncol=2, byrow = T)))

  dim23 = table(discr$D[,2], discr$D[,3])
  expect_equivalent(dim23, as.table(matrix(c(36, 43,
                                             3, 18),
                                           nrow=2, ncol=2, byrow = T)))

  # test ARI score
  expect_equivalent(round(discr$csimilarity, digits = 3), 1)

  # test 2
  # y = f(x)
  # z = f(x)
  # k = variable (determined by silhouette)
  dqset.seed(321)

  x = dqrnorm(n = 100, mean=10, sd=2)
  y = log(x)
  z = tan(x)
  data = cbind(x, y, z)
  discr = discretize.jointly(data, k=c(3:10))

  # test marginal levels
  expect_equivalent(length(unique(discr$D[,1])), 7)
  expect_equivalent(length(unique(discr$D[,2])), 7)
  expect_equivalent(length(unique(discr$D[,3])), 4)

  # test marginal distribution
  expect_equivalent(table(discr$D[,1]), as.table(c(17, 6, 11, 34, 6, 4, 22)))
  expect_equivalent(table(discr$D[,2]), as.table(c(17, 6, 11, 34, 6, 4, 22)))
  expect_equivalent(table(discr$D[,3]), as.table(c(2, 3, 82, 13)))

  # test 2d joint distributions
  dim12 = table(discr$D[,1], discr$D[,2])
  expect_equivalent(dim12, as.table(matrix(c(17, 0, 0, 0, 0, 0, 0,
                                             0, 6, 0, 0, 0, 0, 0,
                                             0, 0, 11, 0, 0, 0, 0,
                                             0, 0, 0, 34, 0, 0, 0,
                                             0, 0, 0, 0, 6, 0, 0,
                                             0, 0, 0, 0, 0, 4, 0,
                                             0, 0, 0, 0, 0, 0, 22),
                                           nrow=7, ncol=7, byrow = T)))
  dim13 = table(discr$D[,1], discr$D[,3])
  expect_equivalent(dim13, as.table(matrix(c(0, 0, 16, 1,
                                             1, 0, 0, 5,
                                             0, 0, 11, 0,
                                             0, 0, 34, 0,
                                             0, 0, 0, 6,
                                             1, 3, 0, 0,
                                             0, 0, 21, 1),
                                           nrow=7, ncol=4, byrow = T)))
  dim23 = table(discr$D[,2], discr$D[,3])
  expect_equivalent(dim23, as.table(matrix(c(0, 0, 16, 1,
                                             1, 0, 0, 5,
                                             0, 0, 11, 0,
                                             0, 0, 34, 0,
                                             0, 0, 0, 6,
                                             1, 3, 0, 0,
                                             0, 0, 21, 1),
                                           nrow=7, ncol=4, byrow = T)))

  # test ARI score
  expect_equivalent(round(discr$csimilarity, digits = 3), 0.921)


  # test 3
  # y != f(x)
  # z = f(x, y)
  # k = variable (determined by silhouette)
  dqset.seed(1234)

  x = dqrexp(n=50, rate = 0.6)
  y = dqrnorm(50, mean=2, sd=0.5)
  z = sin(x) + cos(y)
  data = cbind(x, y, z)
  discr = discretize.jointly(data, k=c(3:10))

  # test marginal levels
  expect_equivalent(length(unique(discr$D[,1])), 3)
  expect_equivalent(length(unique(discr$D[,2])), 2)
  expect_equivalent(length(unique(discr$D[,3])), 2)

  # test marginal distribution
  expect_equivalent(table(discr$D[,1]), as.table(c(34, 12, 4)))
  expect_equivalent(table(discr$D[,2]), as.table(c(21, 29)))
  expect_equivalent(table(discr$D[,3]), as.table(c(9, 41)))

  # test 2d joint distributions
  dim12 = table(discr$D[,1], discr$D[,2])
  expect_equivalent(dim12, as.table(matrix(c(13, 21,
                                             5, 7,
                                             3, 1),
                                           nrow=3, ncol=2, byrow = T)))

  dim13 = table(discr$D[,1], discr$D[,3])
  expect_equivalent(dim13, as.table(matrix(c(3, 31,
                                             3, 9,
                                             3, 1),
                                           nrow=3, ncol=2, byrow = T)))

  dim23 = table(discr$D[,2], discr$D[,3])
  expect_equivalent(dim23, as.table(matrix(c(3, 18,
                                             6, 23),
                                           nrow=2, ncol=2, byrow = T)))

  # test ARI score
  expect_equivalent(round(discr$csimilarity, digits = 3), 0.992)

  # test 4
  # y = f(x)
  # z = f(x)
  # k = fixed
  # using an alternate clustering strategy
  dqset.seed(2468)

  x = dqrnorm(n = 1000, mean = 10, sd = 2)
  y = sin(x)
  z = cos(y)
  data = cbind(x, y, z)
  # use PAM to cluster
  alt.cluster = pam(x = data, k = 5, diss = FALSE, metric = "euclidean", cluster.only = TRUE)

  discr = discretize.jointly(data = data, cluster_label = alt.cluster)

  # test marginal levels
  expect_equivalent(length(unique(discr$D[,1])), 5)
  expect_equivalent(length(unique(discr$D[,2])), 3)
  expect_equivalent(length(unique(discr$D[,3])), 2)

  # test marginal distribution
  expect_equivalent(table(discr$D[,1]), as.table(c(82, 197, 252, 293, 176)))
  expect_equivalent(table(discr$D[,2]), as.table(c(321, 454, 225)))
  expect_equivalent(table(discr$D[,3]), as.table(c(472, 528)))

  # test 2d joint distributions
  dim12 = table(discr$D[,1], discr$D[,2])
  expect_equivalent(dim12, as.table(matrix(c(18, 59, 5,
                                             0, 14, 183,
                                             0, 252, 0,
                                             292, 1, 0,
                                             11, 128, 37),
                                           nrow=5, ncol=3, byrow = T)))

  dim13 = table(discr$D[,1], discr$D[,3])
  expect_equivalent(dim13, as.table(matrix(c(14, 68,
                                             165, 32,
                                             0, 252,
                                             265, 28,
                                             28, 148),
                                           nrow=5, ncol=2, byrow = T)))

  dim23 = table(discr$D[,2], discr$D[,3])
  expect_equivalent(dim23, as.table(matrix(c(279, 42,
                                             0, 454,
                                             193, 32),
                                           nrow=3, ncol=2, byrow = T)))

  # test ARI score
  expect_equivalent(round(discr$csimilarity, digits = 3), 0.915)

})
