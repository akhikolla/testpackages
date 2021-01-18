context("uspr.R")
library('TreeTools')

test_that("Bad tree input is handled correctly", {
  expect_error(USPRDist(PectinateTree(1:8), PectinateTree(2:9)))
  expect_error(USPRDist(PectinateTree(1:8), PectinateTree(1:9)))
  expect_error(USPRDist(PectinateTree(1:9), PectinateTree(1:8)))
  list2 <- list(PectinateTree(1:8), BalancedTree(1:8))
  list3 <- list(PectinateTree(1:8), BalancedTree(1:8), BalancedTree(1:8))
  expect_error(USPRDist(list2, list3))
  expect_error(USPRDist(list2, list3, checks = FALSE))

  nwk2 <- as.Newick(list2)
  nwk3 <- as.Newick(list3)
  expect_error(replug_dist(nwk2, nwk3))
  expect_error(tbr_dist(nwk2, nwk3, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE))
  expect_error(uspr_dist(nwk2, nwk3, FALSE, FALSE, FALSE))
})

test_that("TBR options", {
  bal8 <- BalancedTree(8)
  pec8 <- PectinateTree(8)
  expect_warning(TBRDist(bal8, pec8, maf = TRUE, exact = FALSE))
  expect_warning(
    expect_equal(list(tbr_exact = 0L, n_maf = 1L,
                      maf_1 = "(((t1,t2),(t3,t4)),(t7,t8),(t5,t6));",
                      maf_2 = "(((t1,t2),(t3,t4)),(t7,t8),(t5,t6));"),
                 TBRDist(bal8, maf = TRUE, countMafs = TRUE)))
  expect_warning(expect_equal(0L, TBRDist(bal8, exact = TRUE)))
  expect_equal(list(tbr_min = 1, tbr_max = 3, n_maf = 13),
               TBRDist(bal8, pec8, optimize = FALSE, countMafs = TRUE))
})

test_that("SPR distances are calculated correctly", {
  tree1 <- BalancedTree(10)
  tree2 <- PectinateTree(10)
  treeR <- ape::read.tree(text="(t1, (((t5, t7), (t9, (t3, t2))), (t4, ((t6, t8), t10))));")
  list1 <- list(one = tree1, oneAgain = tree1, two = tree2, three = treeR)
  list2 <- list(tree1, tree2)
  expect_equivalent(2L, USPRDist(tree1, tree2))
  expect_equivalent(c(0, 2L), USPRDist(list(tree1, tree2), tree1))
  expect_equivalent(c(0, 2L), USPRDist(tree1, list(tree1, tree2)))

  goodRet <- structure(c(0, 2, 4, 2, 4, 4),
                       Size = 4L,
                       Labels = names(list1),
                       Diag = FALSE,
                       Upper = FALSE,
                       class = 'dist')
  expect_equal(goodRet, USPRDist(list1))

  expect_equal(goodRet, ReplugDist(list1))
  first <- ReplugDist(list1, list1[[1]], maf = TRUE)
  each <- ReplugDist(list1, maf = TRUE)
  expect_equivalent(first[[1]], as.matrix(each[[1]])[, 1])
  expect_equivalent(first[[2]][-1], each[[2]][-1, 1])
  expect_equivalent(first[[3]][-1], each[[3]][-1, 1])

  expect_equivalent(c(0, 1, 4, 1, 4, 4),
                    as.integer(TBRDist(list1, exact = TRUE)))
  expect_equivalent(c(0, 1, 3, 1, 3, 3),
                    as.integer(TBRDist(list1, exact = FALSE)$tbr_min))
  first <- TBRDist(list1, list1[[1]], maf = TRUE)
  each <- TBRDist(list1, maf = TRUE)
  expect_equivalent(first[[1]], as.matrix(each[[1]])[1, ])
  # Non-numerics will be replaced with NAs, so [-1]:
  expect_equivalent(first[[2]][-1], as.matrix(each[[2]])[-1, 1])

  expect_equal(matrix(c(0, 2, 0, 2, 2, 0, 4, 4), 4, 2, byrow = TRUE,
                      dimnames = list(names(list1), names(list2))),
               ReplugDist(list1, list2, allPairs = TRUE))
  expect_equivalent(ReplugDist(list1, tree2, maf = TRUE)$maf_2,
                    ReplugDist(list1, list2, allPairs = TRUE, maf = TRUE)$maf_2[, 2])


  Test <- function (tree1, tree2) {
    td <- TBRDist(tree1, tree2, exact = TRUE, approximate = TRUE)
    expect_true(USPRDist(tree1, tree2) >= TBRDist(tree1, tree2, exact = TRUE))
    expect_true(td$tbr_exact >= td$tbr_min)
    expect_true(td$tbr_exact <= td$tbr_max)

    td4 <- TBRDist(list(tree1, tree1, tree2, tree2),
                   list(tree1, tree2, tree1, tree2), exact = TRUE)
    expect_equal(0L, td4[[1]])
    expect_equal(td4[[2]], td4[[3]])
    expect_equal(0L, td4[[4]])

    sd4 <- USPRDist(list(tree1, tree1, tree2, tree2),
                    list(tree1, tree2, tree1, tree2))
    expect_equal(0L, sd4[[1]])
    expect_equal(sd4[[2]], sd4[[3]])
    expect_equal(0L, sd4[[4]])

    rd4 <- ReplugDist(list(tree1, tree1, tree2, tree2),
                      list(tree1, tree2, tree1, tree2))
    expect_equal(0L, rd4[[1]])
    expect_equal(rd4[[2]], rd4[[3]])
    expect_equal(0L, rd4[[4]])

  }

  Test(tree1, tree2)
  Test(PectinateTree(13), BalancedTree(13))

  expect_warning(
    expect_equal(invisible(),
               TBRDist(tree1, tree2, exact = FALSE, approximate = FALSE))
  )
})

test_that("MAF info is calculated", {
  tree1 <- BalancedTree(8)
  tree2 <- PectinateTree(8)
  expect_equal(LnUnrooted(8) / log(2), MAFInfo(tree1, tree1))
  expect_lt(MAFInfo(tree1, tree2), LnUnrooted(8) / log(2))
})
