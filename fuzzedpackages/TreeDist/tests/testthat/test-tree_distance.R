context('tree_distance.R')

# Labels in different order to confound as.Splits
treeSym8 <- ape::read.tree(text='((e, (f, (g, h))), (((a, b), c), d));')
treeBal8 <- ape::read.tree(text='(((e, f), (g, h)), ((a, b), (c, d)));')
treeOpp8 <- ape::read.tree(text='(((a, f), (c, h)), ((g, b), (e, d)));')
treesSBO8 <- structure(list(treeSym8, treeBal8, treeOpp8), 
                            class = 'multiPhylo')
treesSSBB8 <- structure(list(treeSym8, treeSym8, treeBal8, treeBal8), 
                            class = 'multiPhylo')

treeCat8 <- ape::read.tree(text='((((h, g), f), e), (d, (c, (b, a))));')
treeTac8 <- ape::read.tree(text='((((e, c), g), a), (h, (b, (d, f))));')
treeStar8 <- ape::read.tree(text='(e, c, g, h, b, a, d, f);')

treeAb.Cdefgh <- ape::read.tree(text='((a, b), (c, d, e, f, g, h));')
treeAbc.Defgh <- ape::read.tree(text='((a, b, c), (d, e, f, g, h));')
treeAcd.Befgh <- ape::read.tree(text='((a, c, d), (b, e, f, g, h));')
treeAbcd.Efgh <- ape::read.tree(text='((a, b, c, d), (e, f, g, h));')
treeTwoSplits <- ape::read.tree(text="(((a, b), c, d), (e, f, g, h));")

test_that("Split compatibility is correctly established", {
  expect_true(SplitsCompatible(as.logical(c(0,0,1,1,0)), 
                               as.logical(c(0,0,1,1,0))))
  expect_true(SplitsCompatible( as.logical(c(0,0,1,1,0)), 
                               !as.logical(c(0,0,1,1,0))))
  expect_true(SplitsCompatible(as.logical(c(0,0,1,1,0)), 
                               as.logical(c(1,0,1,1,0))))
  expect_true(SplitsCompatible(!as.logical(c(0,0,1,1,0)), 
                                as.logical(c(1,0,1,1,0))))
  expect_false(SplitsCompatible(as.logical(c(0,0,1,1,0)), 
                                as.logical(c(1,1,0,1,0))))
})

methodsToTest <- list(
  SharedPhylogeneticInfo,
  DifferentPhylogeneticInfo,
  MatchingSplitInfo,
  MatchingSplitInfoDistance,
  MutualClusteringInfo,
  ClusteringInfoDistance,
  NyeSimilarity,
  JaccardRobinsonFoulds,
  MatchingSplitDistance,
  RobinsonFoulds,
  InfoRobinsonFoulds,
  KendallColijn # List last: requires rooted trees.
)

NormalizationTest <- function (FUNC, ...) {
  expect_equal(c(1L, 1L), 
               FUNC(treesSSBB8, normalize = TRUE, ...)[c(1, 6)],
               tolerance = 1e-7)
}

test_that('Bad labels cause error', {
  treeBadLabel8 <- ape::read.tree(text='((a, b, c, D), (e, f, g, h));')
  lapply(methodsToTest, function(Func) 
    expect_error(Func(treeSym8, treeBadLabel8)))
})

test_that('Size mismatch causes error', {
  treeSym7 <- ape::read.tree(text='((e, (f, g)), (((a, b), c), d));')
  splits7 <- as.Splits(treeSym7)
  splits8 <- as.Splits(treeSym8)
  
  lapply(methodsToTest, function(Func) 
    expect_error(Func(treeSym8, treeSym7)))
  
  lapply(methodsToTest, function(Func) 
    expect_error(Func(treeSym7, treeSym8)))
  
  expect_error(MeilaVariationOfInformation(splits7, splits8))
  
  Test <- function (Func) {
    expect_error(Func(splits8, as.Splits(BalancedTree(9)), 8))
  }
  Test(cpp_robinson_foulds_distance)
  Test(cpp_robinson_foulds_info)
  Test(cpp_matching_split_distance)
  Test(cpp_jaccard_similarity)
  Test(cpp_mmsi_distance)
  Test(cpp_mutual_clustering)
  Test(cpp_shared_phylo)
})

test_that('Metrics handle polytomies', {
  polytomy8 <- ape::read.tree(text='(a, b, c, d, e, f, g, h);')
  lapply(list(SharedPhylogeneticInfo, MutualClusteringInfo,
              MatchingSplitDistance, NyeSimilarity),
         function (Func) expect_equal(0, Func(treeSym8, polytomy8)))
})

#Func <- ClusteringInfoDistance # FUNC =
test_that('Output dimensions are correct', {
  list1 <- list(sym = treeSym8, bal = treeBal8)
  list2 <- list(sym = treeSym8, abc = treeAbc.Defgh, abcd = treeAbcd.Efgh)
  dimNames <- list(c('sym', 'bal'), c('sym', 'abc', 'abcd'))
  
  Test <- function (Func) {
    allPhylo <- matrix(
      c(Func(treeSym8, treeSym8),      Func(treeBal8, treeSym8),
        Func(treeSym8, treeAbc.Defgh), Func(treeBal8, treeAbc.Defgh),
        Func(treeSym8, treeAbcd.Efgh), Func(treeBal8, treeAbcd.Efgh)),
      2L, 3L, dimnames = dimNames)
    phylo1 <- matrix(c(Func(treeSym8, list2), Func(treeBal8, list2)),
                     byrow = TRUE, 2L, 3L, dimnames = dimNames)
    phylo2 <- matrix(c(Func(list1, treeSym8), Func(list1, treeAbc.Defgh),
                       Func(list1, treeAbcd.Efgh)), 2L, 3L, dimnames = dimNames)
    noPhylo <- Func(list1, list2)
    expect_equal(allPhylo, phylo1)
    expect_equal(allPhylo, phylo2)
    expect_equal(allPhylo, noPhylo)
  }
  
  lapply(methodsToTest, Test)
})

test_that('Robinson Foulds Distance is correctly calculated', {
  RFTest <- function (t1, t2) {
    expect_equal(suppressMessages(phangorn::RF.dist(t1, t2)),
                 RobinsonFoulds(t1, t2))
    
    expected <- RobinsonFoulds(t1, t2, reportMatching = TRUE, similarity = TRUE)
    attr(expected, 'pairScores') <- attr(expected, 'pairScores') == 0L
    expect_equal(expected, RobinsonFouldsMatching(t1, t2))
  }
  RFTest(treeSym8, treeSym8)
  RFTest(treeSym8, treeStar8)
  RFTest(treeStar8, treeStar8)
  RFTest(treeAb.Cdefgh, treeAbc.Defgh)
  RFTest(treeAb.Cdefgh, treeAbcd.Efgh)
  
  # Invariant to tree description order
  sq_pectinate <- ape::read.tree(text='((((((1, 2), 3), 4), 5), 6), (7, (8, (9, (10, 11)))));')
  shuffle1 <- ape::read.tree(text='(((((1, 5), 2), 6), (3, 4)), ((8, (7, 9)), (10, 11)));')
  shuffle2 <- ape::read.tree(text='(((8, (7, 9)), (10, 11)), ((((1, 5), 2), 6), (3, 4)));')
  RFTest(shuffle1, sq_pectinate)
  RFTest(sq_pectinate, shuffle1)
  RFTest(shuffle1, shuffle2)
  RFTest(shuffle1, sq_pectinate)
  RFTest(shuffle2, sq_pectinate)
})

test_that('Shared Phylogenetic Info is correctly calculated', {
  
  expect_equal(5.529821, tolerance = 1e-7,
               cpp_shared_phylo(
                 as.Splits(as.logical(c(1, 1, 1, 1, 0, 0, 0, 0))),
                 as.Splits(as.logical(c(1, 1, 1, 1, 0, 0, 0, 0))),
                 8L)$score)
  
  expect_equal(0.2895066, tolerance = 1e-7,
               cpp_shared_phylo(
                 as.Splits(as.logical(c(1, 1, 0, 0, 0, 0, 0, 0))),
                 as.Splits(as.logical(c(0, 0, 1, 1, 0, 0, 0, 0))),
                 8L)$score)
  
  expect_equal(1.137504, tolerance = 1e-6,
               cpp_shared_phylo(
                 as.Splits(as.logical(c(1, 1, 0, 0, 0, 0, 0, 0))),
                 as.Splits(as.logical(c(1, 1, 1, 1, 0, 0, 0, 0))),
                 8L)$score)
  
  expect_equal(3.45943, tolerance = 1e-6,
               cpp_shared_phylo(
                 as.Splits(as.logical(c(1, 1, 0, 0, 0, 0, 0, 0))),
                 as.Splits(as.logical(c(1, 1, 0, 0, 0, 0, 0, 0))),
                 8L)$score)
  
  expect_equal(22.53747, tolerance = 1e-05,
               SharedPhylogeneticInfo(treeSym8, treeSym8, normalize = FALSE))
  
  expect_equal(1, tolerance = 1e-05,
               SharedPhylogeneticInfo(treeSym8, treeSym8, normalize = TRUE))
  
  expect_equal(0,
               SharedPhylogeneticInfo(treeSym8, treeStar8, normalize = TRUE))

  expect_equal(0,
               SharedPhylogeneticInfo(treeStar8, treeStar8, normalize = FALSE))

  expect_equal(NaN, # Division by zero
               SharedPhylogeneticInfo(treeStar8, treeStar8, normalize = TRUE))
  
  expect_equal(13.75284, SharedPhylogeneticInfo(treeSym8, treeBal8), tolerance=1e-05)
  
  expect_equal(DifferentPhylogeneticInfo(treeSym8, treeAcd.Befgh),
               DifferentPhylogeneticInfo(treeAcd.Befgh, treeSym8), tolerance=1e-05)
  
  expect_equal(0, DifferentPhylogeneticInfo(treeSym8, treeSym8, normalize = TRUE))
  
  infoSymBal <- SplitwiseInfo(treeSym8) + SplitwiseInfo(treeBal8)
  expect_equal(infoSymBal - 13.75284 - 13.75284, tolerance = 1e-05,
    DifferentPhylogeneticInfo(treeSym8, treeBal8, normalize = TRUE) * infoSymBal)
  
  expect_equal(22.53747 + SharedPhylogeneticInfo(treeAcd.Befgh, treeAcd.Befgh) - 
                 (2 * SharedPhylogeneticInfo(treeSym8, treeAcd.Befgh)), 
               DifferentPhylogeneticInfo(treeSym8, treeAcd.Befgh), 
               tolerance=1e-06)
  
