library(testthat)
context("PCMLik, OU")

library(PCMBase)
library(PCMBaseCpp)

if(PCMBaseCppIsADevRelease()) {
  
  load("testobjects.RData")
  
  
  test_that("Equal R and Cpp likelihood on a fixed model, single regime (a)", {
    expect_silent(metaICpp <- PCMInfoCpp(X = traits.a.123[, 1:length(tree.a$tip.label)],
                                         tree = tree.a, model.a.123))
    expect_equal(PCMLik(traits.a.123, tree.a, model.a.123),
                 PCMLik(traits.a.123, tree.a, model.a.123, metaI = metaICpp))
    expect_equal(PCMLik(traits.ab.123, tree.ab, model.ab.123),
                 PCMLik(traits.ab.123, tree.ab, model.ab.123,
                        metaI = PCMInfoCpp(X = traits.ab.123[, 1:length(tree.ab$tip.label)],
                                           tree = tree.ab,
                                           model.ab.123)))
  })
  
  set.seed(1)
  
  test_that("Equal R and Cpp likelihood on a random model, single regime (a)", {
    expect_silent(model.a.123.OU <- PCM("OU", k = 3, regimes = "a"))
    expect_silent(PCMParamLoadOrStore(model.a.123.OU, 
                                      PCMParamRandomVecParams(model.a.123.OU), 
                                      offset = 0, k = 3, load = TRUE))
    expect_silent(metaICpp <- PCMInfoCpp(X = traits.a.123[, 1:length(tree.a$tip.label)],
                                         tree = tree.a, model.a.123.OU))
    expect_equal(PCMLik(traits.a.123, tree.a, model.a.123.OU),
                 PCMLik(traits.a.123, tree.a, model.a.123.OU, metaI = metaICpp))
    
  })
  
  test_that("Equal R and Cpp likelihood on a random model, multiple regimes (ab)", {
    expect_silent(model.ab.123.OU <- PCM("OU", k = 3, regimes = c("a", "b")))
    expect_silent(PCMParamLoadOrStore(model.ab.123.OU, 
                                      PCMParamRandomVecParams(model.ab.123.OU), 
                                      offset = 0, k = 3, load = TRUE))
    expect_silent(metaICpp <- PCMInfoCpp(X = traits.ab.123[, 1:length(tree.ab$tip.label)],
                                         tree = tree.ab, model.ab.123.OU))
    expect_equal(PCMLik(traits.ab.123, tree.ab, model.ab.123.OU),
                 PCMLik(traits.ab.123, tree.ab, model.ab.123.OU, metaI = metaICpp))
    
    expect_silent(metaICpp <- PCMInfoCpp(
      X = traits.ab.123[, 1:length(tree.ab$tip.label)],
      tree = tree.ab, model.ab.123.OU,
      SE = 0.05*abs(traits.ab.123[, 1:length(tree.ab$tip.label)])))
    expect_equal(PCMLik(traits.ab.123, tree.ab, model.ab.123.OU, SE = 0.05*abs(traits.ab.123[, 1:length(tree.ab$tip.label)])),
                 PCMLik(traits.ab.123, tree.ab, model.ab.123.OU, metaI = metaICpp))
  })
  
}

