library(testthat)
context("PCMLik, White")

library(PCMBase)
library(PCMBaseCpp)

if(PCMBaseCppIsADevRelease()) {
  
  load("testobjects.RData")
  
  set.seed(1)
  
  test_that("Equal R and Cpp likelihood on a random model, single regime (a)", {
    expect_silent(model.a.123.White <- PCM("White", k = 3, regimes = "a"))
    expect_silent(PCMParamLoadOrStore(model.a.123.White, 
                                      PCMParamRandomVecParams(model.a.123.White), 
                                      offset = 0, k = 3, load = TRUE))
    expect_silent(metaICpp <- PCMInfoCpp(X = traits.a.123[, 1:length(tree.a$tip.label)],
                                         tree = tree.a, model.a.123.White))
    expect_equal(PCMLik(traits.a.123, tree.a, model.a.123.White),
                 PCMLik(traits.a.123, tree.a, model.a.123.White, metaI = metaICpp))
    
  })
  
  test_that("Equal R and Cpp likelihood on a random model, multiple regimes (ab)", {
    expect_silent(model.ab.123.White <- PCM("White", k = 3, regimes = c("a", "b")))
    expect_silent(PCMParamLoadOrStore(model.ab.123.White, 
                                      PCMParamRandomVecParams(model.ab.123.White), 
                                      offset = 0, k = 3, load = TRUE))
    expect_silent(metaICpp <- PCMInfoCpp(X = traits.ab.123[, 1:length(tree.ab$tip.label)],
                                         tree = tree.ab, model.ab.123.White))
    expect_equal(PCMLik(traits.ab.123, tree.ab, model.ab.123.White),
                 PCMLik(traits.ab.123, tree.ab, model.ab.123.White, metaI = metaICpp))
    
  })
}
