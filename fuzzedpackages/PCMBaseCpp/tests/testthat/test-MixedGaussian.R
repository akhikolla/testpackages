library(testthat)
context("PCMLik, MixedGaussian")

library(PCMBase)
library(PCMBaseCpp)
library(data.table)

if(PCMBaseCppIsADevRelease()) {
  
  load("testobjects.RData")
  
  test_that("Calling PCMGenerateParameterizations()", {
    expect_silent(tableParametrizationsOU <- PCMTableParameterizations(structure(0.0, class="OU")))
    expect_true(is.data.table(tableParametrizationsOU))
    expect_silent(
      PCMGenerateParameterizations(
        model = structure(0.0, class="OU"),
        
        # note that I am not using data.table but data.frame syntax for subsetting
        # tableParameterizationsOU. This to avoid a problem with devtools::test
        # see https://github.com/r-lib/devtools/issues/192
        # Another work-around would be to add data.table to Depends:, but I don't 
        # want this now.
        tableParameterizations = tableParametrizationsOU[
          sapply(tableParametrizationsOU$X0, function(type)
            identical(type, c("VectorParameter", "_Global")) ||
              identical(type, c("VectorParameter", "_Omitted"))
          ) &
            sapply(tableParametrizationsOU$H, function(type)
              identical(type, c("MatrixParameter"))) &
            sapply(tableParametrizationsOU$Theta, function(type)
              identical(type, "VectorParameter") ), ])
    )
    expect_silent(tableParametrizationsBM <- PCMTableParameterizations(structure(0.0, class="BM")))
    expect_true(is.data.table(tableParametrizationsBM))
    expect_silent(
      PCMGenerateParameterizations(
        model = structure(0.0, class="BM"),
        tableParameterizations = tableParametrizationsBM[
          sapply(tableParametrizationsBM$X0, function(type)
            identical(type, c("VectorParameter", "_Global")) ||
              identical(type, c("VectorParameter", "_Omitted")) ), ])
    )
  })
  
  
  test_that("Equal OU and Cpp likelihoods on a fixed MixedGaussian model", {
    expect_equal(
      PCMLik(traits.ab.123, tree.ab, model_MixedGaussian_ab_globalSigmae_x),
      PCMLik(traits.ab.123, tree.ab, model_MixedGaussian_ab_globalSigmae_x, 
             metaI = PCMInfoCpp(traits.ab.123, tree.ab, model_MixedGaussian_ab_globalSigmae_x)
             )
    )
    
  })
  
  
  
  set.seed(1)
  
  test_that("Equal R and Cpp likelihood on a random MixedGaussian model", {
    expect_silent(model.ab.123.MG <- MixedGaussian(
      k = 3, 
      modelTypes = c("BM__Omitted_X0__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x", 
                     "OU__Omitted_X0__H__Theta__UpperTriangularWithDiagonal_WithNonNegativeDiagonal_Sigma_x__Omitted_Sigmae_x"), 
      mapping = c(a=2L, b=1L), className = "MG"))
    expect_silent(PCMParamLoadOrStore(model.ab.123.MG, 
                                      PCMParamRandomVecParams(model.ab.123.MG), 
                                      offset = 0, k = 3, load = TRUE))
    expect_silent(metaICpp <- PCMInfoCpp(X = traits.ab.123[, 1:length(tree.ab$tip.label)],
                                         tree = tree.ab, model.ab.123.MG))
    expect_equal(PCMLik(traits.ab.123, tree.ab, model.ab.123.MG),
                 PCMLik(traits.ab.123, tree.ab, model.ab.123.MG, metaI = metaICpp))
    
  })
}

