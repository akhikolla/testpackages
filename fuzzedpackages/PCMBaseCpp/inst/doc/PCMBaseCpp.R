## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ------------------------------------------------------------------------
library(PCMBase)
library(PCMBaseCpp)

system.time(llR <- PCMLik(
  X = PCMBaseTestObjects$traits.ab.123, 
  tree = PCMBaseTestObjects$tree.ab,
  model = PCMBaseTestObjects$model_MixedGaussian_ab))

system.time(llCpp <- PCMLik(
  X = PCMBaseTestObjects$traits.ab.123, 
  tree = PCMBaseTestObjects$tree.ab,
  model = PCMBaseTestObjects$model_MixedGaussian_ab, 
  metaI = PCMInfoCpp))

print(llR)
print(llCpp)

## ------------------------------------------------------------------------
logLikFunR <- PCMCreateLikelihood(
  X = PCMBaseTestObjects$traits.ab.123, 
  tree = PCMBaseTestObjects$tree.ab,
  model = PCMBaseTestObjects$model_MixedGaussian_ab)

logLikFunCpp <- PCMCreateLikelihood(
  X = PCMBaseTestObjects$traits.ab.123, 
  tree = PCMBaseTestObjects$tree.ab,
  model = PCMBaseTestObjects$model_MixedGaussian_ab, 
  metaI = PCMInfoCpp)

metaICpp <- PCMInfoCpp(
  X = PCMBaseTestObjects$traits.ab.123, 
  tree = PCMBaseTestObjects$tree.ab,
  model = PCMBaseTestObjects$model_MixedGaussian_ab)

logLikFunCpp2 <- PCMCreateLikelihood(
  X = PCMBaseTestObjects$traits.ab.123, 
  tree = PCMBaseTestObjects$tree.ab,
  model = PCMBaseTestObjects$model_MixedGaussian_ab, 
  metaI = metaICpp)

set.seed(1, kind = "Mersenne-Twister", normal.kind = "Inversion")
randParam <- PCMParamRandomVecParams(PCMBaseTestObjects$model_MixedGaussian_ab)

system.time(llR <- logLikFunR(randParam))

system.time(llCpp <- logLikFunCpp(randParam))

system.time(llCpp2 <- logLikFunCpp2(randParam))

print(llR)
print(llCpp)
print(llCpp2)

## ------------------------------------------------------------------------
metaIR <- PCMInfo(
  X = PCMBaseTestObjects$traits.ab.123, 
  tree = PCMBaseTestObjects$tree.ab,
  model = PCMBaseTestObjects$model_MixedGaussian_ab)

metaICpp <- PCMInfoCpp(
  X = PCMBaseTestObjects$traits.ab.123, 
  tree = PCMBaseTestObjects$tree.ab,
  model = PCMBaseTestObjects$model_MixedGaussian_ab)

system.time(llR <- PCMLik(
  X = PCMBaseTestObjects$traits.ab.123, 
  tree = PCMBaseTestObjects$tree.ab,
  model = PCMBaseTestObjects$model_MixedGaussian_ab, 
  metaI = metaIR))

system.time(llCpp <- PCMLik(
  X = PCMBaseTestObjects$traits.ab.123, 
  tree = PCMBaseTestObjects$tree.ab,
  model = PCMBaseTestObjects$model_MixedGaussian_ab, 
  metaI = metaICpp))

print(llR)
print(llCpp)

## ------------------------------------------------------------------------
logLikFunR <- PCMCreateLikelihood(
  X = PCMBaseTestObjects$traits.ab.123, 
  tree = PCMBaseTestObjects$tree.ab,
  model = PCMBaseTestObjects$model_MixedGaussian_ab, 
  metaI = metaIR)

logLikFunCpp <- PCMCreateLikelihood(
  X = PCMBaseTestObjects$traits.ab.123, 
  tree = PCMBaseTestObjects$tree.ab,
  model = PCMBaseTestObjects$model_MixedGaussian_ab, 
  metaI = metaICpp)

system.time(llR <- PCMLik(
  X = PCMBaseTestObjects$traits.ab.123, 
  tree = PCMBaseTestObjects$tree.ab,
  model = PCMBaseTestObjects$model_MixedGaussian_ab, 
  metaI = metaIR))

system.time(llCpp <- PCMLik(
  X = PCMBaseTestObjects$traits.ab.123, 
  tree = PCMBaseTestObjects$tree.ab,
  model = PCMBaseTestObjects$model_MixedGaussian_ab, 
  metaI = metaICpp))

print(llR)
print(llCpp)

