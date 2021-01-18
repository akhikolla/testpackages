context("attachEpiNet")
library(epinetr)

test_that("Sole epistatic effects give correct mean and variance", {
  pop <- Population(
    popSize = 10, map = map100snp, QTL = 20, alleleFrequencies = runif(100), broadH2 = 1,
    narrowh2 = 0, traitVar = 40
  )
  pop <- attachEpiNet(pop)
  components <- getComponents(pop)
  expect_equal(mean(components$Epistatic), 0)
  expect_equal(var(components$Epistatic), 40)
  expect_equal(mean(components$Phenotype), 0)
  expect_equal(var(components$Phenotype), 40)
})

test_that("Epistatic and environmental effects give correct means and variances", {
  pop <- Population(
    popSize = 10, map = map100snp, QTL = 20, alleleFrequencies = runif(100), broadH2 = 0.8,
    narrowh2 = 0, traitVar = 40
  )
  pop <- attachEpiNet(pop)
  components <- getComponents(pop)
  expect_equal(mean(components$Epistatic), 0)
  expect_equal(var(components$Epistatic), 32)
  expect_equal(mean(components$Environmental), 0)
  expect_equal(mean(components$Phenotype), 0)
  # expect_equal(var(components$Phenotype), 40)
})

test_that("Epistatic, environmental and additive effects give correct means and variances", {
  pop <- Population(
    popSize = 10, map = map100snp, QTL = 20, alleleFrequencies = runif(100), broadH2 = 0.8,
    narrowh2 = 0.5, traitVar = 40
  )
  pop <- attachEpiNet(pop)
  pop <- addEffects(pop)
  components <- getComponents(pop)
  expect_equal(mean(components$Additive), 0)
  expect_equal(var(components$Additive), 20)
  expect_equal(mean(components$Epistatic), 0)
  expect_equal(var(components$Epistatic), 12)
  # expect_equal(var(components$Additive + components$Epistatic), 32)
  expect_equal(mean(components$Environmental), 0)
  expect_equal(mean(components$Phenotype), 0)
  # expect_equal(var(components$Phenotype), 40)
})

test_that("Changing broad-sense heritability changes variances accordingly", {
  pop <- Population(
    popSize = 10, map = map100snp, QTL = 20, alleleFrequencies = runif(100), broadH2 = 0.8,
    narrowh2 = 0.5, traitVar = 40
  )
  pop <- attachEpiNet(pop)
  pop <- addEffects(pop)
  pop <- Population(pop, broadH2 = 0.7)
  components <- getComponents(pop)
  expect_equal(mean(components$Additive), 0)
  expect_equal(var(components$Additive), 20)
  expect_equal(mean(components$Epistatic), 0)
  expect_equal(var(components$Epistatic), 8)
  # expect_equal(var(components$Additive + components$Epistatic), 28)
  expect_equal(mean(components$Environmental), 0)
  expect_equal(mean(components$Phenotype), 0)
  # expect_equal(var(components$Phenotype), 40)
})

test_that("Changing narrow-sense heritability changes variances accordingly", {
  pop <- Population(
    popSize = 10, map = map100snp, QTL = 20, alleleFrequencies = runif(100), broadH2 = 0.8,
    narrowh2 = 0.5, traitVar = 40
  )
  pop <- attachEpiNet(pop)
  pop <- addEffects(pop)
  pop <- Population(pop, narrowh2 = 0.3)
  components <- getComponents(pop)
  expect_equal(mean(components$Additive), 0)
  expect_equal(var(components$Additive), 12)
  expect_equal(mean(components$Epistatic), 0)
  expect_equal(var(components$Epistatic), 20)
  # expect_equal(var(components$Additive + components$Epistatic), 32)
  expect_equal(mean(components$Environmental), 0)
  expect_equal(mean(components$Phenotype), 0)
  # expect_equal(var(components$Phenotype), 40)
})

test_that("Changing narrow/broad-sense heritability changes variances accordingly", {
  pop <- Population(
    popSize = 10, map = map100snp, QTL = 20, alleleFrequencies = runif(100), broadH2 = 0.8,
    narrowh2 = 0.5, traitVar = 40
  )
  pop <- attachEpiNet(pop)
  pop <- addEffects(pop)
  pop <- Population(pop, broadH2 = 0.4, narrowh2 = 0.3)
  components <- getComponents(pop)
  expect_equal(mean(components$Additive), 0)
  expect_equal(var(components$Additive), 12)
  expect_equal(mean(components$Epistatic), 0)
  expect_equal(var(components$Epistatic), 4)
  # expect_equal(var(components$Additive + components$Epistatic), 16)
  expect_equal(mean(components$Environmental), 0)
  expect_equal(mean(components$Phenotype), 0)
  # expect_equal(var(components$Phenotype), 40)
})

test_that("Changing broad-sense heritability & trait variance changes variances accordingly", {
  pop <- Population(
    popSize = 10, map = map100snp, QTL = 20, alleleFrequencies = runif(100), broadH2 = 0.8,
    narrowh2 = 0.5, traitVar = 40
  )
  pop <- attachEpiNet(pop)
  pop <- addEffects(pop)
  pop <- Population(pop, broadH2 = 0.7, traitVar = 20)
  components <- getComponents(pop)
  expect_equal(mean(components$Additive), 0)
  expect_equal(var(components$Additive), 10)
  expect_equal(mean(components$Epistatic), 0)
  expect_equal(var(components$Epistatic), 4)
  # expect_equal(var(components$Additive + components$Epistatic), 14)
  expect_equal(mean(components$Environmental), 0)
  expect_equal(mean(components$Phenotype), 0)
  # expect_equal(var(components$Phenotype), 20)
})

test_that("Changing narrow-sense heritability & trait variance changes variances accordingly", {
  pop <- Population(
    popSize = 10, map = map100snp, QTL = 20, alleleFrequencies = runif(100), broadH2 = 0.8,
    narrowh2 = 0.5, traitVar = 40
  )
  pop <- attachEpiNet(pop)
  pop <- addEffects(pop)
  pop <- Population(pop, narrowh2 = 0.3, traitVar = 20)
  components <- getComponents(pop)
  expect_equal(mean(components$Additive), 0)
  expect_equal(var(components$Additive), 6)
  expect_equal(mean(components$Epistatic), 0)
  expect_equal(var(components$Epistatic), 10)
  # expect_equal(var(components$Additive + components$Epistatic), 16)
  expect_equal(mean(components$Environmental), 0)
  expect_equal(mean(components$Phenotype), 0)
  # expect_equal(var(components$Phenotype), 20)
})

test_that("Changing narrow/broad-sense heritability & trait variance changes variances accordingly", {
  pop <- Population(
    popSize = 10, map = map100snp, QTL = 20, alleleFrequencies = runif(100), broadH2 = 0.8,
    narrowh2 = 0.5, traitVar = 40
  )
  pop <- attachEpiNet(pop)
  pop <- addEffects(pop)
  pop <- Population(pop, broadH2 = 0.4, narrowh2 = 0.3, traitVar = 20)
  components <- getComponents(pop)
  expect_equal(mean(components$Additive), 0)
  expect_equal(var(components$Additive), 6)
  expect_equal(mean(components$Epistatic), 0)
  expect_equal(var(components$Epistatic), 2)
  # expect_equal(var(components$Additive + components$Epistatic), 8)
  expect_equal(mean(components$Environmental), 0)
  expect_equal(mean(components$Phenotype), 0)
  # expect_equal(var(components$Phenotype), 20)
})

test_that("Changing narrow or broad-sense heritability so they swap places throws an error", {
  pop <- Population(
    popSize = 10, map = map100snp, QTL = 20, alleleFrequencies = runif(100), broadH2 = 0.8,
    narrowh2 = 0.5, traitVar = 40
  )
  expect_error(Population(pop, broadH2 = 0.4), "Narrow-sense heritability cannot exceed broad-sense heritability")
  expect_error(Population(pop, narrowh2 = 0.9), "Narrow-sense heritability cannot exceed broad-sense heritability")
})
