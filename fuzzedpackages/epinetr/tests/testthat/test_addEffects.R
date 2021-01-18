context("addEffects")
library(epinetr)

test_that("Sole additive effects give correct mean and variance", {
  pop <- Population(
    popSize = 10, map = map100snp, QTL = 20, alleleFrequencies = runif(100), broadH2 = 1,
    narrowh2 = 1, traitVar = 40
  )
  pop <- addEffects(pop)
  components <- getComponents(pop)
  expect_equal(mean(components$Additive), 0)
  expect_equal(var(components$Additive), 40)
  expect_equal(mean(components$Phenotype), 0)
  expect_equal(var(components$Phenotype), 40)

  pop <- addEffects(pop, distrib = runif)
  components <- getComponents(pop)
  expect_equal(mean(components$Additive), 0)
  expect_equal(var(components$Additive), 40)
  expect_equal(mean(components$Phenotype), 0)
  expect_equal(var(components$Phenotype), 40)

  pop <- addEffects(pop, effects = runif(20) * 1000)
  components <- getComponents(pop)
  expect_equal(mean(components$Additive), 0)
  expect_equal(var(components$Additive), 40)
  expect_equal(mean(components$Phenotype), 0)
  expect_equal(var(components$Phenotype), 40)

  pop <- Population(pop, traitVar = 20)
  pop <- addEffects(pop)
  components <- getComponents(pop)
  expect_equal(mean(components$Additive), 0)
  expect_equal(var(components$Additive), 20)
  expect_equal(mean(components$Phenotype), 0)
  expect_equal(var(components$Phenotype), 20)

  pop <- addEffects(pop, distrib = runif)
  components <- getComponents(pop)
  expect_equal(mean(components$Additive), 0)
  expect_equal(var(components$Additive), 20)
  expect_equal(mean(components$Phenotype), 0)
  expect_equal(var(components$Phenotype), 20)

  pop <- addEffects(pop, effects = runif(20) * 1000)
  components <- getComponents(pop)
  expect_equal(mean(components$Additive), 0)
  expect_equal(var(components$Additive), 20)
  expect_equal(mean(components$Phenotype), 0)
  expect_equal(var(components$Phenotype), 20)
})

test_that("Additive and environmental effects give the correct means and variances", {
  pop <- Population(
    popSize = 10, map = map100snp, QTL = 20, alleleFrequencies = runif(100), broadH2 = 0.7,
    narrowh2 = 0.7, traitVar = 40
  )

  pop <- addEffects(pop)
  components <- getComponents(pop)
  expect_equal(mean(components$Additive), 0)
  expect_equal(var(components$Additive), 28)
  expect_equal(mean(components$Phenotype), 0)
  # expect_equal(var(components$Phenotype), 40)

  pop <- addEffects(pop, distrib = runif)
  components <- getComponents(pop)
  expect_equal(mean(components$Additive), 0)
  expect_equal(var(components$Additive), 28)
  expect_equal(mean(components$Phenotype), 0)
  # expect_equal(var(components$Phenotype), 40)

  pop <- addEffects(pop, effects = runif(20) * 1000)
  components <- getComponents(pop)
  expect_equal(mean(components$Additive), 0)
  expect_equal(var(components$Additive), 28)
  expect_equal(mean(components$Phenotype), 0)
  # expect_equal(var(components$Phenotype), 40)
})

test_that("Changing narrow-sense heritability changes variances accordingly", {
  pop <- Population(
    popSize = 10, map = map100snp, QTL = 20, alleleFrequencies = runif(100), broadH2 = 0.7,
    narrowh2 = 0.7, traitVar = 40
  )
  pop <- addEffects(pop)
  pop <- Population(pop, narrowh2 = 0.75, broadH2 = 0.75)

  components <- getComponents(pop)
  expect_equal(mean(components$Additive), 0)
  expect_equal(var(components$Additive), 30)
  expect_equal(mean(components$Phenotype), 0)
  # expect_equal(var(components$Phenotype), 40)
})

test_that("Changing trait variance changes variances accordingly", {
  pop <- Population(
    popSize = 10, map = map100snp, QTL = 20, alleleFrequencies = runif(100), broadH2 = 0.5,
    narrowh2 = 0.5, traitVar = 40
  )
  pop <- addEffects(pop)
  pop <- Population(pop, traitVar = 30)

  components <- getComponents(pop)
  expect_equal(mean(components$Additive), 0)
  expect_equal(var(components$Additive), 15)
  expect_equal(mean(components$Phenotype), 0)
  # expect_equal(var(components$Phenotype), 30)
})

test_that("Changing narrow-sense heritability and trait variance changes variances accordingly", {
  pop <- Population(
    popSize = 10, map = map100snp, QTL = 20, alleleFrequencies = runif(100), broadH2 = 0.5,
    narrowh2 = 0.5, traitVar = 40
  )
  pop <- addEffects(pop)
  pop <- Population(pop, broadH2 = 0.2, narrowh2 = 0.2, traitVar = 50)

  components <- getComponents(pop)
  expect_equal(mean(components$Additive), 0)
  expect_equal(var(components$Additive), 10)
  expect_equal(mean(components$Phenotype), 0)
  # expect_equal(var(components$Phenotype), 50)
})

test_that("Removing additive effects works", {
  pop <- Population(
    popSize = 10, map = map100snp, QTL = 20, alleleFrequencies = runif(100), broadH2 = 0.5,
    narrowh2 = 0.5, traitVar = 40
  )
  pop <- addEffects(pop)
  pop <- Population(pop, broadH2 = 0, narrowh2 = 0, traitVar = 50)

  components <- getComponents(pop)
  expect_equal(mean(components$Additive), 0)
  expect_equal(var(components$Additive), 0)
  expect_equal(mean(components$Phenotype), 0)
  expect_equal(var(components$Phenotype), 50)
})

test_that("Additive effects are merely scaled", {
  pop <- Population(
    popSize = 10, map = map100snp, QTL = 20, alleleFrequencies = runif(100), broadH2 = 0.5,
    narrowh2 = 0.5, traitVar = 40
  )
  effects <- rnorm(20)
  pop <- addEffects(pop, effects = effects)
  reteffects <- getAddCoefs(pop)
  effects <- effects / reteffects
  expect_equal(abs(max(effects) - min(effects)) < 1e-10, TRUE)
})
