context("runSim")
library(epinetr)

test_that("Lack of necessary additive effects produces error", {
  pop <- Population(
    popSize = 10, map = map100snp, QTL = 20, alleleFrequencies = runif(100), broadH2 = 0.9,
    narrowh2 = 0.6, traitVar = 40
  )
  expect_error(runSim(pop), "Additive effects have not been attached")
  pop <- attachEpiNet(pop)
  expect_error(runSim(pop), "Additive effects have not been attached")

  pop <- Population(
    popSize = 10, map = map100snp, QTL = 20, alleleFrequencies = runif(100), broadH2 = 0.9,
    narrowh2 = 0.9, traitVar = 40
  )
  expect_error(runSim(pop), "Additive effects have not been attached")

  pop <- Population(
    popSize = 10, map = map100snp, QTL = 20, alleleFrequencies = runif(100), broadH2 = 1,
    narrowh2 = 1, traitVar = 40
  )
  expect_error(runSim(pop), "Additive effects have not been attached")
})

test_that("Lack of necessary epistatic effects produces error", {
  pop <- Population(
    popSize = 10, map = map100snp, QTL = 20, alleleFrequencies = runif(100), broadH2 = 0.9,
    narrowh2 = 0, traitVar = 40
  )
  expect_error(runSim(pop), "Epistatic effects have not been attached")

  pop <- Population(
    popSize = 10, map = map100snp, QTL = 20, alleleFrequencies = runif(100), broadH2 = 1,
    narrowh2 = 0, traitVar = 40
  )
  expect_error(runSim(pop), "Epistatic effects have not been attached")
})

test_that("Running with a single chromosome works", {
  map1 <- data.frame(snp = 1:100, chr = rep(1, 100), pos = sort(sample(100000, 100)))
  pop <- Population(
    popSize = 250, map = map1, QTL = 20, alleleFrequencies = runif(100), broadH2 = 0.9,
    narrowh2 = 0.5, traitVar = 40
  )
  pop <- addEffects(pop)
  pop <- attachEpiNet(pop)
  skip_on_cran()
  expect_error(runSim(pop, generations = 150), NA)
})
