context("Population constructor")
library(epinetr)

test_that("New popSize overrides genotypes", {
  popSizeGiven <- nrow(geno100snp) - 5
  pop <- Population(
    popSize = popSizeGiven, map = map100snp, QTL = 20, genotypes = geno100snp, broadH2 = 0.9,
    narrowh2 = 0.6, traitVar = 40
  )
  popSize <- nrow(getPhased(pop))
  expect_equal(popSizeGiven, popSize)
})

test_that("New population size overrides old", {
  pop1 <- Population(
    popSize = 10, map = map100snp, QTL = 20, alleleFrequencies = runif(100), broadH2 = 0.9,
    narrowh2 = 0.6, traitVar = 40
  )
  pop2 <- Population(pop1, popSize = 8)
  popSize <- nrow(getPhased(pop2))
  expect_equal(popSize, 8)
})

test_that("Number of rows in genotypes overrides old popSize if literal is TRUE", {
  pop1 <- Population(
    popSize = 10, map = map100snp, QTL = 20, alleleFrequencies = runif(100), broadH2 = 0.9,
    narrowh2 = 0.6, traitVar = 40
  )
  pop2 <- Population(pop1, genotypes = geno100snp, literal = TRUE)
  popSize <- nrow(getPhased(pop2))
  expect_equal(popSize, nrow(geno100snp))
})

test_that("QTLs are selected correctly", {
  ids <- character(100)
  for (i in 0:9)
    for (j in 1:10)
      ids[i * 10 + j] <- paste(letters[i + 1], letters[j], sep = "")

  map1 <- map100snp
  map1$V1 <- ids

  qtls <- sample(ids, 20)

  pop <- Population(
    popSize = 2, map = map1, QTL = qtls, alleleFrequencies = runif(100), broadH2 = 0.9,
    narrowh2 = 0.6, traitVar = 40
  )

  retqtls <- getQTL(pop)$ID

  # Check that all qtls are in retqtls and vice versa
  expect_equal(all(qtls %in% retqtls) && all(retqtls %in% qtls), TRUE)

  # Check that indices are correct
  indices <- sort(which(ids %in% qtls))

  expect_equal(indices, getQTL(pop)$Index)
})

test_that("Correct number of QTLs are selected randomly", {
  pop <- Population(
    popSize = 2, map = map100snp, QTL = 20, alleleFrequencies = runif(100), broadH2 = 0.9,
    narrowh2 = 0.6, traitVar = 40
  )

  expect_equal(nrow(getQTL(pop)), 20)
})

test_that("No heritability gives correct environmental mean and variance", {
  pop <- Population(
    popSize = 10, map = map100snp, QTL = 20, alleleFrequencies = runif(100), broadH2 = 0,
    narrowh2 = 0, traitVar = 40
  )
  components <- getComponents(pop)
  expect_equal(nrow(components), 10)
  expect_equal(mean(components$Environmental), 0)
  expect_equal(var(components$Environmental), 40)
  expect_equal(mean(components$Phenotype), 0)
  expect_equal(var(components$Phenotype), 40)
})

test_that("Pedigree is the size of population", {
  pop <- Population(
    popSize = 10, map = map100snp, QTL = 20, alleleFrequencies = runif(100), broadH2 = 0,
    narrowh2 = 0, traitVar = 40
  )
  expect_equal(nrow(getPedigree(pop)), 10)
})
