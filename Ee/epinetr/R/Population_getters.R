#' Epistatic network retrieval.
#'
#' Get an epistatic network from a \code{Population} object.
#'
#' \code{getEpiNet()} is merely a function for retrieving an
#' epistatic network object. The common purpose is to plot the
#' network.
#'
#' @param pop An object of class \code{'Population'} which has an
#'   \code{EpiNet} object attached.
#'
#' @return An object of class \code{'EpiNet'} is returned.
#' @author Dion Detterer, Paul Kwan, Cedric Gondro
#' @export
#'
#' @examples
#' # Create population and attach an epistatic network
#' pop <- Population(
#'   popSize = 200, map = map100snp, QTL = 20,
#'   alleleFrequencies = runif(100),
#'   broadH2 = 0.9, narrowh2 = 0, traitVar = 40
#' )
#' pop <- attachEpiNet(pop)
#'
#' # Plot epistatic network
#' epiNet <- getEpiNet(pop)
#' plot(epiNet)
#' @seealso \code{\link{Population}}, \code{\link{plot.EpiNet}},
#' \code{\link{attachEpiNet}}
getEpiNet <- function(pop) {
  testPop(pop)

  pop$epiNet
}


#' Incidence matrix retrieval.
#'
#' Get an incidence matrix from a \code{Population}.
#'
#' \code{getIncMatrix()} retrieves the incidence matrix used in
#' epistatic interactions within the given \code{Population} object.
#' This is most useful for copying the network structure to a new
#' \code{Population} object.
#'
#' @param pop An object of class \code{'Population'} which has an
#'   \code{EpiNet} object attached.
#'
#' @return An incidence matrix representing the epistatic network
#'   within the given \code{Population} object.
#' @author Dion Detterer, Paul Kwan, Cedric Gondro
#' @export
#'
#' @examples
#' # Create population
#' pop <- Population(
#'   popSize = 200, map = map100snp, QTL = 20,
#'   alleleFrequencies = runif(100),
#'   broadH2 = 0.9, narrowh2 = 0, traitVar = 40
#' )
#'
#' # Attach random epistatic network and retrieve incidence matrix
#' pop <- attachEpiNet(pop)
#' inc <- getIncMatrix(pop)
#'
#' # Create second population
#' pop2 <- Population(
#'   popSize = 200, map = map100snp, QTL = 20,
#'   alleleFrequencies = runif(100),
#'   broadH2 = 0.8, narrowh2 = 0.6, traitVar = 40
#' )
#'
#' # Attach epistatic network to second population
#' # using incidence matrix from first
#' pop2 <- attachEpiNet(pop2, incmat = inc)
#' @seealso \code{\link{attachEpiNet}}
getIncMatrix <- function(pop) {
  testPop(pop)

  if (is.null(pop$epiNet)) {
    stop("There is no epistatic network attached")
  }

  mm <- pop$epiNet$Incidence
  mm[mm > 0] <- 1

  return(mm)
}


#' QTL retrieval.
#'
#' Retrieve the QTLs being used for this population.
#'
#' \code{getQTL} retrieves the IDs and indices of the SNPs being used
#' as QTLs for a given \code{Population} object.
#'
#' @param pop An object of class \code{'Population'}.
#'
#' @return A \code{data.frame} containing the IDs and indices of all
#' QTLs is returned.
#'
#' @author Dion Detterer, Paul Kwan, Cedric Gondro
#'
#' @export
#'
#' @examples
#' # Create population
#' pop <- Population(
#'   popSize = 200, map = map100snp, QTL = 20,
#'   alleleFrequencies = runif(100),
#'   broadH2 = 0.9, narrowh2 = 0.6, traitVar = 40
#' )
#'
#' # Get the SNP IDs of the QTLs
#' getQTL(pop)
#' @seealso \code{\link{Population}}
getQTL <- function(pop) {
  testPop(pop)

  foo <- pop$map$SNP[pop$qtl]
  bar <- pop$qtl

  return(data.frame(ID = foo, Index = bar))
}


#' Get population pedigree.
#'
#' Retrieve the pedigree of a \code{Population} object.
#'
#' \code{getPedigree()} can be used to retrieve the pedigree of a
#' population, including the phenotypic components of all individuals
#' in the pedigree.
#'
#' @param pop a valid object of class \code{Population}.
#'
#' @return A \code{data.frame} containing vectors for each
#' individual's ID, its sire and dam IDs, and its phenotypic
#' components, across the pedigree.
#' @author Dion Detterer, Paul Kwan, Cedric Gondro
#' @export
#'
#' @examples
#' \donttest{
#' # Construct a population with additive and epistatic effects
#' pop <- Population(
#'   popSize = 200, map = map100snp, QTL = 20,
#'   alleleFrequencies = runif(100),
#'   broadH2 = 0.9, narrowh2 = 0.6, traitVar = 40
#' )
#' pop <- addEffects(pop)
#' pop <- attachEpiNet(pop)
#'
#' # Run the simulator
#' pop2 <- runSim(pop, generations = 150)
#'
#' # Retrieve the population pedigree from the simulation
#' ped <- getPedigree(pop2)
#'
#' # Re-run the simulation using the same pedigree
#' pop3 <- runSim(pop, ped)
#' }
#' @seealso \code{\link{runSim}}
getPedigree <- function(pop) {
  testPop(pop)

  return(pop$ped)
}


#' Get allele frequencies.
#'
#' Get allele frequencies across a simulation run.
#'
#' Retrieves the allele frequencies in each generation across a
#' simulation run.
#'
#' @param pop a \code{Population} object as returned by a previous
#' simulation run.
#'
#' @return Returns a matrix of allele frequencies, one generation per
#' row.
#' @author Dion Detterer, Paul Kwan, Cedric Gondro
#' @export
#'
#' @examples
#' \donttest{
#' # Construct a population with additive and epistatic effects
#' pop <- Population(
#'   popSize = 200, map = map100snp, QTL = 20,
#'   alleleFrequencies = runif(100),
#'   broadH2 = 0.9, narrowh2 = 0.6, traitVar = 40
#' )
#' pop <- addEffects(pop)
#' pop <- attachEpiNet(pop)
#'
#' # Run the simulator
#' pop2 <- runSim(pop, generations = 150)
#'
#' af <- getAlleleFreqRun(pop2)
#' }
#' @seealso \code{\link{runSim}}
getAlleleFreqRun <- function(pop) {
  testPop(pop)

  return(pop$af)
}


#' Get population phased genotypes.
#'
#' Retrieves the current phased genotypes in the population.
#'
#' \code{getPhased} retrieves the current phased genotypes in the population,
#' returning a single matrix with one individual per row and two
#' columns per SNP.
#'
#' @param pop a valid \code{Population} object.
#'
#' @return Returns a phased genotypes matrix.
#' @author Dion Detterer, Paul Kwan, Cedric Gondro
#' @export
#'
#' @examples
#' # Construct a population
#' pop <- Population(
#'   popSize = 200, map = map100snp, QTL = 20,
#'   alleleFrequencies = runif(100),
#'   broadH2 = 0.9, narrowh2 = 0.6, traitVar = 40
#' )
#'
#' # Retrieve genotypes
#' geno <- getPhased(pop)
#' @seealso \code{\link{Population}}, \code{\link{getHaplo}}, \code{\link{getGeno}}
getPhased <- function(pop) {
  testPop(pop)

  hap <- pop$hap

  return(hap2geno(hap))
}


#' Get population unphased genotypes.
#'
#' Retrieves the current unphased genotypes in the population.
#'
#' \code{getGeno} retrieves the current unphased genotypes in the population,
#' returning a single matrix with one individual per row and one SNP per
#' column.
#'
#' @param pop a valid \code{Population} object.
#'
#' @return Returns an unphased genotypes matrix.
#' @author Dion Detterer, Paul Kwan, Cedric Gondro
#' @export
#'
#' @examples
#' # Construct a population
#' pop <- Population(
#'   popSize = 200, map = map100snp, QTL = 20,
#'   alleleFrequencies = runif(100),
#'   broadH2 = 0.9, narrowh2 = 0.6, traitVar = 40
#' )
#'
#' # Retrieve genotypes
#' geno <- getGeno(pop)
#' @seealso \code{\link{Population}}, \code{\link{getPhased}}, \code{\link{getHaplo}}
getGeno <- function(pop) {
  testPop(pop)

  hap <- pop$hap

  return(hap[[1]] + hap[[2]])
}


#' Retrieve haplotypes.
#'
#' Retrieve haplotypes from the population.
#'
#' Retrieves the haplotypes from both the population which, when
#' added together, form the genotypes.
#'
#' @param pop a valid \code{Population} object
#'
#' @return A list is returned with two elements, corresponding to the
#'   two haplotype matrices.
#'
#' @author Dion Detterer, Paul Kwan, Cedric Gondro
#'
#' @export
#'
#' @examples
#' # Construct a new population with additive effects
#' pop <- Population(
#'   popSize = 20, map = map100snp, QTL = 20,
#'   broadH2 = 0.4, narrowh2 = 0.4, traitVar = 40,
#'   alleleFrequencies = runif(100, 0.05, 0.5)
#' )
#' pop <- addEffects(pop)
#'
#' # Find the additive contribution to the individuals' phenotypes
#' hap <- getHaplo(pop)
#' hap <- (hap[[1]] + hap[[2]])[, getQTL(pop)$Index]
#' (hap %*% getAddCoefs(pop))[, 1] + getAddOffset(pop)
#' @seealso \code{\link{getAddCoefs}}, \code{\link{getAddOffset}}, \code{\link{getPhased}}, \code{\link{getGeno}}
getHaplo <- function(pop) {
  testPop(pop)

  return(pop$hap)
}


#' Get phenotypic components.
#'
#' Get pedigree and phenotypic data from current population.
#'
#' Retrieves the pedigree and phenotypic data components from all
#' individuals in the current population.
#'
#' @param pop a \code{Population} object
#'
#' @return Returns a \code{data.frame} giving the individual's ID,
#' the ID of the sire, the ID of the dam, the additive, epistatic and
#' environmental components of the phenotype and the overall
#' phenotypic value.
#' @author Dion Detterer, Paul Kwan, Cedric Gondro
#' @export
#'
#' @examples
#' # Construct a population with additive and epistatic effects
#' pop <- Population(
#'   popSize = 200, map = map100snp, QTL = 20,
#'   alleleFrequencies = runif(100),
#'   broadH2 = 0.9, narrowh2 = 0.6, traitVar = 40
#' )
#' pop <- addEffects(pop)
#' pop <- attachEpiNet(pop)
#'
#' # Retrieve phenotypic components from population
#' components <- getComponents(pop)
#' @seealso \code{\link{Population}}, \code{\link{addEffects}},
#'   \code{\link{attachEpiNet}}
getComponents <- function(pop) {
  testPop(pop)

  # Find indices of rows which match current IDs
  comp <- pop$ped[pop$ped$ID %in% pop$ID, ]

  # Sort the data frame according to current IDs
  indices <- match(pop$ID, comp$ID)
  comp <- comp[indices, ]

  return(comp)
}


#' Get additive coefficients.
#'
#' Retrieve additive coefficients from population.
#'
#' \code{getAddCoefs} retrieves the additive coefficients currently
#' in use in a \code{Population} object, assuming additive effects
#' have been attached.
#'
#' @param pop a \code{Population} object with additive effects
#'   attached
#'
#' @return \code{getAddCoefs} returns the additive coefficients
#'   currently in use by the population.
#' @author Dion Detterer, Paul Kwan, Cedric Gondro
#' @export
#'
#' @examples
#' # Construct a population with additive effects
#' pop <- Population(
#'   popSize = 200, map = map100snp, QTL = 20,
#'   alleleFrequencies = runif(100),
#'   broadH2 = 0.6, narrowh2 = 0.6, traitVar = 40
#' )
#' pop <- addEffects(pop)
#'
#' # Get additive coefficients
#' additive <- getAddCoefs(pop)
#' @seealso \code{\link{addEffects}}
getAddCoefs <- function(pop) {
  testPop(pop)

  return(pop$additive)
}


#' Retrieve interaction values.
#'
#' Retrieve  values for an interaction within an epistatic network.
#'
#' This function returns a \eqn{k}-dimensional array for a particular
#' interaction, where \eqn{k} is the order of interaction and the
#' array holds \eqn{3^k} entries . This means that a 5-way
#' interaction, for example, will return a 5-dimensional array
#' consisting of \eqn{3^5 = 243} entries.
#'
#' Within each dimension, the three indices (1, 2 and 3) correspond
#' to the homozygous genotype coded 0/0, the heterozygous genotype
#' and the homozygous genotype coded 1/1, respectively. Each entry
#' is drawn from a normal distribution. (Any offset needs to be
#' applied manually using \code{getEpiOffset}.)
#'
#' @param pop a valid population object
#'
#' @param n the interaction to return
#'
#' @author Dion Detterer, Paul Kwan, Cedric Gondro
#'
#' @examples
#' # Construct a new population
#' pop <- Population(
#'   popSize = 150, map = map100snp, QTL = 20,
#'   alleleFrequencies = runif(100), broadH2 = 0.7,
#'   narrowh2 = 0.45, traitVar = 40
#' )
#'
#' # Attach additive effects
#' pop <- addEffects(pop)
#'
#' # Attach a network of epistatic effects
#' pop <- attachEpiNet(pop)
#'
#' # Retrieve the possible values for the first two-way interaction
#' getInteraction(pop, 1)
#'
#' # Retrieve the value for the case where, in the fourth two-way
#' # interaction, the first QTL in the interaction is heterozygous
#' # and the second QTL in the interaction is the homozygous
#' # reference genotype.
#' getInteraction(pop, 4)[2, 1]
#'
#' # Retrieve the value for the case where, in the second two-way
#' # interaction, the first QTL in the interaction is the homozygous
#' # reference genotype and the second QTL in the interaction is the
#' # homozygous alternative genotype.
#' getInteraction(pop, 2)[1, 3]
#' @seealso \code{\link{attachEpiNet}}, \code{\link{getEpiOffset}}
#'
#' @export
getInteraction <- function(pop, n) {
  testPop(pop)

  nqtl <- sum(pop$epiNet$Incidence[, n] > 0)
  # temp <- saveRNG()
  # set.seed(pop$epiNet$Seeds[n])
  # vals <- rnorm(3^nqtl) * pop$epiScale # + pop$epiOffset / ncol(pop$epiNet$Incidence)
  # restoreRNG(temp)
  vals <- rng(3^nqtl, pop$epiNet$Seeds[n]) * pop$epiScale
  return(array(vals, dim = rep(3, nqtl)))
}


#' Calculate epistatic interactions.
#'
#' Calculate epistatic interactions for members of the population.
#'
#' This function calculates the values of each epistatic interaction
#' per individual, returning a matrix whose rows are the individuals
#' in the population and whose columns are the epistatic
#' interactions. The sums of these rows represent the total epistatic
#' contribution to the phenotype, once the epistatic offset is added.
#'
#' @param pop a valid \code{Population} object with epistatic effects
#'   attached
#' @param scale a boolean value indicating whether to scale the
#'   values to bring them in line with the desired initial variance
#' @param geno by default, the function uses the current genotypes
#'   in the population; alternatively, \code{geno} is a user-supplied
#'   set of genotypes (limited to QTLs) on which to base the
#'   calculation
#'
#' @author Dion Detterer, Paul Kwan, Cedric Gondro
#'
#' @examples
#' # Construct a new population with epistatic effects
#' pop <- Population(
#'   popSize = 20, map = map100snp, QTL = 20,
#'   broadH2 = 0.4, narrowh2 = 0, traitVar = 40,
#'   alleleFrequencies = runif(100, 0.05, 0.5)
#' )
#' pop <- attachEpiNet(pop)
#'
#' # Find the epistatic contribution to the individuals' phenotypes
#' rowSums(getEpistasis(pop)) + getEpiOffset(pop)
#'
#' # Compare with epistatic component from getComponents()
#' getComponents(pop)$Epistatic
#' @seealso \code{\link{getEpiOffset}}
#'
#' @export
getEpistasis <- function(pop, scale = TRUE, geno = NULL) {
  testPop(pop)

  # Get the QTL values
  if (is.null(geno)) {
    geno <- (pop$hap[[1]] + pop$hap[[2]])[, pop$qtl]
  }

  # Get the indices per individual per interaction
  indices <- (geno %*% pop$epiNet$Incidence) + 1

  # Get the seeds for each interaction
  seeds <- pop$epiNet$Seeds

  # Save the RNG
  # oldseed <- saveRNG()

  # Get epistatic values per individual per interaction
  for (i in 1:ncol(indices)) {
    # set.seed(seeds[i])
    # rnvals <- rnorm(max(indices[, i]))
    rnvals <- rng(max(indices[, i]), seeds[i])
    indices[, i] <- rnvals[indices[, i]]
  }

  # Restore the RNG
  # restoreRNG(oldseed)

  if (!is.null(pop$epiScale) && scale) {
    indices <- indices * pop$epiScale
  }

  return(indices)
}


#' Retrieve additive offset.
#'
#' Retrieve offset used for calculating additive component.
#'
#' In order for the initial population to have an additive
#' component with a mean of 0 for its phenotype, an offset is added,
#' and it remains fixed across generations. This function retrieves
#' that offset.
#'
#' @param pop A valid \code{Population} object with additive effects
#'   attached
#'
#' @return The additive offset is returned.
#' @author Dion Detterer, Paul Kwan, Cedric Gondro
#'
#' @export
#'
#' @examples
#' # Construct a new population with additive effects
#' pop <- Population(
#'   popSize = 20, map = map100snp, QTL = 20,
#'   broadH2 = 0.4, narrowh2 = 0.4, traitVar = 40,
#'   alleleFrequencies = runif(100, 0.05, 0.5)
#' )
#' pop <- addEffects(pop)
#'
#' # Find the additive contribution to the individuals' phenotypes
#' hap <- getHaplo(pop)
#' hap <- (hap[[1]] + hap[[2]])[, getQTL(pop)$Index]
#' (hap %*% getAddCoefs(pop))[, 1] + getAddOffset(pop)
#'
#' # Compare with additive component from getComponents()
#' getComponents(pop)$Additive
#' @seealso \code{\link{getAddCoefs}}, \code{\link{addEffects}}
getAddOffset <- function(pop) {
  testPop(pop)

  return(pop$addOffset)
}


#' Retrieve epistatic offset.
#'
#' Retrieve offset used for calculating epistatic component.
#'
#' In order for the initial population to have an epistatic
#' component with a mean of 0 for its phenotype, an offset is added,
#' and it remains fixed across generations. This function retrieves
#' that offset.
#'
#' @param pop A valid \code{Population} object with epistatic effects
#'   attached
#'
#' @author Dion Detterer, Paul Kwan, Cedric Gondro
#'
#' @export
#'
#' @examples
#' # Construct a new population with epistatic effects
#' pop <- Population(
#'   popSize = 20, map = map100snp, QTL = 20,
#'   broadH2 = 0.4, narrowh2 = 0, traitVar = 40,
#'   alleleFrequencies = runif(100, 0.05, 0.5)
#' )
#' pop <- attachEpiNet(pop)
#'
#' # Find the epistatic contribution to the individuals' phenotypes
#' rowSums(getEpistasis(pop)) + getEpiOffset(pop)
#'
#' # Compare with epistatic component from getComponents()
#' getComponents(pop)$Epistatic
#' @seealso \code{\link{getEpistasis}}
getEpiOffset <- function(pop) {
  testPop(pop)

  return(pop$epiOffset)
}


#' Get subpopulation
#'
#' Retrieve a subset of a \code{Population} without recalculating effects
#'
#' \code{getSubPop()} returns a new \code{Population} object using the
#' individuals with IDs specified by the vector \code{ID}.
#'
#' Any additive and epistatic effects will be copied as-is to the new
#' \code{Population} object, with heritability parameters recalculated.
#'
#' Any IDs given but not present will be discarded.
#'
#' @param pop a valid object of class \code{Population}.
#' @param ID a vector giving the IDs of individuals to include in the
#'   subset
#'
#' @return A new \code{Population} object containing the specified
#'   individuals is returned.
#'
#' @author Dion Detterer, Paul Kwan, Cedric Gondro
#' @export
#'
#' @examples
#' \donttest{
#' # Construct a population with additive and epistatic effects
#' pop <- Population(
#'   popSize = 2000, map = map100snp, QTL = 20,
#'   alleleFrequencies = runif(100)
#' )
#' pop <- addEffects(pop)
#' pop <- attachEpiNet(pop)
#'
#' # Run the simulator
#' pop2 <- runSim(pop, generations = 10)
#'
#' # Create a new subpopulation of 500 individuals
#' ID <- getComponents(pop2)$ID
#' ID <- sample(ID, 500)
#' pop3 <- getSubPop(pop2, ID)
#' }
#' @seealso \code{\link{Population}}, \code{\link{getComponents}}
getSubPop <- function(pop, ID) {
  pop2 <- list()
  class(pop2) <- "Population"

  # Discard invalid IDs
  ID <- ID[ID %in% pop$ID]

  # Get the indices of the IDs
  indices <- match(ID, pop$ID)

  components <- getComponents(pop)[indices, 1:7]

  pop2$popSize <- length(ID)

  if (pop2$popSize == 0) {
    stop("IDs not found in population")
  }

  pop2$map <- pop$map
  pop2$qtl <- pop$qtl
  pop2$hap <- list()
  pop2$hap[[1]] <- pop$hap[[1]][indices, ]
  pop2$hap[[2]] <- pop$hap[[2]][indices, ]
  pop2$alleleFreq <- pop$alleleFreq
  pop2$VarP <- var(components$Phenotype)
  pop2$VarG <- var(components$Additive + components$Epistatic)
  pop2$VarE <- pop$VarE
  pop2$H2 <- pop2$VarG / pop2$VarP
  pop2$VarA <- var(components$Additive)
  pop2$h2 <- pop2$VarA / pop2$VarP
  pop2$isMale <- pop$isMale[indices]
  pop2$ID <- pop$ID[indices]
  pop2$additive <- pop$additive
  pop2$addOffset <- pop$addOffset
  pop2$epiNet <- pop$epiNet
  pop2$epiOffset <- pop$epiOffset
  pop2$epiScale <- pop$epiScale
  pop2$phenotype <- pop$phenotype[indices]

  indices <- match(ID, pop$ped$ID)
  pop2$ped <- pop$ped[indices, 1:7]
  pop2$ped$Sire <- rep(0, pop2$popSize)
  pop2$ped$Dam <- rep(0, pop2$popSize)

  return(pop2)
}
