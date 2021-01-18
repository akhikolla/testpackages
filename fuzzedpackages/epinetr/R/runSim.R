#' Run simulation on population.
#'
#' Run a forward-time simulation on a \code{Population} object.
#'
#' \code{runSim} is the forward-time simulation engine of the
#' \code{epinetr} package. A \code{Population} object with necessary
#' additive and epistatic effects must be supplied; all other
#' arguments are optional, though either \code{pedigree} or
#' \code{generations} must be supplied.
#'
#' \code{pedigree} should be a \code{data.frame} where the first
#' three columns are the ID, sire ID and dam ID respectively. Sire
#' and dam IDs of 0 indicate that the individual is in the first
#' generation; each ID in the first generation should match an ID in
#' the given \code{Population} object. The pedigree will be sorted
#' into generations before running, where a 'generation' in this case
#' is defined as the set of individuals whose parents are both from a
#' previous generation. If a pedigree is supplied, all further
#' arguments (which pertain to selection) will be ignored.
#'
#' \code{generations} is the number of generations through which the
#' simulation will iterate. The supplied population represents the
#' first generation: the default value of 2 for this argument thus
#' means that the simulator will simply return the next generation.
#'
#' \code{selection} is a string specifying 'ranking' for linear
#' ranking selection; any other string is interpreted as 'random' for
#' random selection.
#'
#' Linear ranking selection mimics natural selection: if the
#' individuals in a population of size \eqn{n} are each given a rank
#' \eqn{r} based on descending order of phenotypic value (i.e. the
#' individual with the highest phenotypic value is given the rank
#' \ifelse{html}{\out{<i>r</i><sub>1</sub> = 1}}{{\eqn{r_1 = 1}}{r_1 = 1}} while the individual with the lowest phenotypic value
#' is given the rank \ifelse{html}{\out{<i>r<sub>n</sub></i> = <i>n</i>}}{{\eqn{r_n = n}}{r_n = n}}), the probability of an individual
#' \eqn{i} being selected for mating is given by:
#'
#' \ifelse{html}{
#' \out{
#' <i>P</i>(<i>i</i> is selected) = 2(<i>n</i> - <i>r<sub>i</sub></i> + 1) / <i>n</i>(<i>n</i> + 1)
#' }}{{
#' \deqn{P(i \textrm{ is selected}) = \frac{2(n - r_i + 1)}{n(n + 1)}}
#' }{
#' P(i is selected) = 2(n - r_i + 1) / n(n + 1)
#' }}
#'
#' Selection occurs by the population first being split into male and female
#' sub-populations. Next, if the round is outside any initial burn-in period,
#' each sub-population is truncated to a proportion of its original size per
#' the values of \code{truncSire} and \code{truncDam}, respectively.
#'
#' When linear ranking selection is used, females are exhaustively
#' sampled, without replacement, for each mating pair using their linear
#' ranking probabilities, as given above; males are sampled for each mating
#' pair using their linear ranking probabilities but with replacement, where
#' they are each only replaced a maximum number of times as specified by
#' \code{breedSire}. Random selection occurs in the same manner, but all
#' probabilities are uniform. During any initial \code{burnIn} period, random
#' selection is enforced.
#'
#' Each mating pair produces a number of full-sibling offspring by sampling
#' once from the litter-size probability mass function given by \code{litterDist}
#' (with the default guaranteeing two full-sibling offspring per mating pair).
#' The PMF is specified via a vector giving the probabilities for each
#' litter size, starting with a litter size of 0. For example,
#' \code{c(0.2, 0.0, 0.1, 0.4, 0.3)} gives a 20\% chance of a litter
#' size of 0, a 10\% chance of litter size of 2, a 40\% chance of a
#' litter size of 3, a 30\% chance of a litter size of 4 and a 0\%
#' chance of a litter size of 1 or greater than 4.
#'
#' Half-siblings occur when sires can mate more than once per round (as given by
#' \code{breedSire}) or when sires or dams survive beyond one round (as given by
#' \code{roundsSire} and \code{roundsDam}, respectively). It is important to note
#' that \code{roundsSire} and \code{roundsDam}, which specify the maximum number
#' of generations for males and females to survive, respectively, will be ignored
#' in the case where an insufficient number of offspring are produced to replace
#' the individuals who have nonetheless survived the maximum number of rounds: in
#' this case, younger individuals will be preserved in order to meet the
#' population size.
#'
#' \code{recombination} is a vector of recombination rates between
#' SNPs. The length of this vector should be equal to the number of
#' SNPs in the population's \code{map} minus the number of
#' chromosomes. The order of the chromosomes is as per the
#' \code{map}.
#'
#' \code{allGenoFileName} is the name of a file in which the
#' phased genotype for every individual will be stored. The output
#' is serialised and can be read using \code{loadGeno}. If the
#' \code{allGenoFileName} argument is not given, no genotypes
#' will be written to file.
#'
#' @param pop a valid \code{Population} object with all necessary
#'   additive and epistatic effects attached
#' @param pedigree an optional \code{data.frame} giving the pedigree
#'   to follow for selection
#' @param generations an optional integer giving the number of
#'   generations to iterate through in the simulation
#' @param selection an optional string specifying random (the default) or
#'   linear ranking selection
#' @param burnIn an optional integer giving the initial number of
#'   generations in which to use random selection without truncation,
#'   even when linear ranking selection or truncation is otherwise employed
#' @param truncSire an optional value giving the proportion of the
#'   males in the population with the highest phenotypic value to
#'   select within
#' @param truncDam an optional value giving the proportion of the
#'   females in the population with the highest phenotypic value to
#'   select within
#' @param roundsSire an optional integer giving the maximum number of
#'   generations for a male to survive; see details below
#' @param roundsDam an optional integer giving the maximum number of
#'   generations for a female to survive; see details below
#' @param litterDist an optional vector giving the probability mass
#'   function for the litter sizes, starting with a litter of 0
#' @param breedSire an optional integer indicating the maximum number
#'   of times that a sire can breed per generation
#' @param mutation an optional value giving the rate of mutation per
#'   SNP
#' @param recombination an optional vector giving the probabilities
#'   for recombination events between each SNP
#' @param allGenoFileName a string giving a file name,
#'   indicating that all genotypes will be outputted to the file during the run
#'
#' @return A new \code{Population} object is returned.
#'
#' @author Dion Detterer, Paul Kwan, Cedric Gondro
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Create population
#' pop <- Population(
#'   popSize = 200, map = map100snp, QTL = 20,
#'   alleleFrequencies = runif(100),
#'   broadH2 = 0.9, narrowh2 = 0.6, traitVar = 40
#' )
#'
#' # Attach additive effects using a normal distribution
#' pop <- addEffects(pop)
#'
#' # Attach epistatic effects
#' pop <- attachEpiNet(pop)
#'
#' # Run simulation for 150 generations
#' pop <- runSim(pop, generations = 150)
#'
#' # Display results
#' pop
#'
#' # Plot results
#' plot(pop)
#' }
#' @seealso \code{\link{Population}}, \code{\link{addEffects}},
#' \code{\link{attachEpiNet}}, \code{\link{print.Population}},
#' \code{\link{plot.Population}}, \code{\link{loadGeno}}
runSim <- function(pop, pedigree = NULL, generations = 2, selection = "random",
                   burnIn = 0, truncSire = 1, truncDam = 1, roundsSire = 1, roundsDam = 1,
                   litterDist = c(0, 0, 1), breedSire = 10, mutation = 10^-9, recombination = NULL,
                   allGenoFileName = NULL) {
  testPop(pop)

  if (pop$h2 > 0 && is.null(pop$additive)) {
    stop("Additive effects have not been attached")
  }

  if (pop$h2 < pop$H2 && is.null(pop$epiNet)) {
    stop("Epistatic effects have not been attached")
  }

  pedDropper <- FALSE

  # If a separate pedigree has been supplied...
  if (!is.null(pedigree)) {
    pedDropper <- TRUE

    message("Validating pedigree data...")
    pedigree <- prepPed(pedigree)
    # message("Done.\n")

    # Sort the pedigree
    message("Sorting pedigree data into generations...")
    foo <- sortPed(pedigree)
    pedigree <- pedigree[foo$sort, ]

    # Get number of matings per generation
    matings <- foo$matings
    message(paste("Found", length(matings), "generations."))

    # Ensure that all first-generation pedigree IDs are present in the
    # population
    if (!all(pedigree$ID[1:matings[1]] %in% pop$ID)) {
      stop("All pedigree IDs must be present in the population")
    }

    # Remove any individuals in the population not in the pedigree's first
    # generation
    cull <- which(!(pop$ID %in% pedigree$ID[1:matings[1]]))
    if (length(cull) != 0) {
      pop$ID <- pop$ID[-cull]
      pop$hap[[1]] <- pop$hap[[1]][-cull, ]
      pop$hap[[2]] <- pop$hap[[2]][-cull, ]
      pop$phenotype <- pop$phenotype[-cull]
      pop$isMale <- pop$isMale[-cull]
      pop$popSize <- matings[1]
      pop$alleleFreq <- calcAF(pop$hap)

      if (!is.null(pop$additive) && pop$h2 > 0) {
        pop <- addEffects(pop, effects = pop$additive)
      }

      if (!is.null(pop$epiNet) && pop$h2 < pop$H2) {
        pop <- calcEpiScale(pop)
      }

      pop$ped <- pedigree
      pop <- updatePedigree(pop)
    } else {
      # Copy phenotypic components to matching pedigree IDs
      components <- getComponents(pop)
      pedigree[, 4:7] <- components[
        match(pedigree$ID, components$ID),
        4:7
      ]
      pop$ped <- pedigree
    }

    # Zero-out all phenotypic data beyond the pedigree's first generation
    pop$ped[(matings[1] + 1):nrow(pop$ped), 4:7] <- 0

    # Add round for each individual
    pop$ped$Round <- rep(1:length(matings), matings)

    # Get number of generations
    pop$numGen <- length(matings)
    if (pop$numGen < 2) {
      stop("Number of generations must be greater than 1")
    }
  } else {
    # Validate the arguments and pass parameters onto population
    pop <- prepPop(
      pop, generations, selection, burnIn, truncSire,
      truncDam, roundsSire, roundsDam, litterDist, breedSire
    )

    # Counter for animal IDs
    IDcounter <- pop$popSize + 1

    # Store allele frequencies
    pop$af[1, ] <- calcAF(pop$hap)
  }

  pop$allGenotypes <- (!is.null(allGenoFileName) && is.character(allGenoFileName))

  # Validate mutation rate
  if (length(mutation) != 1 || !is.numeric(mutation) || mutation < 0 ||
    mutation >= 1) {
    stop("Mutation probability must be greater than or equal to 0 and less than 1")
  }
  pop$mutProb <- 10^-9

  # Get recombination probabilies or generate if not given
  if (is.null(recombination)) {
    numChr <- length(levels(as.factor(pop$map$chr)))
    recombination <- rnorm(nrow(pop$map) - numChr, mean = 2 / 3000, sd = 1 / 7000)
    recombination[recombination < 0] <- 0
  }
  pop$recProb <- getRecProb(recombination, pop$map)

  # Storage for phenotype summary
  pop$summaryData <- matrix(0, pop$numGen, 6)
  pop$summaryData[1, ] <- summary(pop$phenotype)

  # Output sire and dam haplotypes to file
  if (pop$allGenotypes) {
    serialMat(hap2geno(pop$hap), allGenoFileName, append = FALSE)
  }

  # Print relevant statistics
  message("Completed generation 1")
  message(paste("Mean phenotypic value:", mean(pop$phenotype)))
  message(paste("Maximum phenotypic value:", max(pop$phenotype), "\n"))

  # Generation loop
  for (i in 2:pop$numGen) {
    if (pedDropper) {
      # Remove any individuals who do not have offspring in future
      # generations
      cull <- cullIndices(pop, pop$ped, matings, i - 1)
      if (length(cull) > 0) {
        pop$hap[[1]] <- pop$hap[[1]][-cull, ]
        pop$hap[[2]] <- pop$hap[[2]][-cull, ]
        pop$phenotype <- pop$phenotype[-cull]
        pop$ID <- pop$ID[-cull]
        pop$isMale <- pop$isMale[-cull]
        pop$popSize <- length(pop$ID)
      }

      # First and last pedigree indices for this generation
      index1 <- sum(matings[1:(i - 1)]) + 1
      index2 <- sum(matings[1:i])

      # Get indices of mating pairs
      males <- match(pop$ped$Sire[index1:index2], pop$ID)
      females <- match(pop$ped$Dam[index1:index2], pop$ID)
    } else {
      # Age population
      pop$roundsCount <- pop$roundsCount - 1

      # Select the sires and dams according to selection criteria

      # Do not use truncation during burn-in period
      truncm <- pop$truncMale
      truncf <- pop$truncFemale
      if (i <= pop$burnIn) {
        truncm <- 1
        truncf <- 1
      }

      # How many males should we select from?
      numMales <- round(sum(pop$isMale) * truncm)

      # How many females should we select from?
      numFemales <- round(sum(!pop$isMale) * truncf)

      if (numFemales < 2) {
        i <- i - 1
        warning("Exhausted female population")
        break
      }

      if (numMales == 0) {
        i <- i - 1
        warning("Exhausted male population")
        break
      }

      # Get the indices of the top males and females
      sorted <- sort(pop$phenotype, decreasing = TRUE, index.return = TRUE)$ix
      males <- sorted[pop$isMale[sorted]][1:numMales]
      females <- sorted[!pop$isMale[sorted]][1:numFemales]

      # Select based on phenotype ranking
      if (pop$ranking & (i > pop$burnIn)) {
        # Select sires for each pair within breeding count restriction
        males <- getMales(males, numFemales, rep(
          pop$sireBreed,
          pop$popSize
        ))

        # Select dams for each pair
        females <- sample(females, length(males),
          prob = ranking(numFemales),
          replace = FALSE
        )

        # Select randomly
      } else {
        # Randomly select within male population
        # males <- sample(which(pop$isMale), numMales, replace = FALSE)

        # Select sires for each pair within breeding count restriction
        males <- getMales(males, numFemales, rep(
          pop$sireBreed,
          pop$popSize
        ), random = TRUE)

        # Randomly select dams for sires
        females <- sample(which(!pop$isMale), length(males), replace = FALSE)
      }

      # Determine how many offspring for each pairing; remove pairings where
      # no offspring result
      valid <- (pop$litterDist != 0)
      outcomes <- (0:(length(pop$litterDist) - 1))[valid]
      probs <- pop$litterDist[valid]
      if (max(probs) == 1) {
        numOffspring <- rep(outcomes, length(females))
      } else {
        numOffspring <- sample(outcomes, length(females),
          prob = probs,
          replace = TRUE
        )
      }
      offindex <- (numOffspring > 0)
      males <- males[offindex]
      females <- females[offindex]
      numOffspring <- numOffspring[offindex]

      # Repeat pairs for each offspring they produce
      males <- rep(males, numOffspring)
      females <- rep(females, numOffspring)

      # Create offspring pedigree data frame
      zeroes <- rep(0, length(males))
      offped <- data.frame(
        ID = zeroes, Sire = zeroes, Dam = zeroes,
        Additive = zeroes, Epistatic = zeroes, Environmental = zeroes,
        Phenotype = zeroes, Sex = sample(c("M", "F"), length(males), replace = TRUE), Round = zeroes + i,
        stringsAsFactors = FALSE
      )
    }

    # Allocate space for this generation
    offspring <- list()
    offspring[[1]] <- matrix(0, length(males), nrow(pop$map))
    offspring[[2]] <- offspring[[1]]

    # Generate offspring
    for (j in 1:length(males)) {
      sire <- males[j]
      dam <- females[j]

      # Get the sire and dam haplotypes
      sireHap1 <- pop$hap[[1]][sire, ]
      sireHap2 <- pop$hap[[2]][sire, ]
      damHap1 <- pop$hap[[1]][dam, ]
      damHap2 <- pop$hap[[2]][dam, ]

      # Recombine sire and dam haplotypes
      offspring[[1]][j, ] <- recombine(sireHap1, sireHap2, pop$recProb)
      offspring[[2]][j, ] <- recombine(damHap1, damHap2, pop$recProb)

      # Mutate haplotypes
      mut <- (runif(length(sireHap1)) <= pop$mutProb)
      offspring[[1]][j, mut] <- abs(offspring[[1]][j, mut] - 1)
      mut <- (runif(length(sireHap1)) <= pop$mutProb)
      offspring[[2]][j, mut] <- abs(offspring[[2]][j, mut] - 1)

      # Update offspring pedigree
      if (!pedDropper) {
        offped[j, 1:3] <- c(IDcounter, pop$ID[sire], pop$ID[dam])
        IDcounter <- IDcounter + 1
      }
    }

    # Output sire and dam haplotypes
    if (pop$allGenotypes) {
      serialMat(hap2geno(list(offspring[[1]], offspring[[2]])), allGenoFileName, append = TRUE)
    }

    # Evaluate this generation
    pheno <- getPheno(pop, offspring[[1]] + offspring[[2]])

    # Append offspring to population
    pop$hap[[1]] <- rbind(pop$hap[[1]], offspring[[1]])
    pop$hap[[2]] <- rbind(pop$hap[[2]], offspring[[2]])

    # Append phenotype to pedigree and ID to population
    if (pedDropper) {
      pop$ped[index1:index2, 4:7] <- pheno
      pop$ID <- c(pop$ID, pop$ped[index1:index2, 1])
      pop$isMale <- c(pop$isMale, (pop$ped$Sex[index1:index2] == "M"))

      # Update population size
      pop$popSize <- length(pop$ID)
    } else {
      offped[, 4:7] <- pheno
      pop$ID <- c(pop$ID, offped$ID)
      sexOffspring <- (offped$Sex == "M")
      pop$isMale <- c(pop$isMale, sexOffspring)
    }

    # Append phenotype to population
    pop$phenotype <- c(pop$phenotype, pheno[, 4])

    # Append sex to population
    # sexOffspring <- sample(c(TRUE, FALSE), length(males), replace = TRUE)
    # pop$isMale <- c(pop$isMale, sexOffspring)

    if (!pedDropper) {
      # Append rounds to survive
      roundsCount <- rep(pop$damRounds, length(males))
      roundsCount[sexOffspring] <- pop$sireRounds
      pop$roundsCount <- c(pop$roundsCount, roundsCount)

      pedindex1 <- which(pop$ped$ID == 0)

      # Expand pedigree dataframe as needed
      if (length(pedindex1) == 0 || length(pop$ped$ID) - pedindex1[1] +
        1 < nrow(offped)) {
        zeroes <- rep(0, ceiling(pop$expectedOffspring))
        foo <- data.frame(
          ID = zeroes, Sire = zeroes, Dam = zeroes,
          Additive = zeroes, Environmental = zeroes, Epistatic = zeroes,
          Phenotype = zeroes, Sex = rep("X", ceiling(pop$expectedOffspring)), Round = zeroes
        )
        pop$ped <- rbind(pop$ped, foo)
      }

      # Store pedigree data for generated offspring
      pedindex1 <- which(pop$ped$ID == 0)[1]
      pedindex2 <- pedindex1 + nrow(offped) - 1
      pop$ped[pedindex1:pedindex2, ] <- offped

      # Get the indices of individuals to be culled due to age
      oldIndex <- pop$roundsCount < 1
      cull <- which(oldIndex)

      # If there aren't enough offspring to make up the difference,
      # we'll need to truncate the number culled
      if (length(cull) > length(males)) {
        cull <- cull[sort(pop$roundsCount[oldIndex], index.return = TRUE)$ix]
        cull <- cull[1:length(males)]
      }

      # Perform cull
      if (length(cull) > 0) {
        pop$hap[[1]] <- pop$hap[[1]][-cull, ]
        pop$hap[[2]] <- pop$hap[[2]][-cull, ]
        pop$ID <- pop$ID[-cull]
        pop$phenotype <- pop$phenotype[-cull]
        pop$isMale <- pop$isMale[-cull]
        pop$roundsCount <- pop$roundsCount[-cull]
      }

      # If population is now too large, further cull individuals accordingly
      if (length(pop$isMale) > pop$popSize) {
        # If we're within the burn-in period or selection is random, remove
        # individuals randomly, else remove those with lowest phenotype
        if (i <= pop$burnIn || !pop$ranking) {
          index <- sample(1:length(pop$isMale), pop$popSize, replace = FALSE)
        } else {
          index <- sort(pop$phenotype, decreasing = TRUE, index.return = TRUE)$ix[1:pop$popSize]
        }

        pop$hap[[1]] <- pop$hap[[1]][index, ]
        pop$hap[[2]] <- pop$hap[[2]][index, ]
        pop$ID <- pop$ID[index]
        pop$phenotype <- pop$phenotype[index]
        pop$isMale <- pop$isMale[index]
        pop$roundsCount <- pop$roundsCount[index]
      }

      # Store allele frequencies
      pop$af[i, ] <- calcAF(pop$hap)
    }

    # Store phenotype summary
    pop$summaryData[i, ] <- summary(pop$phenotype)

    # Print relevant statistics
    message(paste("Completed generation", i))
    message(paste("Mean phenotypic value:", mean(pop$phenotype)))
    message(paste("Maximum phenotypic value:", max(pop$phenotype), "\n"))
  }

  # Trim pedigree data
  if (min(pop$ped$ID) == 0) {
    index <- which(pop$ped$ID == 0)[1] - 1
    pop$ped <- pop$ped[1:index, ]
  }

  # Trim summary data
  pop$summaryData <- pop$summaryData[1:i, ]

  # Trim allele frequency data
  pop$af <- pop$af[1:i, ]

  pop$hasRun <- TRUE

  return(pop)
}
