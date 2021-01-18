# Select from a set of males, taking into account how many times they
# can breed per round
getMales <- function(males, numPairs, breedCount, random = FALSE) {
  maleIndex <- list()

  # Remove any males whose breed count is less than one
  index <- which(breedCount < 1)
  if (length(index) > 0) {
    males <- males[-which(males %in% index)]
  }

  for (i in 1:numPairs) {
    # Terminate if no males remain
    if (length(males) < 1) {
      break
    }

    # If one male remains, select it
    if (length(males) == 1) {
      maleIndex[[i]] <- males
    } else {

      # Select a male
      if (random) {
        maleIndex[[i]] <- sample(males, 1)
      } else {
        maleIndex[[i]] <- sample(males, 1, prob = ranking(length(males)))
      }
    }

    # Decrease this male's breed count and remove if necessary
    breedCount[maleIndex[[i]]] <- breedCount[maleIndex[[i]]] - 1
    if (breedCount[maleIndex[[i]]] < 1) {
      males <- males[males != maleIndex[[i]]]
    }
  }

  return(unlist(maleIndex))
}


# Given a pair of haplotypes and a set of recombination probabilities
# for each chromosome, produces a new haplotype
recombine <- function(hap1, hap2, recProb) {
  recPoints <- (runif(length(recProb)) < recProb)
  mask <- recMask(recPoints)
  hap1[mask] <- hap2[mask]

  return(hap1)
}


# Linear ranking function
ranking <- function(n) {
  return(2 * (n - 1:n + 1) / (n * (n + 1)))
}


prepPed <- function(pedigree) {
  # Ensure pedigree is a data frame with at least three columns
  if (!is.data.frame(pedigree) || ncol(pedigree) < 3) {
    stop("Pedigree must be a data frame with at least three columns")
  }

  # Check that all IDs are non-negative integers
  if (!(is.numeric(pedigree[, 1])) || !(is.numeric(pedigree[, 2])) ||
    !(is.numeric(pedigree[, 3])) || any(pedigree[, 1:3] %% 1 != 0)) {
    stop("Pedigree IDs must be integers.")
  }

  # Check that no offspring has only one parent
  if (sum(xor(pedigree[, 2] > 0, pedigree[, 3] > 0)) > 0) {
    stop("Individuals must have two parents or be from initial generation.")
  }

  # Ensure that the same ID does not appear more than once in the same
  # row
  for (i in 1:nrow(pedigree)) if (anyDuplicated(pedigree[i, 1:3]) &&
      pedigree[i, 2] != 0) {
      stop("No IDs can be repeated within the same row.")
    }

  # Ensure that no animal appears as both a sire and a dam
  sire <- pedigree[, 2]
  sire <- sire[sire != 0]
  dam <- pedigree[, 3]
  dam <- dam[dam != 0]
  if (any(sire %in% dam) || any(dam %in% sire)) {
    stop("No individual can be both a sire and a dam")
  }

  # Tag sires as male and dams as female
  males <- unique(match(sire, pedigree[, 1]))
  females <- unique(match(dam, pedigree[, 1]))
  sex <- rep("X", nrow(pedigree))
  sex[males] <- "M"
  sex[females] <- "F"
  sexless <- which(sex == "X")
  sex[sexless] <- sample(c("M", "F"), length(sexless), replace = TRUE)

  pedigree <- pedigree[, 1:3]
  names(pedigree) <- c("ID", "Sire", "Dam")
  pedigree$Additive <- rep(0, length(pedigree$ID))
  pedigree$Epistatic <- rep(0, length(pedigree$ID))
  pedigree$Environmental <- rep(0, length(pedigree$ID))
  pedigree$Phenotype <- rep(0, length(pedigree$ID))
  pedigree$Sex <- sex

  return(pedigree)
}


prepPop <- function(pop, generations, selection, burnIn, truncSire, truncDam,
                    roundsSire, roundsDam, litterDist, breedSire) {

  # Validate number of generations
  if (!is.numeric(generations) || length(generations) != 1 || generations %% 1 !=
    0 || generations < 2) {
    stop("Number of generations must be an integer greater than 1")
  }
  pop$numGen <- generations

  # Validate selection
  pop$ranking <- (tolower(trimws(selection)) == "ranking")

  # Validate burn-in period
  if (length(burnIn) != 1 || !is.numeric(burnIn) || burnIn %% 1 != 0 ||
    burnIn < 0 || burnIn > pop$numGen) {
    stop("Burn-in length must be an integer between 0 and the number of generations to run")
  }
  pop$burnIn <- burnIn

  # Validate truncation rates
  if (length(truncSire) != 1 || !is.numeric(truncDam) || truncSire >
    1 || truncSire <= 0) {
    stop("Sire truncation rate incorrectly specified")
  }
  pop$truncMale <- truncSire

  if (length(truncDam) != 1 || !is.numeric(truncDam) || truncDam > 1 ||
    truncDam <= 0) {
    stop("Dam truncation rate incorrectly specified")
  }
  pop$truncFemale <- truncDam

  # Validate rounds per sire and dam
  if (length(roundsSire) != 1 || length(roundsDam) != 1 || !is.numeric(roundsSire) ||
    !is.numeric(roundsDam) || roundsSire %% 1 != 0 || roundsDam %% 1 !=
    0 || roundsSire < 1 || roundsDam < 1 || roundsSire > pop$numGen ||
    roundsDam > pop$numGen) {
    stop("Sire and dam rounds incorrectly specified")
  }
  pop$sireRounds <- roundsSire
  pop$damRounds <- roundsDam

  # Validate litter distribution
  if (length(litterDist) < 2 || !is.numeric(litterDist) || any(litterDist <
    0)) {
    stop("Invalid litter distribution")
  }
  pop$litterDist <- litterDist / sum(litterDist)

  # Validate number of times a sire can breed per round
  if (length(breedSire) != 1 || !is.numeric(breedSire) || breedSire <
    1) {
    stop("Sires must each have a maximum breeding rate of at least once per round")
  }
  pop$sireBreed <- breedSire

  # Set counter for each animal
  pop$roundsCount <- rep(pop$damRounds, pop$popSize)
  pop$roundsCount[pop$isMale] <- pop$sireRounds

  # Calculate expected number of offspring
  pop$expectedOffspring <- (pop$numGen - 1) * 0.5 * pop$truncFemale *
    pop$popSize * sum(pop$litterDist * 0:(length(pop$litterDist) -
      1))

  # Expand pedigree data frame
  zeroes <- rep(0, pop$popSize + pop$expectedOffspring)
  tempped <- data.frame(
    ID = zeroes, Sire = zeroes, Dam = zeroes, Additive = zeroes,
    Epistatic = zeroes, Environmental = zeroes, Phenotype = zeroes, Sex = as.character(rep("X", pop$popSize + pop$expectedOffspring)), Round = zeroes,
    stringsAsFactors = FALSE
  )
  components <- getComponents(pop)
  components$Sex <- as.character(components$Sex)
  tempped[1:pop$popSize, ] <- components
  tempped$Round <- zeroes
  pop$ped <- tempped

  pop$ped$Round[1:pop$popSize] <- 1

  # Storage for allele frequencies
  pop$af <- matrix(0, pop$numGen, length(pop$alleleFreq))

  return(pop)
}


sortPed <- function(ped) {
  # Get the indices of initial individuals
  indices <- which(ped[, 2] == 0)

  # Size of each generation
  matings <- length(indices)

  # Get indices of individuals who have both parents already found in
  # previous rounds
  indices2 <- which(ped[, 2] %in% ped[indices, 1] & ped[, 3] %in% ped[
    indices,
    1
  ])
  indices2 <- indices2[!(indices2 %in% indices)]
  while (length(indices2) > 0) {
    indices <- c(indices, indices2)
    matings <- c(matings, length(indices2))
    indices2 <- which(ped[, 2] %in% ped[indices, 1] & ped[, 3] %in%
      ped[indices, 1])
    indices2 <- indices2[!(indices2 %in% indices)]
  }

  retval <- list()

  # Return sorted indices and size of each generation
  retval$sort <- indices
  retval$matings <- matings

  return(retval)
}


cullIndices <- function(pop, ped, matings, gen) {
  index1 <- sum(matings[1:gen]) + 1
  index2 <- nrow(ped)
  parents <- unique(c(ped$Sire[index1:index2], ped$Dam[index1:index2]))
  cull <- unique(which(!(pop$ID %in% parents)))

  return(cull)
}


# Derives a data structure for recombination probabilities from a
# probability vector and a map
getRecProb <- function(rp, map) {
  chr <- unique(map$chr)

  if (length(rp) != nrow(map) - length(chr)) {
    stop("Number of recombination probabilities must match number of possible recombination points")
  }

  if (any(rp < 0 | rp > 1)) {
    stop("Recombination probabilities must be between 0 and 1")
  }

  if (length(chr) <= 0) {
    stop("There must be at least 1 chromosome")
  }

  index <- 1
  index2 <- 1

  recProb <- numeric(nrow(map))

  for (i in chr) {
    recProb[index] <- 0.5
    numSNP <- sum(map$chr == i)
    recProb[(index + 1):(index + numSNP - 1)] <- rp[index2:(index2 +
      numSNP - 2)]
    index <- index + numSNP
    index2 <- index2 + numSNP - 1
  }

  return(recProb)
}


# Create a numbered filename for use with outputting all haplotypes
getFilename <- function(pop, generation) {
  # Get width of number
  numwidth <- nchar(paste(pop$numGen))

  # Get formatted number
  filenumber <- formatC(generation, width = numwidth, format = "d", flag = "0")

  filename <- tools::file_path_sans_ext(pop$hapout)
  fileext <- tools::file_ext(pop$hapout)

  filename <- paste(filename, filenumber, sep = "-")

  return(paste(filename, fileext, sep = "."))
}


# Fitness function
getPheno <- function(pop, geno = NULL) {
  if (is.null(geno)) {
    geno <- (pop$hap[[1]] + pop$hap[[2]])
  }

  geno <- geno[, pop$qtl]

  dimensions <- dim(geno)

  zeros <- rep(0, dimensions[1])

  if (pop$h2 > 0 && !is.null(pop$additive)) {
    add <- (geno %*% t(t(pop$additive))) + pop$addOffset
  } else {
    add <- zeros
  }

  if (pop$h2 < pop$H2 && !is.null(pop$epiNet)) {
    epi <- getEpi(pop, geno = geno) * pop$epiScale + pop$epiOffset
  } else {
    epi <- zeros
  }

  if (pop$H2 < 1) {
    env <- rnorm(nrow(geno), sd = sqrt(pop$VarE))
  } else {
    env <- zeros
  }

  tot <- add + epi + env

  return(matrix(c(add, epi, env, tot), length(add), 4))
}
