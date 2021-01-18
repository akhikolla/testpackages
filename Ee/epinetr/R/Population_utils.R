# Add size parameter to population
addPopSize <- function(pop, pop2, popSize, genotypes, literal) {
  # Check that genotypes is a matrix, if given
  if (!is.null(genotypes) && !is.matrix(genotypes)) {
    stop("genotypes must be supplied via a matrix")
  }

  # Use explicit popSize if given
  if (!is.null(popSize)) {
    if (length(popSize) != 1 || !is.numeric(popSize) || popSize %% 1 != 0 ||
      popSize < 2) {
      stop("popSize must be a single integer greater than 1")
    }

    pop2$popSize <- popSize

    # Use literal genotypes to infer population size
  } else if (!is.null(genotypes) && literal) {
    pop2$popSize <- nrow(genotypes)

    # Use previous population's size
  } else if (!is.null(pop$popSize)) {
    pop2$popSize <- pop$popSize
  } else {
    stop("Population size not given")
  }

  return(pop2)
}


# Add map to population
addMap <- function(pop, pop2, map) {
  # Grab from previous population if not explicitly given
  if (is.null(map)) {
    if (is.null(pop$map)) {
      stop("Map not given")
    } else {
      map <- pop$map
    }
  }

  # Ensure map is a data frame
  if (!is.data.frame(map)) {
    stop("Map must be a data frame")
  }

  # Check number of columns
  if (ncol(map) != 3) {
    stop("Map must have 3 columns")
  }

  # Name columns
  names(map) <- c("SNP", "chr", "pos")

  # Positions of SNP must be positive integers
  if (!is.numeric(map$pos) || any(map$pos %% 1 != 0) || any(map$pos < 1)) {
    stop("Map in incorrect format")
  }

  # Sort map by position within each chromosome
  pop2$map <- sortMap(map)

  return(pop2)
}


# Add QTL to population
addQTL <- function(pop, pop2, QTL) {
  if (!is.null(QTL)) {
    if (length(QTL) == 1) {

      # Use a number of randomly selected SNP
      if (!is.numeric(QTL) || QTL %% 1 != 0 || QTL < 1 || QTL > nrow(pop2$map)) {
        stop("Number of QTLs must be a positive integer no greater than the number of SNP")
      }

      QTL <- sort(sample(1:nrow(pop2$map), QTL))
    } else if (length(QTL) > 1) {

      # Use a vector of SNP IDs
      if (!all(QTL %in% pop2$map$SNP)) {
        stop("Invalid QTLs")
      }
      QTL <- which(pop2$map$SNP %in% QTL)
      if (length(QTL) == 0) {
        stop("No QTLs found")
      }
    } else {
      stop("QTLs vector cannot have length 0")
    }
  } else if (!is.null(pop$qtl)) {

    # Use existing QTL
    if (!is.numeric(pop$qtl) || length(pop$qtl) > nrow(pop2$map) ||
      length(pop$qtl) < 1 || any(pop$qtl %% 1 != 0)) {
      stop("Invalid number of QTLs")
    }
    QTL <- pop$qtl
  } else {
    stop("QTLs not given")
  }

  pop2$qtl <- QTL

  return(pop2)
}


# Add haplotypes to the population
addHaplotypes <- function(pop, pop2, genotypes, literal, alleleFrequencies) {

  # Cannot specify both allele frequencies and genotypes
  if (!is.null(genotypes) && !is.null(alleleFrequencies)) {
    stop("Cannot use both genotypes and allele frequencies")
  }

  # If genotypes are given
  if (!is.null(genotypes)) {
    # Ensure genotypes is a matrix with an even number of columns
    if (!is.matrix(genotypes) || ncol(genotypes) %% 2 != 0) {
      stop("Genotypes must be in a matrix with an even number of columns")
    }

    # Ensure there is the right number of SNP
    if (ncol(genotypes) / 2 != nrow(pop2$map)) {
      stop("Discrepancy between number of SNP in map and genotypes")
    }

    # Convert to pairs of haplotypes
    hap <- geno2hap(genotypes)

    # Add to the population either literally or using derived allele
    # frequencies
    pop2 <- hapLogic(pop2, hap, literal)
  } else if (!is.null(alleleFrequencies)) {
    # Validate allele frequencies
    if (length(alleleFrequencies) != nrow(pop2$map)) {
      stop("Discrepancy between number of SNP in map and number of allele frequencies")
    }

    if (!is.numeric(alleleFrequencies) || any(alleleFrequencies > 1) ||
      any(alleleFrequencies < 0)) {
      stop("Invalid allele frequencies")
    }

    # Use allele frequencies to generate haplotypes
    pop2$hap <- generateHap(alleleFrequencies, pop2$popSize)
    pop2$alleleFreq <- calcAF(pop2$hap)
  } else if (!is.null(pop) && !is.null(pop$hap)) {
    # Validate current haplotypes
    if (!is.list(pop$hap)) {
      stop("Current haplotypes not in proper format")
    }

    for (i in 1:2) {
      if (!is.matrix(pop$hap[[i]]) || ncol(pop$hap[[i]]) != nrow(pop2$map)) {
        stop("Current haplotypes must be matrices with appropriate number of SNP")
      }
    }

    # Use existing haplotypes
    pop2 <- hapLogic(pop2, geno2hap(getPhased(pop)), literal)
  } else {
    stop("Genotypes or allele frequencies not specified.")
  }

  # Sort haplotypes and allele frequencies
  # pop2$hap[[1]] <- pop2$hap[[1]][, pop2$mapIndices]
  # pop2$hap[[2]] <- pop2$hap[[2]][, pop2$mapIndices]
  # pop2$alleleFreq <- pop2$alleleFreq[pop2$mapIndices]

  return(pop2)
}


# Add variances to the population
addVariances <- function(pop, pop2, traitVar, broadH2, narrowh2) {
  # Grab traitVar from previous population if not explicitly given
  if (is.null(traitVar)) {
    if (is.null(pop$VarP)) {
      # stop("Trait variance not given")
      traitVar <- 1
    } else {
      traitVar <- pop$VarP
    }
  }

  # Validate trait variance
  if (!is.numeric(traitVar) || length(traitVar) != 1 || traitVar <= 0) {
    stop("Invalid trait variance")
  }

  pop2$VarP <- traitVar

  # Add variances based on broad-sense heritability
  pop2 <- addBroadH2(pop, pop2, broadH2)

  # Add variances based on narrow-sense heritability
  pop2 <- addNarrowh2(pop, pop2, narrowh2)

  # Verify all are in agreement
  if (pop2$VarA > pop2$VarG) {
    stop("Narrow-sense heritability cannot exceed broad-sense heritability")
  }

  return(pop2)
}


# Add sexes to the population
addSex <- function(pop, pop2) {

  # if the old pop has sexes, and the sexes are the same as popsize, copy
  # them over
  if (!is.null(pop$isMale) && length(pop$isMale) == pop2$popSize) {
    pop2$isMale <- pop$isMale
  } else {

    # Generate sexes
    pop2$isMale <- sample(c(TRUE, FALSE), pop2$popSize, replace = TRUE)

    # Ensure there is at least one male and one female
    if (sum(pop2$isMale) == 0) {
      pop2$isMale[sample(pop2$popSize, 1)] <- TRUE
    } else if (sum(pop2$isMale) == pop2$popSize) {
      pop2$isMale[sample(pop2$popSize, 1)] <- FALSE
    }
  }

  return(pop2)
}


# Sort map by position within each chromosome
sortMap <- function(map) {

  # Get a vector of chromosome IDs
  chr <- unique(map$chr)

  newmap <- data.frame()

  for (i in chr) {
    # Find all SNPs for this chromosome
    index <- which(map$chr == i)

    # Store the submap for this chromosome
    tempmap <- map[index, ]

    # Sort the submap by position
    tempmap <- tempmap[sort(tempmap$pos, index.return = TRUE)$ix, ]

    # Append to the sorted map
    newmap <- rbind(newmap, tempmap)
  }

  return(newmap)
}


# Process a genotype matrix into two haplotype matrices
geno2hap <- function(geno) {
  # Reject is there are any NAs
  if (any(is.na(geno))) {
    stop("Cannot process genotype matrix with NAs")
  }

  # Transform matrix into one column per SNP
  geno <- matrix(geno, nrow = nrow(geno) * 2, ncol = ncol(geno) / 2)

  # For each column, code alleles as 0 and 1
  for (i in 1:ncol(geno)) {
    # Get allele codes
    alleles <- unique(geno[, i])

    # Ensure SNP is biallelic
    if (length(alleles) == 0 || length(alleles) > 2) {
      stop("SNP must be biallelic.")
    }

    # Code alleles as 0 and 1
    index1 <- (geno[, i] == alleles[1])
    if (length(alleles == 2)) {
      index2 <- (geno[, i] == alleles[2])
    }
    geno[index, i] <- 0
    if (length(alleles == 2)) {
      geno[index, i] <- 1
    }
  }

  # Convert to integers
  geno <- matrix(as.integer(geno), nrow = nrow(geno), ncol = ncol(geno))

  # Get allele frequencies
  # alleleFrequencies <- .colMeans(geno, m = nrow(geno), n = ncol(geno))

  # Return list of haplotypes
  hap <- list()
  hap[[1]] <- geno[1:(nrow(geno) / 2), ]
  hap[[2]] <- geno[(nrow(geno) / 2 + 1):nrow(geno), ]

  return(hap)
}


# Either use a list of haplotypes literally or generate based on allele
# frequencies
hapLogic <- function(pop, hap, literal) {
  af <- calcAF(hap)

  if (nrow(hap[[1]]) != pop$popSize || !literal) {
    hap <- generateHap(af, pop$popSize)
  }

  pop$hap <- hap
  pop$alleleFreq <- af

  return(pop)
}


# Generate haplotypes for the population
generateHap <- function(af, popSize) {
  nsnp <- length(af)

  # Empty matrix
  mtemp <- matrix(0, nrow = popSize * 2, ncol = nsnp)

  # Ensure we have minor allele frequencies
  # index <- af > 0.5
  # af[index] <- 1 - af[index]

  # Populate SNP
  for (i in 1:nsnp) mtemp[, i] <- sample(0:1, popSize * 2,
      replace = TRUE,
      prob = c(1 - af[i], af[i])
    )

  # Return haplotypes
  hap <- list()
  hap[[1]] <- mtemp[1:(popSize), ]
  hap[[2]] <- mtemp[(popSize + 1):(popSize * 2), ]

  return(hap)
}


# Add variances based on broad-sense heritability
addBroadH2 <- function(pop, pop2, broadH2) {
  if (is.null(broadH2)) {

    # Grab broad-sense heritability from previous population if not
    # explicitly given
    if (is.null(pop$H2)) {
      pop2$VarG <- pop2$VarP * 0.5
    } else {
      pop2$VarG <- pop2$VarP * pop$H2
    }
  } else {

    # Calculate genetic and environmental variance
    if (!is.numeric(broadH2) || length(broadH2) != 1 || broadH2 < 0 ||
      broadH2 > 1) {
      stop("Invalid broad-sense heritability")
    } else {
      pop2$VarG <- pop2$VarP * broadH2
    }
  }

  pop2$VarE <- pop2$VarP - pop2$VarG
  pop2$H2 <- pop2$VarG / pop2$VarP

  return(pop2)
}


# Add variances based on narrow-sense heritability
addNarrowh2 <- function(pop, pop2, narrowh2) {
  if (is.null(narrowh2)) {

    # Grab narrow-sense heritability from previous population if not
    # explicitly given
    if (is.null(pop$h2)) {
      pop2$VarA <- pop2$VarP * 0.3
    } else {
      pop2$VarA <- pop2$VarP * pop$h2
    }
  } else {

    # Calculate additive variance
    if (!is.numeric(narrowh2) || length(narrowh2) != 1 || narrowh2 <
      0 || narrowh2 > 1) {
      stop("Invalid narrow-sense heritability")
    } else {
      pop2$VarA <- pop2$VarP * narrowh2
    }
  }

  pop2$h2 <- pop2$VarA / pop2$VarP

  return(pop2)
}
