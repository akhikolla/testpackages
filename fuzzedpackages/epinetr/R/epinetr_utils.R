# Test to see if object is of class 'Population'
testPop <- function(pop) {
  if (!is(pop, "Population")) {
    stop("pop must be of class Population")
  }
}


# Calculate scaling factor for epistasis
calcEpiScale <- function(pop) {
  epiVals <- getEpi(pop)

  pop$epiScale <- sqrt(pop$VarG - pop$VarA) / sd(epiVals)
  pop$epiOffset <- -mean(epiVals) * pop$epiScale

  return(pop)
}


# Update the pedigree data frame with new phenotypic components
updatePedigree <- function(pop) {
  # Evaluate all phenotypic components across population
  pheno <- getPheno(pop)

  # Ensure that epistasis is orthogonal to additive effects
  if (pop$H2 > pop$h2 && !is.null(pop$epiNet) && pop$h2 > 0 && !is.null(pop$additive)) {
    pop <- getSeeds(pop)
    pop$epiOffset <- 0
    pop$epiScale <- 1
    epi <- getEpi(pop)

    pop$epiScale <- sqrt(pop$VarG - pop$VarA) / sd(epi)
    epi <- epi * pop$epiScale
    pop$epiOffset <- -mean(epi)
    epi <- epi + pop$epiOffset
    pheno[, 2] <- epi
    pheno[, 4] <- rowSums(pheno[, 1:3])
  }

  if (pop$H2 < 1) {
    pheno[, 3] <- fitRnd(pheno[, 1] + pheno[, 2], normalisedrnorm, length(pheno[, 1]), sqrt(pop$VarE))
    pheno[, 4] <- rowSums(pheno[, 1:3])
  }

  # Store total phenotypic value for each member of the population
  pop$phenotype <- pheno[, 4]

  sex <- pop$isMale
  sex[pop$isMale] <- "M"
  sex[!pop$isMale] <- "F"

  pop$ped <- data.frame(ID = pop$ID, Sire = rep(0, pop$popSize), Dam = rep(
    0,
    pop$popSize
  ), Additive = pheno[, 1], Epistatic = pheno[, 2], Environmental = pheno
  [
    ,
    3
  ], Phenotype = pheno[, 4], Sex = sex)

  return(pop)
}


# Calculate allele frequencies from haplotypes
calcAF <- function(hap) {
  return(.colMeans(hap[[1]] + hap[[2]], m = nrow(hap[[1]]), n = ncol(hap[[1]])) / 2)
}


# Transform a set of haplotypes into a genotypes matrix
hap2geno <- function(hap) {
  geno <- rbind(hap[[1]], hap[[2]])
  geno <- matrix(geno, nrow = nrow(geno) / 2, ncol = ncol(geno) * 2)
  return(geno)
}


# Transform a genotypes matrix into a set of haplotypes
geno2hap <- function(geno) {
  hap <- list()

  geno <- matrix(geno, nrow = nrow(geno) * 2, ncol = ncol(geno) / 2)
  hap[[1]] <- geno[1:(nrow(geno) / 2), ]
  hap[[2]] <- geno[(nrow(geno) / 2 + 1):nrow(geno), ]

  return(hap)
}


# Save the state of the random number generator
# saveRNG <- function() {
#   if (exists(".Random.seed", .GlobalEnv)) {
#     oldseed <- .GlobalEnv$.Random.seed
#   } else {
#     oldseed <- NULL
#   }
#
#   return(oldseed)
# }


# Restore the state of the random number generator
# restoreRNG <- function(oldseed) {
#   if (!is.null(oldseed)) {
#     .GlobalEnv$.Random.seed <- oldseed
#   } else {
#     rm(".Random.seed", envir = .GlobalEnv)
#   }
# }


# Get the epistasis for each individual
getEpi <- function(pop, geno = NULL) {
  epi <- getEpistasis(pop, scale = FALSE, geno = geno)

  # Get dimension of indices
  dimensions <- dim(epi)

  # Get epistatic values per individual
  epi <- .rowSums(epi, dimensions[1], dimensions[2])

  return(epi)
}


# Select seeds to minimise covariance with additive effects
getSeeds <- function(pop) {
  geno <- pop$hap[[1]][, pop$qtl] + pop$hap[[2]][, pop$qtl]
  add <- (geno %*% t(t(pop$additive)))[, 1]
  n <- ncol(pop$epiNet$Incidence)
  retval <- GA::ga(type = "real-valued", fitness = selectSeeds, pop, geno, add, lower = rep(10^5, n), upper = rep(10^8, n), run = 15)
  pop$epiNet$Seeds <- floor(retval@solution[1, ])

  return(pop)
}


# Objective function for selecting seeds to minimise covariance
selectSeeds <- function(x, pop, geno, add) {
  x <- floor(x)
  pop$epiScale <- 1
  pop$epiOffset <- 0
  pop$epiNet$Seeds <- floor(x)
  epi <- getEpi(pop, geno)

  return(-abs(cov(add, epi)))
}


# Returns a set of normalised values with a mean of 0 and a standard
# deviation of sdev.
normalisedrnorm <- function(n, sdev) {
  env <- rnorm(n)
  env <- (env - mean(env)) / sd(env) * sdev

  return(env)
}


# Returns a random set of values taken from a normal distribution,
# with a mean of 0 and a standard deviation of sdev.  Covariance
# with a second vector, vec, will be minimised within a tolerance
# given, with iter.max iterations.
fitRnd <- function(vec, fun, ..., tolerance = 0.005, iter.max = 1000) {
  best <- fun(...)

  if (abs(cov(vec, best)) > 0) {
    for (i in 1:iter.max) {
      vec2 <- fun(...)

      if (abs(cov(vec, vec2)) < abs(cov(vec, best))) {
        best <- vec2
      }

      if (abs(cov(vec, best)) <= tolerance) {
        break
      }
    }
  }

  return(best)
}


#' Load epinetr genotype file.
#'
#' Load genotypes from a previous epinetr session.
#'
#' When outputting all genotypes during an epinetr simulation run, the genotypes
#' will be written to a serialised format unique to epinetr. The \code{loadGeno}
#' function will load these genotypes into memory as a single matrix object.
#'
#' @param filename the filename for the epinetr genotypes file.
#'
#' @return a numeric matrix holding the genotypes
#'
#' @author Dion Detterer, Paul Kwan, Cedric Gondro
#'
#' @export
#'
#' @examples
#' # Load genotype file
#' filename <- system.file("extdata", "geno.epi", package = "epinetr")
#' geno <- loadGeno(filename)
#'
#' # Use genotypes as basis for new population
#' pop <- Population(
#'   map = map100snp, QTL = 20, genotypes = geno,
#'   broadH2 = 0.8, narrowh2 = 0.6, traitVar = 40
#' )
#' @seealso \code{\link{Population}}
loadGeno <- function(filename) {
  geno <- getSerialMat(filename)

  if (is.null(geno)) {
    stop("Could not open file")
  }

  return(geno)
}
