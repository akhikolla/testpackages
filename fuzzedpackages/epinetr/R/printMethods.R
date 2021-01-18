#' Print function for population.
#'
#' Print a summary of the population object.
#'
#' This is an S3 method for printing a summary of a \code{Population}
#' object. Displayed are the initial parameters for the population
#' (i.e. population size, phenotypic variance, broad-sense
#' heritability, narrow-sense heritability and the SNPs used as
#' QTLs), followed by any current additive and epistatic variance in
#' the population.
#'
#' @param x a valid \code{Population} object
#' @param ... additional parameters (ignored)
#'
#' @author Dion Detterer, Paul Kwan, Cedric Gondro
#'
#' @examples
#' # Build a population
#' pop <- Population(
#'   popSize = 100, map = map100snp, QTL = 20,
#'   alleleFrequencies = runif(100), broadH2 = 0.9,
#'   narrowh2 = 0.5, traitVar = 40
#' )
#' pop <- addEffects(pop)
#' pop <- attachEpiNet(pop)
#' 
#' # Print the initial population
#' pop
#' 
#' # Run population in simulation
#' pop2 <- runSim(pop, generations = 150)
#' 
#' # Print the population following the simulation
#' pop2
#' @seealso \code{\link{Population}}, \code{\link{addEffects}},
#' \code{\link{attachEpiNet}}, \code{\link{runSim}}
#'
#' @export
print.Population <- function(x, ...) {
  testPop(x)

  geno <- (x$hap[[1]] + x$hap[[2]])[, x$qtl]

  if (!is.null(x$additive) && x$h2 > 0) {
    add <- (geno %*% t(t(x$additive))) + x$addOffset
  }

  if (!is.null(x$epiNet) && x$H2 > x$h2) {
    epi <- getEpi(x) * x$epiScale + x$epiOffset
  }

  cat("Population of size", x$popSize, "\n")
  cat("Specified initial trait variance:", x$VarP, "\n")
  cat("Initial broad-sense heritability:", x$H2, "\n")
  cat("Initial narrow-sense heritability:", x$h2, "\n")

  if (!is.null(x$additive) && x$h2 > 0) {
    cat(
      paste("Additive variance in population:", round(var(add), digits = 5)),
      "\n"
    )
  }

  if (!is.null(x$epiNet) && x$H2 > x$h2) {
    cat(paste("Epistatic variance in population:", round(var(epi),
      digits = 5
    )), "\n")
  }

  if (length(x$qtl) <= 100) {
    cat(
      "Using", length(x$alleleFreq), "SNPs with", length(x$qtl),
      "QTLs:\n"
    )

    for (i in x$qtl) cat(x$map$SNP[i], " ")

    cat("\n")
  }
}
