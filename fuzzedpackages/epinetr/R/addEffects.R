#' Attach additive effects to population.
#'
#' Attach additive effects to a \code{Population} object.
#'
#' \code{addEffects()} is a function for attaching additive effects to a given
#' population, ensuring that the initial additive variance is as given in the
#' population parameters.
#'
#' If additive effect coefficients are directly supplied via the \code{effects}
#' vector, these may be scaled in order to comply with the initial additive
#' variance.
#'
#' @param pop an object of class \code{Population}.
#' @param effects an optional vector of additive effect coefficients.
#' @param distrib an optional random number generator function for a
#'   distribution, defaulting to \code{\link{rnorm}}.
#'
#' @return A copy of the supplied \code{Population} object is returned, with
#'   additive effects attached.
#'
#' @author Dion Detterer, Paul Kwan, Cedric Gondro
#'
#' @examples
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
#' # Attach additive effects using a uniform distribution
#' pop2 <- addEffects(pop, distrib = runif)
#'
#' # Attach additive effects using a vector of coefficients
#' effects <- c(
#'   1.2, 1.5, -0.3, -1.4, 0.8,
#'   2.4, 0.2, -0.8, -0.4, 0.8,
#'   -0.2, -1.4, 1.4, 0.2, -0.9,
#'   0.4, -0.8, 0.0, -1.1, -1.3
#' )
#' pop3 <- addEffects(pop, effects = effects)
#'
#' # Print first population
#' pop
#'
#' # Print second population
#' pop2
#'
#' # Print third population
#' pop3
#' @export
#'
#' @seealso \code{\link{Population}}, \code{\link{attachEpiNet}}
addEffects <- function(pop, effects = NULL, distrib = rnorm) {
  testPop(pop)

  if (pop$h2 == 0) {
    stop("No additive variance in population")
  }

  # Initial coefficients
  if (!is.null(effects)) {
    if (length(effects) != length(pop$qtl) || !is.numeric(effects)) {
      stop("Additive effects vector does not match QTL")
    } else {
      coef <- effects
    }
  } else {
    coef <- distrib(length(pop$qtl))
  }

  if (!is.numeric(coef) || length(coef) != length(pop$qtl)) {
    stop("distrib is not a valid distribution")
  }

  # Get additive components of phenotypic trait within population
  addComp <- as.vector((pop$hap[[1]] + pop$hap[[2]])[, pop$qtl] %*% coef)

  # Scale the coefficients to match user additive variance
  pop$additive <- coef * sqrt(pop$VarA / var(addComp))

  # Add an offset to set the initial mean to 0
  addComp <- as.vector((pop$hap[[1]] + pop$hap[[2]])[, pop$qtl] %*% pop$additive)
  pop$addOffset <- -mean(addComp)

  # Update pedigree data frame if epiNet is present or unneeded (but not
  # both)
  if ((!is.null(pop$epiNet) || pop$h2 == pop$H2) && (pop$h2 < pop$H2 ||
    is.null(pop$epiNet))) {
    pop <- updatePedigree(pop)
  }

  if (is.null(pop$epiNet) && pop$H2 > pop$h2) {
    message("Run attachEpiNet() to attach epistatic effects to population.")
  }

  return(pop)
}
