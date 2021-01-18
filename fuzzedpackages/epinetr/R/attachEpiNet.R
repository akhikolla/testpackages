#' Attach epistatic network to population.
#'
#' Constructs and attaches a new epistatic network between QTLs to a given
#' \code{Population} object.
#'
#' \code{attachEpiNet()} can be used to construct a new epistatic network based
#' either on stochastic processes or by adapting a user-supplied incidence
#' matrix.
#'
#' If \code{scaleFree} is \code{FALSE}, a network is constructed at random; if
#' it is \code{TRUE}, however, a scale-free network is constructed between QTLs
#' using the Barabasi-Albert model. (The random network is constructed using
#' the same algorithm but without preferential attachment.)
#'
#' A minimal initial network of randomly selected QTLs is first constructed,
#' before then growing the network using randomly selected QTLs, preferentially
#' attaching each QTL to at least \code{m} other QTLs (in the scale-free case).
#' For increasing orders of interaction, degrees from lower orders of
#' interaction are used when determining preferential attachment.
#'
#' If a user-supplied incidence matrix is given via the \code{incmat}
#' argument, the \code{scaleFree}, \code{k} and \code{m} arguments are ignored,
#' and the network structure is instead derived from the incidence matrix. The
#' matrix should be such that the rows represent QTLs and the columns represent
#' interactions between QTLs: a \code{1} at both \code{incmat[i1,j]} and
#' \code{incmat[i2,j]} means there is an interaction between the QTLs
#' \code{i1} and \code{i2}.
#'
#' The contributions that interactions make towards the phenotypic value are
#' generated randomly when this function is called, based on a normal
#' distribution. By default, each order of interaction has the same variance;
#' however, using the \code{varfun} argument, a function can be supplied to
#' alter the variance per order of interaction.
#'
#' Internally, each QTL is assigned a randomly generated value per
#' interaction per genotype (i.e. the heterozygous genotype and the two
#' homozygous genotypes), and these values are summed according to the value
#' of the genotype.
#'
#' @param pop the \code{Population} object to which the epistatic network will
#'   be attached
#' @param scaleFree an optional logical value indicating whether to construct
#'   a network using the Barabasi-Albert model for generating scale-free
#'   networks
#' @param k an optional vector listing orders of interaction to include
#' @param m an optional integer indicating the minimum number of interactions
#'   for each QTL; see details below
#' @param additive an optional integer indicating the number of QTLs which
#'   will remain purely additive
#' @param incmat an optional incidence matrix specifying a user-defined
#'   network
#'
#' @return A copy of the supplied \code{Population} is returned, with the new
#'   epistatic network attached.
#' @author Dion Detterer, Paul Kwan, Cedric Gondro
#' @references Barabasi AL, Albert R, 'Emergence of scaling in random networks,'
#'   \emph{Science} 286(5439): 509-12, 15 October 1999.
#' @export
#'
#' @examples
#' \donttest{
#' # Generate a population and attach additive effects
#' pop <- Population(
#'   popSize = 200, map = map100snp, QTL = 20,
#'   alleleFrequencies = runif(100),
#'   broadH2 = 0.9, narrowh2 = 0.6, traitVar = 40
#' )
#' pop <- addEffects(pop)
#'
#' # Attach a random epistatic network with two- to four-way
#' # interactions between QTLs
#' popRnd <- attachEpiNet(pop, k = 2:4)
#'
#' # Plot random network
#' plot(getEpiNet(popRnd))
#'
#' # Attach a scale-free epistatic network with two-way interactions
#' # between QTLs and a minimum of three interactions per QTL
#' popSF <- attachEpiNet(pop, scaleFree = TRUE, m = 3)
#'
#' # Plot scale-free network
#' plot(getEpiNet(popSF))
#'
#' # Attach user-defined epistatic network
#' popUser <- attachEpiNet(pop, incmat = rincmat100snp)
#'
#' # Plot user-defined network
#' plot(getEpiNet(popUser))
#'
#' # Attach a random epistatic network with two- to ten-way
#' # interactions between QTLs and decaying variance
#' popDecay <- attachEpiNet(pop, k = 2:10)
#' }
#' @seealso \code{\link{Population}}, \code{\link{getEpiNet}},
#' \code{\link{plot.EpiNet}}, \code{\link{addEffects}}
attachEpiNet <- function(pop, scaleFree = FALSE, k = 2, m = 1, additive = 0,
                         incmat = NULL) {
  testPop(pop)

  if (pop$h2 == pop$H2) {
    stop("No epistatic interaction in population")
  }

  n <- length(pop$qtl)

  if (is.matrix(incmat)) {
    if (nrow(incmat) != n) {
      stop("Number of rows in incmat must match number of QTL in pop")
    }

    # If incmat is a matrix, use that as a user-defined network
    pop$epiNet <- userNetwork(incmat)
  } else {
    pop$epiNet <- buildNetwork(n, k, m, additive, scaleFree, pop)
  }

  # Set up powers of 3 in the incidence matrix
  pow <- function(x) {
    indices <- which(x != 0)
    for (i in 1:length(indices)) x[indices[i]] <- 3^(i - 1)
    return(x)
  }
  pop$epiNet$Incidence <- apply(pop$epiNet$Incidence, 2, pow)

  # Calculate offset and scaling factor for epistasis
  pop <- calcEpiScale(pop)

  # Calculate phenotypic components for pedigree data frame if additive
  # effects are present or unneeded (but not both)
  if ((!is.null(pop$additive) || pop$h2 == 0) && (pop$h2 > 0 || is.null(pop$additive))) {
    pop <- updatePedigree(pop)
  }

  if (is.null(pop$additive) && pop$h2 > 0) {
    message("Run addEffects() to attach additive effects to population.")
  }

  return(pop)
}
