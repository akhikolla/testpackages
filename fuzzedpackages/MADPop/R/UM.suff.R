#' Sufficient statistics for the Unconstrained-Multinomial model
#'
#' Converts a matrix or data.frame of genotype data into the sufficient statistics required to fit a Dirichlet-Multinomial hierarchical model.
#'
#' @param X Genotype adata.  Either a \code{N x 4} matrix with \code{NA}'s indicating duplicates or a \code{N x 5} column data.frame with the first column being the \code{popId}.
#' @param popId grouping variable for \code{X}. Must be supplied if \code{X} has 4 columns.
#' @return A list with elements:
#' \itemize{
#' \item \code{A}: Vector of unique alleles names.  The allele numbers in the following quantities correspond to the indices of \code{A}.
#' \item \code{G}: 4-column matrix of unique genotype combinations.  The presence of 0's indicates that less than four alleles were amplified indicating that a given genotype either has less than 4 distinct alleles or that some alleles are duplicated.
#' \item \code{tab}: Observed data in a simplified numerical format.  This is a contingency table with rows given by the unique elements of \code{popId} and columns given by each row of \code{G}.
#' }
#' @examples
#' # sufficient statistics in 3 lakes
#'
#' X <- fish215[fish215$Lake %in% c("Hogan", "Manitou", "Simcoe"),]
#' suff <- UM.suff(X)
#'
#' suff$A # allele names
#' suff$G # unique genotypes
#' suff$tab # contingency table
#'@export
UM.suff <- function(X, popId) {
  # get popId
  if(ncol(X) == 5) {
    popId <- X[,1]
    X <- X[,-1]
  } else if(ncol(X) == 4) {
    if(missing(popId)) {
      stop("Need to specify popId, either explicitly or implicitly as first column of X.")
    }
  } else {
    stop("X must be a matrix or data.frame with 4 or 5 columns.")
  }
  # format gene data
  af <- geno.format(X, Y.only = FALSE)
  # create contingency table
  # gen <- apply(af$Y, 1,
  #              function(x) paste0(sort(x), collapse = "."))
  gen <- apply(af$Y, 1,
               function(x) paste0(sort(x[x>0]), collapse = "."))
  tab <- cbind(popId = as.character(popId), genotype = gen)
  tab <- table(tab[,1], tab[,2])
  list(A = af$A, G = af$G, tab = tab)
}
