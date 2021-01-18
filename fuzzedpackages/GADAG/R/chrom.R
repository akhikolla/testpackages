###########################################################
### returning the permutation matrix associated to a permutation ###############

##' @keywords internal
# @examples
#  ########################################################
#  # Loading toy data
#  ########################################################
#  c <- sample(10)
#
#  ########################################################
#  # Transforming in a matrix form
#  ########################################################
#  P <- chrom(c=c)
#
# @title Transform a permutation in a matrix form
# @description Internal function of the genetic algorithm to transform a permutation in a matrix form.
# @param c A permutation from [1,p].
# @return The pxp matrix of permutation associated to c.
# @author \packageAuthor{GADAG}
# @seealso \code{\link{GADAG}}, \code{\link{GADAG_Run}}.
#

chrom <- function(c){
  # INPUTS
  # c: permutation vector (size p)
  #
  # OUTPUTS
  # P: permutation matrix (p*p)

  p <- length(c)
  P <- matrix(0,p,p)
  I <- p*seq(from = 0, to = (p-1), by = 1)+c
  P[I] <- 1
  return(P)
}
