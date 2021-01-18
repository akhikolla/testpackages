#' @title Parameters for variable block structure search.
#' @name vbSearchControl
#' @description This function creates a list with parameters for the search of a variable 
#' block structure used as an argument for \code{hmmvbTrain}.
#' @param nperm The number of variable permutations. This parameter is ignored 
#' if permutations are provided in \code{perm} argument.
#' @param dim Data dimensionality. Must be provided with \code{perm} argument,
#' otherwise is ignored. 
#' @param perm A list of integer vectors specifying variable permutations. If
#' provided, the argument \code{dim} must be supplied.
#' @param minDim Minimum variable block dimension. Should be an integer equal to 1 or 2.
#' @param maxDim Maximum variable block dimension. 
#' @param numstPerDim An integer vector of length \code{maxDim} specifying a map from
#' the variable block dimensionality to the number of states in the block. 
#' \emph{k}th value in the vector corresponds to number of states for dimensionality \emph{k}.
#' @param relax A logical value indicating whether or not variable block structure search 
#' will be performed under less restricting conditions.
#' @return The named list with parameters.
#' @seealso \code{\link{hmmvbTrain}}
#' @examples 
#' # setting up permutations
#' perm <- list(c(1,2,3), c(1,3,2), c(3,2,1))
#' searchControl <- vbSearchControl(perm=perm, dim=3)
#' 
#' # setting up a map between block dimensionality and number of states
#' searchControl <- vbSearchControl(maxDim=5, numstPerDim=c(3,4,5,6,7))
vbSearchControl <- function(perm=NULL, numstPerDim=NULL, dim=NULL,
                              maxDim=10, minDim=1, nperm=1, relax=FALSE){

  if (!is.null(perm)){
    if (!is.list(perm))
      stop('Permutation argument should be a list!\n')
    if (is.null(dim))
      stop('dim argument should be set in case permutations are provided!\n')
    dim = as.integer(dim)
    if ((length(dim)!=1) || (dim <= 0))
      stop('Invalid dim argument provided! dim should be a scalar positive integer\n')

    for (i in 1:length(perm)){
      if (!isTRUE(all.equal(perm[[i]], as.integer(perm[[i]]))))
        stop('Invalid permutation ',i,' is provided! Permutations should be integer arrays\n')
      if (length(perm[[i]]) != dim)
        stop('Length of permutation ',i,' does not match dim!\n')
      if (length(unique(perm[[i]])) != length(perm[[i]]))
        stop('Permutation ',i,' contains repeats!\n')
    }

    message('Since permutations were provided, argument nperm will be ignored!\n')
  }

  if (!is.null(dim) && is.null(perm))
    message('dim argument is ignored if no permutations are provided!\n')

  if (!is.null(numstPerDim)){
    numstPerDim = as.integer(numstPerDim)
    if (length(numstPerDim) != maxDim || any(numstPerDim < 0) )
      stop('Invalid numstatesPerDim argument provided! numstatesPerDim should
           be a positive integer array with length maxDim')
  }

  maxDim = as.integer(maxDim)
  if ((length(maxDim)>1) || (maxDim <= 0))
    stop('Invalid maxDim argument provided! maxDim should be a scalar positive integer\n')

  
  if ( (minDim != 1L) && (minDim != 2L) )
    stop('Only values 1 and 2 are supported for minDim argument!\n')

  nperm = as.integer(nperm)
  if ((length(nperm)!=1) || (nperm <= 0))
    stop('Invalid nperm argument provided! nperm should be a scalar positive integer\n')

  if (!is.logical(relax) || length(relax)!=1)
    stop('Invalid relax argument provided! relax should be a scalar logical\n')


  return(list(perm=perm, numstPerDim=numstPerDim, maxDim=maxDim,
              minDim=minDim, nperm=nperm, relax=relax))
}
