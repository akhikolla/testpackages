#' Prepare a S3 Class Object 'riemdata'
#' 
#' Most of the functions for \code{RiemBase} package require data to be wrapped as a \code{riemdata} class. 
#' Since manifolds of interests endow data points with specific constraints, the function \code{riemfactory}
#' first checks the requirements to characterize the manifold and then wraps the data into 
#' \code{riemdata} class, which is simply a list of manifold-valued data and the name of manifold. 
#' Manifold name input is, fortunately, \emph{case-insensitive}.
#' 
#' @param data data to be wrapped as \code{riemdata} class. Following input formats are considered,
#' \describe{
#' \item{2D array}{an \eqn{(m\times p)} matrix where data are stacked in columns over 2nd dimension. Appropriate for vector-valued \code{Euclidean} or \code{Sphere} manifold case.}
#' \item{3D array}{an \eqn{(m\times n\times p)} matrix where data are stacked in slices over 3rd dimension.}
#' \item{list}{unnamed list where each element of the list is a single data point. Sizes of all elements must match.}
#' }
#' @param name the name of Riemmanian manifold for data to which data belong.
#' 
#' @return a named \code{riemdata} S3 object containing
#' \describe{
#'   \item{data}{a list of manifold-valued data points.}
#'   \item{size}{size of each data matrix.}
#'   \item{name}{name of the manifold of interests.}
#' }
#' 
#' @examples 
#' \donttest{
#' # Test with Sphere S^2 in R^3 example
#' ## Prepare a matrix and list of 20 samples on S^2
#' sp.mat  = array(0,c(3,20)) # each vector will be recorded as a column
#' sp.list = list()
#' for (i in 1:20){
#'   tgt = rnorm(3)               # sample random numbers
#'   tgt = tgt/sqrt(sum(tgt*tgt)) # normalize
#'   
#'   sp.mat[,i]   = tgt   # record it as column vector
#'   sp.list[[i]] = tgt   # record it as an element in a list
#' }
#' 
#' ## wrap it using 'riemfactory'
#' rsp1 = riemfactory(sp.mat, name="Sphere")
#' rsp2 = riemfactory(sp.list, name="spHeRe")
#' }
#' 
#' @export
riemfactory <- function(data, name=c("euclidean","grassmann","spd","sphere","stiefel")){
  ############################################################
  ## MANIFOLD TYPE MATCHING
  allnames = tolower(c("Euclidean","Grassmann","SPD","Sphere","Stiefel"))
  name     = match.arg(tolower(name), allnames)
  ## I Believe Inputs should also be an array of length 2 or 3
  if (is.array(data)){
    datasize = dim(data)
    if (length(datasize)==2){
      newdata = list()
      for (i in 1:datasize[2]){
        newdata[[i]] = data[,i]
      }
      data = newdata
    } else if (length(datasize)==3){
      newdata = list()
      for (i in 1:datasize[3]){
        newdata[[i]] = data[,,i]
      }
      data = newdata
    } else {
      stop("* riemfactory : once 'data' is given as an array, it should be either 2D matrix or 3D array.")
    }
  }
  if (length(names(data))>1){
    names(data)=NULL
  }
  ## 0.2.0 : for matrix with one row or column, convert it into 
  for (i in 1:length(data)){
    tgt = data[[i]]
    if (is.matrix(tgt)){
      if ((nrow(tgt)==1)||(ncol(tgt)==1)){
        data[[i]] = as.vector(tgt)
      }
    }
  }
  switch(name,
         euclidean = stopifnot(islist_euclidean(data)),
         spd       = stopifnot(islist_spd(data)),
         sphere    = stopifnot(islist_sphere(data)),
         stiefel   = stopifnot(islist_stiefel(data)),
         grassman  = stopifnot(islist_stiefel(data))
  )
  
  ############################################################
  ## COMMON INFORMATION
  ## 1. transform vectors into matrices and name it as 'data'
  if (is.vector(data[[1]])){
    data = lapply(data, as.matrix)
  } else {
    data = data
  }
  ## 2. size
  size = c(nrow(data[[1]]), ncol(data[[1]]))
  
  output = list()
  output$data = data
  output$size = size
  output$name = name
  
  ############################################################
  ##  RETURN THE S3 CLASS
  return(structure(output, class="riemdata"))
}