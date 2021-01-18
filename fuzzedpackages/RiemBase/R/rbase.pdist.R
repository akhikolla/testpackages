#' Pairwise Geodesic Distances of a Data Set
#' 
#' Geodesic distance \eqn{\rho(x,y)} is the length of (locally) shortest path 
#' connecting two points \eqn{x,y\in\mathcal{M}}. Some manifolds have closed-form expression, while 
#' others need numerical approximation. 
#' 
#' @param input a S3 object of \code{riemdata} class, whose \code{$data} element is of length \eqn{n}. See \code{\link{riemfactory}} for more details.
#' @param parallel a flag for enabling parallel computation.
#' 
#' @return an \eqn{(n\times n)} matrix of pairwise distances.
#' 
#' @examples
#' ### Generate 10 2-frames in R^4
#' ndata = 10
#' data = array(0,c(4,2,ndata))
#' for (i in 1:ndata){
#'   tgt = matrix(rnorm(4*4),nrow=4)
#'   data[,,i] = qr.Q(qr(tgt))[,1:2]
#' }
#' 
#' ## Compute Pairwise Distances as if for Grassmann and Stiefel Manifold
#' A = rbase.pdist(riemfactory(data,name="grassmann"))
#' B = rbase.pdist(riemfactory(data,name="stiefel"))
#' 
#' ## Visual Comparison in Two Cases
#' opar = par(no.readonly=TRUE)
#' par(mfrow=c(1,2))
#' image(A, col=gray((0:100)/100), main="Grassmann")
#' image(B, col=gray((0:100)/100), main="Stiefel")
#' par(opar)
#' 
#' @export
rbase.pdist <- function(input, parallel=FALSE){
  #-------------------------------------------------------
  # must be of 'riemdata' class
  if ((class(input))!="riemdata"){
    stop("* rbase.pdist : the input must be of 'riemdata' class. Use 'riemfactory' first to manage your data.")
  }
  # acquire manifold name
  mfdname = tolower(input$name)
  # stack data as 3d matrices
  newdata = aux_stack3d(input)
  
  #-------------------------------------------------------
  # support of parallel computation using OpenMP and run
  # must be of 'riemdata' class
  nCores = parallel::detectCores()
  
  if ((nCores==1)||(is.na(nCores))||(parallel==FALSE)){
    output = engine_pdist(newdata, mfdname)
  } else {
    output = engine_pdist_openmp(newdata, mfdname, nCores)
  }  
  

  #-------------------------------------------------------
  # return the matrix
  return(output)
}