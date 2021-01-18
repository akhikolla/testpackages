#' Pairwise Geodesic Distances Between Two Sets of Data
#' 
#' Unlike \code{\link[RiemBase]{rbase.pdist}}, \code{rbase.pdist2} takes two sets of data \eqn{X=\{ x_i \}_{i=1}^m} and \eqn{Y=\{ y_j \}_{j=1}^m} 
#' and compute \eqn{mn} number of pairwise distances for all \eqn{i} and \eqn{j}.
#' 
#' @param input1 a S3 object of \code{riemdata} class, whose \code{$data} element is of length \eqn{m}. 
#' @param input2 a S3 object of \code{riemdata} class, whose \code{$data} element is of length \eqn{n}. 
#' @param parallel a flag for enabling parallel computation.
#' 
#' @return an \eqn{(m\times n)} matrix of pairwise distances.
#' 
#' @examples
#' ### Generate 10 2-frames in R^4 : as grassmann points
#' ndata = 10
#' data = array(0,c(4,2,ndata))
#' for (i in 1:ndata){
#'   tgt = matrix(rnorm(4*4),nrow=4)
#'   data[,,i] = qr.Q(qr(tgt))[,1:2]
#' }
#' 
#' gdata = riemfactory(data, name="grassmann")
#' 
#' ## Compute Pairwise Distances using pdist and pdist2
#' A = rbase.pdist(gdata)
#' B = rbase.pdist2(gdata,gdata)
#' 
#' ## Visual Comparison in Two Cases
#' opar = par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' image(A, col=gray((0:100)/100), main="pdist")
#' image(B, col=gray((0:100)/100), main="pdist2")
#' par(opar)
#' 
#' @export
rbase.pdist2 <- function(input1, input2, parallel=FALSE){
  #-------------------------------------------------------
  # must be of 'riemdata' class
  if ((class(input1))!="riemdata"){
    stop("* rbase.pdist2 : the input1 must be of 'riemdata' class. Use 'riemfactory' first to manage your data.")
  }
  if ((class(input2))!="riemdata"){
    stop("* rbase.pdist2 : the input2 must be of 'riemdata' class. Use 'riemfactory' first to manage your data.")
  }
  # acquire manifold name
  mfdname = tolower(input1$name)
  # other conditions : size and name should match
  if (tolower(input1$name)!=(input2$name)){
    stop("* pdist2 : two inputs should be of same manifold type.")
  }
  if (any(input1$size!=input2$size)){
    stop("* pdist2 : two inputs should be of same size and dimension.")
  }
  
  #-------------------------------------------------------
  # stack data as 3d matrices
  newdata1 = aux_stack3d(input1)
  newdata2 = aux_stack3d(input2)
  
  #-------------------------------------------------------
  # support of parallel computation using OpenMP and run
  # must be of 'riemdata' class
  nCores = parallel::detectCores()
  
  if (parallel==FALSE){
    output = engine_pdist2(newdata1, newdata2, mfdname)
  } else {
    if ((nCores==1)||(is.na(nCores))){
      output = engine_pdist2(newdata1, newdata2, mfdname)  
    } else {
      output = engine_pdist2_openmp(newdata1, newdata2, mfdname, nCores)
    }  
  }
  
  return(output)
}