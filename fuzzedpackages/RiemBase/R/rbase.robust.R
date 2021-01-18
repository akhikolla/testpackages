#' Robust Fréchet Mean of Manifold-valued Data
#' 
#' Robust estimator for mean starts from dividing the data \eqn{\{x_i\}_{i=1}^n} into \eqn{k} equally sized 
#' sets. For each subset, it first estimates Fréchet mean. It then follows a step to aggregate 
#' \eqn{k} sample means by finding a geometric median.
#' 
#' @param input a S3 object of \code{riemdata} class. See \code{\link{riemfactory}} for more details.
#' @param k number of subsets for which the data be divided into.
#' @param maxiter maximum number of iterations for gradient descent algorithm and Weiszfeld algorithm.
#' @param eps stopping criterion for the norm of gradient.
#' @param parallel a flag for enabling parallel computation.
#' 
#' @return a named list containing
#' \describe{
#' \item{x}{an estimate geometric median.}
#' \item{iteration}{number of iterations until convergence.}
#' }
#' 
#' @examples
#' \donttest{
#' ### Generate 100 data points on Sphere S^2 near (0,0,1).
#' ndata = 100
#' theta = seq(from=-0.99,to=0.99,length.out=ndata)*pi
#' tmpx  = cos(theta) + rnorm(ndata,sd=0.1)
#' tmpy  = sin(theta) + rnorm(ndata,sd=0.1)
#' 
#' ### Wrap it as 'riemdata' class
#' data  = list()
#' for (i in 1:ndata){
#'   tgt = c(tmpx[i],tmpy[i],1)
#'   data[[i]] = tgt/sqrt(sum(tgt^2)) # project onto Sphere
#' }
#' data = riemfactory(data, name="sphere")
#' 
#' ### Compute Robust Fréchet Mean
#' out1 = rbase.robust(data)
#' out2 = rbase.robust(data,parallel=TRUE) # test parallel implementation
#' }
#' 
#' @references 
#' \insertRef{2011arXiv1112.3914L}{RiemBase}
#' 
#' \insertRef{2013arXiv1308.1334M}{RiemBase}
#' 
#' \insertRef{2014arXiv1409.5937F}{RiemBase}
#' 
#' @seealso \code{\link[RiemBase]{rbase.mean}}, \code{\link[RiemBase]{rbase.median}}
#' @author Kisung You
#' @export
rbase.robust <- function(input, k=5, maxiter=496, eps=1e-6, parallel=FALSE){
  #-------------------------------------------------------
  # must be of 'riemdata' class
  if ((class(input))!="riemdata"){
    stop("* rbase.robust : the input must be of 'riemdata' class. Use 'riemfactory' first to manage your data.")
  }
  k = as.integer(k)
  if (k<=1){
    stop("* rbase.robust : when 'k' <= 1, there is no need to run this.")
  }
  # acquire manifold name
  mfdname = tolower(input$name)
  # stack data as 3d matrices
  newdata = aux_stack3d(input)
  if (is.matrix(newdata)){
    output = list()
    output$x = newdata
    output$iteration = 0
    return(output)
  }
  if (dim(newdata)[3]==1){
    output = list()
    output$x = matrix(newdata,nrow=nrow(newdata))
    output$iteration = 0
    return(output)
  }
  
  #-------------------------------------------------------
  # generate cluster index and separate true data
  nCores   = parallel::detectCores()
  clustidx = aux_rndivide(dim(newdata)[3], k)
  partdata = list()
  for (i in 1:length(clustidx)){
    tmpdata = newdata[,,clustidx[[i]]] 
    if (length(dim(tmpdata))==2){
      tmpdata = rmean_2to3(tmpdata)
    }
    if (length(dim(tmpdata))!=3){
      stop("* rmean : something is wrong.")
    }
    partdata[[i]] = tmpdata
  }
  
  #-------------------------------------------------------
  # let's run parallel for mean
  tmpout = list()
  for (i in 1:k){
    if ((nCores==1)||is.na(nCores)||(parallel=FALSE)){
      tmpoutput = rbase.mean.cube(partdata[[i]], mfdname, maxiter=maxiter, eps=eps)
    } else {
      tmpoutput = rbase.mean.cube(partdata[[i]], mfdname, maxiter=maxiter, eps=eps, parallel=TRUE)
    }
    tmpout[[i]] = tmpoutput$x
  }
  
  #-------------------------------------------------------
  # it's time to run robust median problem
  partmeans = rmean_lto3(tmpout)

  if ((nCores==1)||(is.na(nCores))||(parallel=FALSE)){
    partran  = engine_mean(partmeans, mfdname, as.integer(maxiter), as.double(eps))
    partinit = partran$x
    output   = engine_median(partmeans, mfdname, as.integer(maxiter), as.double(eps), partinit)
  } else {
    partran  = engine_mean_openmp(partmeans, mfdname, as.integer(maxiter), as.double(eps), nCores)
    partinit = partran$x
    output   = engine_median_openmp(partmeans, mfdname, as.integer(maxiter), as.double(eps), nCores, partinit)
  }

  #-------------------------------------------------------
  # return output
  return(output)
}

#' @keywords internal
#' @noRd
rmean_2to3 <- function(A){
  m = nrow(A)
  p = ncol(A)
  
  output = array(0,c(m,1,p))
  for (i in 1:p){
    output[,,i] = A[,i]
  }
  return(output)
}
#' @keywords internal
#' @noRd
rmean_lto3 <- function(dlist){
  p = length(dlist)
  m = nrow(dlist[[1]])
  n = ncol(dlist[[1]])
  
  output = array(0,c(m,n,p))
  for (i in 1:p){
    output[,,i] = dlist[[i]]
  }
  return(output)
}