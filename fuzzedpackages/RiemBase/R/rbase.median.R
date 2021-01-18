#' Geometric Median of Manifold-valued Data
#' 
#' For manifold-valued data, geometric median is the solution of following cost function,
#' \deqn{\textrm{min}_x \sum_{i=1}^n \rho (x, x_i) = \sum_{i=1}^n \| \log_x (x_i) \|,\quad x\in\mathcal{M}}
#' for a given data \eqn{\{x_i\}_{i=1}^n}, \eqn{\rho(x,y)} the geodesic distance 
#' between two points on manifold \eqn{\mathcal{M}}, and \eqn{\| \log_x (y) \|} a logarithmic mapping onto the 
#' tangent space \eqn{T_x \mathcal{M}}. Weiszfeld's algorithms is employed.
#' 
#' @param input a S3 object of \code{riemdata} class. See \code{\link{riemfactory}} for more details.
#' @param maxiter maximum number of iterations for gradient descent algorithm.
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
#' ### Compute Geodesic Median
#' out1 = rbase.median(data)
#' out2 = rbase.median(data,parallel=TRUE) # test parallel implementation
#' }
#' 
#' 
#' @references 
#' \insertRef{fletcher_geometric_2009}{RiemBase}
#' 
#' \insertRef{aftab_generalized_2015}{RiemBase}
#' 
#' @author Kisung You
#' @export
rbase.median <- function(input, maxiter=496, eps=1e-6, parallel=FALSE){
  #-------------------------------------------------------
  # must be of 'riemdata' class
  if ((class(input))!="riemdata"){
    stop("* rbase.median : the input must be of 'riemdata' class. Use 'riemfactory' first to manage your data.")
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
  # calculate initial estimate
  tmpinit = engine_mean(newdata, mfdname, 10, as.double(eps))
  xinit   = tmpinit$x
  
  # must be of 'riemdata' class
  nCores = parallel::detectCores()
  if (parallel==FALSE){
    output = engine_median(newdata, mfdname, as.integer(maxiter), as.double(eps), xinit)
  } else {
    if ((nCores==1)||(is.na(nCores))){
      output = engine_median(newdata, mfdname, as.integer(maxiter), as.double(eps), xinit)
    } else {
      output = engine_median_openmp(newdata, mfdname, as.integer(maxiter), as.double(eps), nCores, xinit)
    }  
  }
  
  return(output)
}
