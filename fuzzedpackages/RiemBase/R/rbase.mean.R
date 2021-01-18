#' Fréchet Mean of Manifold-valued Data
#' 
#' For manifold-valued data, Fréchet mean is the solution of following cost function,
#' \deqn{\textrm{min}_x \sum_{i=1}^n \rho^2 (x, x_i),\quad x\in\mathcal{M}}
#' for a given data \eqn{\{x_i\}_{i=1}^n} and \eqn{\rho(x,y)} is the geodesic distance 
#' between two points on manifold \eqn{\mathcal{M}}. It uses a gradient descent method 
#' with a backtracking search rule for updating.
#' 
#' @param input a S3 object of \code{riemdata} class. See \code{\link{riemfactory}} for more details.
#' @param maxiter maximum number of iterations for gradient descent algorithm.
#' @param eps stopping criterion for the norm of gradient.
#' @param parallel a flag for enabling parallel computation.
#' 
#' @return a named list containing
#' \describe{
#' \item{x}{an estimate Fréchet mean.}
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
#' ### Compute Fréchet Mean
#' out1 = rbase.mean(data)
#' out2 = rbase.mean(data,parallel=TRUE) # test parallel implementation
#' }
#' 
#' @references 
#' \insertRef{karcher_riemannian_1977}{RiemBase}
#' 
#' \insertRef{kendall_probability_1990}{RiemBase}
#' 
#' \insertRef{afsari_convergence_2013}{RiemBase}
#' 
#' @author Kisung You
#' @export
rbase.mean <- function(input, maxiter=496, eps=1e-6, parallel=FALSE){
  #-------------------------------------------------------
  # must be of 'riemdata' class
  if ((class(input))!="riemdata"){
    stop("* rbase.mean : the input must be of 'riemdata' class. Use 'riemfactory' first to manage your data.")
  }
  # acquire manifold name
  mfdname = tolower(input$name)
  # stack data as 3d matrices
  newdata = aux_stack3d(input)
  if (length(input$data)==1){
    output = list()
    output$x = input$data[[1]]
    output$iteration = 0
    return(output)
  }
  
  #-------------------------------------------------------
  # calculate
  nCores = parallel::detectCores()
  if (parallel==FALSE){
    output = engine_mean(newdata, mfdname, as.integer(maxiter), as.double(eps))
  } else {
    if ((nCores==1)||(is.na(nCores))){
      output = engine_mean(newdata, mfdname, as.integer(maxiter), as.double(eps))
    } else {
      output = engine_mean_openmp(newdata, mfdname, as.integer(maxiter), as.double(eps), nCores)
    }  
  }
  
  return(output)
}

#' @keywords internal
#' @noRd
rbase.mean.cube <- function(datacube, mfdname, maxiter=496, eps=1e-6, parallel=FALSE){
  #-------------------------------------------------------
  newdata = datacube
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
  # calculate
  nCores = parallel::detectCores()
  
  # must be of 'riemdata' class
  nCores = parallel::detectCores()
  if (parallel==FALSE){
    output = engine_mean(newdata, mfdname, as.integer(maxiter), as.double(eps))
  } else {
    if ((nCores==1)||(is.na(nCores))){
      output = engine_mean(newdata, mfdname, as.integer(maxiter), as.double(eps))
    } else {
      output = engine_mean_openmp(newdata, mfdname, as.integer(maxiter), as.double(eps), nCores)
    }  
  }
  
  return(output)
}
# mydata = list()
# sdval  = 0.1
# diag8  = diag(8)
# for (i in 1:20){
#   mydata[[i]] = qr.Q(qr(diag8[,1:4] + matrix(rnorm(8*4,sd=sdval),ncol=4)))
# }
# grdata = riemfactory(mydata, name = "stiefel")
# 
# par(mfrow=c(1,3))
# image(rbase.mean(grdata)$x)
# image(rbase.median(grdata)$x)
# image(rbase.robust(grdata)$x)