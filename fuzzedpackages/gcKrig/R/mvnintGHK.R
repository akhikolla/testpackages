


mvnintGHKOcpp_R <- function(mu, R, lower, upper, nrep){
  .Call('mvnintGHKOcpp', PACKAGE = 'gcKrig', mu, R, lower, upper, nrep)
}


mvnintGHKcpp_R <- function(mu, R, lower, upper, nrep){
  .Call('mvnintGHKcpp', PACKAGE = 'gcKrig', mu, R, lower, upper, nrep)
}


mvnintGHK <- function(mean, sigma, lower, upper, nrep = 5000, log = TRUE, reorder = TRUE){
  if(!is.matrix(sigma))
    stop("Input 'sigma' must be of form matrix!")

  if(!isSymmetric(sigma))
    stop("Input covariance matrix 'sigma' must be symmetric!")

  if(length(lower)== 1 & lower[1] == -Inf) lower <- rep(-.Machine$double.xmax, nrow(sigma))
  if(length(upper)== 1 & upper[1] == Inf) upper <- rep(.Machine$double.xmax, nrow(sigma))

  if( inherits(try(chol(sigma),silent=TRUE),"try-error") )
    stop("Cholesky Decomposition failed. Input matrix sigma is not a valid covariance matrix!")

  if(!all.equal(length(mean), nrow(sigma), length(lower), length(upper)))
    stop("Input 'mean', lower' and 'upper' must have same length as dimension of the sigma!")

  lower <- ifelse(lower == -Inf, -.Machine$double.xmax, lower)
  upper <- ifelse(upper == Inf, .Machine$double.xmax, upper)

  if(!all(lower <= upper))
  stop("Elements in 'lower' must be <=  the corresponding elements in 'upper'!")

  if(reorder == TRUE){
    ans <- mvnintGHKOcpp_R(mu = mean, R = sigma, lower = lower, upper = upper, nrep = nrep)
  }else{
    ans <- mvnintGHKcpp_R(mu = mean, R = sigma, lower = lower, upper = upper, nrep = nrep)
  }
  if(ans$value < -.Machine$double.max.exp)
    stop("Computation Failed Due to Numerical Problem or Large Dimensionality!")
  if(log == F) ans$value <- exp(ans$value)
  return(ans)
}






# some comparisons
# unix.time(
# mvnintGHK(mean = rep(0, 121), sigma =  diag(0.2, 121) + matrix(0.8, 121, 121),
# lower = rep(-2,121), upper = rep(2,121), nrep = 10000)
# )
#
# log(pmvnorm(mean = rep(0, 121), sigma = diag(0.2, 121) + matrix(0.8, 121, 121),
#             lower = rep(-2,121), upper = rep(2,121)))
#
# intgrandtmp <- function(x){
#   lower = rep(-2,121)
#   upper = rep(2,121)
#
#   tmp = pnorm((upper - sqrt(0.8)*qnorm(x))/sqrt(1-0.8))-
#     pnorm((lower - sqrt(0.8)*qnorm(x))/sqrt(1-0.8))
#
#   return(prod(tmp))
# }
#
#   tmp2 = Vectorize(intgrandtmp)
#   tmp = log(integrate(tmp2, lower = .Machine$double.eps, upper=1-.Machine$double.eps, subdivisions = 10000,
#                   rel.tol = .Machine$double.eps^0.25, abs.tol = .Machine$double.xmin,
#                   stop.on.error = FALSE)$value)
# tmp







# some comparisons
# unix.time(
# mvnintGHK(mean = rep(0, 50), sigma =  diag(0.2, 50) + matrix(0.8, 50, 50),
# lower = rep(-3,50), upper = rep(2,50), nrep = 10000)
# )
# unix.time(
# log(pmvnorm(mean = rep(0, 50), sigma = diag(0.2, 50) + matrix(0.8, 50, 50),
#             lower = rep(-3,50), upper = rep(2,50)))
# )
# intgrandtmp <- function(x){
#   lower = rep(-3,50)
#   upper = rep(2,50)
#
#   tmp = pnorm((upper - sqrt(0.8)*qnorm(x))/sqrt(1-0.8))-
#     pnorm((lower - sqrt(0.8)*qnorm(x))/sqrt(1-0.8))
#
#   return(prod(tmp))
# }
#
#   tmp2 = Vectorize(intgrandtmp)
#   tmp = log(integrate(tmp2, lower = .Machine$double.eps, upper=1-.Machine$double.eps, subdivisions = 10000,
#                 rel.tol = .Machine$double.eps^0.25, abs.tol = .Machine$double.xmin,
#                   stop.on.error = FALSE)$value)
# tmp
