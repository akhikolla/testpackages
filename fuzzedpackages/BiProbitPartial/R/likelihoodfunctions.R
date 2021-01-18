#' log likelihood of bivariate probit with partial observability
#'
#' @param theta numeric vector of dimension equal to that of the free parameter space
#' @param X1 numeric matrix of covariates for the first equation
#' @param X2 numeric matrix of covariates for the second equation
#' @param Z numeric matrix or column vecotr of response observations
#' @param rho numeric value for rho if fixed
#' @param p numeric precomputed probabilities of Pr(Y1=1,Y2=1)
#' @param summed logical if the log likelihood observations should be summed
#' @param fixrho logical if rho should be fixed
#'
#' @return if summed is TRUE then the function returns the numeric sum of the likelihood vector, else it returns a numeric vector with each entry a value of the likelihood vector
#'
#' @export
#'
llhood1 <- function(theta,X1,X2,Z,rho=0, p = NULL, summed = T, fixrho = F){
  k1 = dim(X1)[2]
  k2 = dim(X2)[2]
  k = k1 + k2
  
  beta1 = theta[1:k1]
  beta2 = theta[(k1+1):k]
  if(!fixrho)
    rho = utils::tail(theta,1)
  
  if(is.null(p))
    p = pbivnorm::pbivnorm(as.vector(X1 %*% beta1), as.vector(X2 %*% beta2), rho)
  
  lpZ1 = log(p[Z==1])
  lpZ0 = log(1-p[Z==0])
  
  if(summed){
    s1 = if(any(is.nan(lpZ1)))-.Machine$double.xmax else sum(lpZ1)
    s2 = if(any(is.nan(lpZ0)))-.Machine$double.xmax else sum(lpZ0)
    out = s1 + s2
  }else{
    out = numeric(length(Z))
    out[Z==1] = lpZ1
    out[Z==0] = lpZ0
  }
    
  return(out)
  
}


#' Gradient of bivariate probit with partial observability
#'
#' @param theta numeric vector of dimension equal to that of the free parameter space
#' @param X1 numeric matrix of covariates for the first equation
#' @param X2 numeric matrix of covariates for the second equation
#' @param Z numeric matrix or column vecotr of response observations
#' @param rho numeric value for rho if fixed
#' @param p numeric precomputed probabilities of Pr(Y1=1,Y2=1)
#' @param summed logical if the gradient observations should be summed
#' @param fixrho logical if rho should be fixed
#'
#' @return if summed is TRUE then the function returns the numeric column sum of the gradient matrix, else it returns a numeric vector with each entry a value of the gradient vector
#' @export
#'
grad1 <- function(theta,X1,X2,Z,rho = 0, p = NULL, summed = T, fixrho = F){
  Z1Dex = Z==1
  Z0Dex = Z==0
  N = length(Z)
  NZ1 = sum(Z1Dex)
  NZ0 = N - NZ1
  k1 = dim(X1)[2]
  k2 = dim(X2)[2]
  k = k1 + k2
  
  beta1 = theta[1:k1]
  beta2 = theta[(k1+1):k]
  if(!fixrho)
    rho = utils::tail(theta,1)
  
  if(is.null(p))
    p = pbivnorm::pbivnorm(as.vector(X1 %*% beta1), as.vector(X2 %*% beta2), rho)
  
  eq1a = X1 * matrix(stats::dnorm(X1 %*% beta1) * 
                       stats::pnorm((X2 %*% beta2 - rho * X1 %*% beta1)/sqrt(1-rho^2)),
                    N,k1,byrow =F)
  eq2a = X2 * matrix(stats::dnorm(X2 %*% beta2) * 
                       stats::pnorm((X1 %*% beta1 - rho * X2 %*% beta2)/sqrt(1-rho^2)),
                    N,k2,byrow =F)

  denom = p
  denom[Z0Dex] = -(1-p[Z0Dex])
  
  eq1b = eq1a/matrix(denom,N,k1,byrow=F)
  eq2b = eq2a/matrix(denom,N,k2,byrow=F)

  if(!fixrho){
    eq3a = mvtnorm::dmvnorm(cbind(X1 %*% beta1,X2 %*% beta2), c(0,0), matrix(c(1,rho,rho,1),2,2))
    eq3b = eq3a
    eq3b= eq3a/matrix(denom,N,1,byrow=F)
    
    if(!summed)
      out = cbind(eq1b,eq2b,eq3b) 
    else
      out = colSums(cbind(eq1b,eq2b,eq3b))
  }else{
    if(!summed)
      out = cbind(eq1b,eq2b) 
    else
      out = colSums(cbind(eq1b,eq2b))
  }
  
  return(out)
  
  
}