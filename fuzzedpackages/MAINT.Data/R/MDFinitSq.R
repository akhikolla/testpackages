# Implementation of Hardin and Rocke (2005) F-distribution approximation of the robust squared raw MCD Mahalanobis distances

# Main function: qHardRoqF wich returns the p-quantile of the relevant F distribution, for a raw MCD based on nobs observations,
#                nvar variables, and h observations kept in the trimmed estimator

# by Pedro Duarte Silva (April 2017)    

# Reference:

# Hardin, J. and Rocke, D.M. (2005) "The distribution of robust distances" Journal of Computational and Graphical Statistics 14, 910-927. 

casy <- function(nobs,nvar,h=floor((nobs+nvar+1)/2))
{
  hovern <- h/nobs
  q <- qchisq(hovern,nvar)

  pchisq(q,nvar+2)/hovern
}

CHmasy <- function(nobs,nvar,h=floor((nobs+nvar+1)/2))
{
  alpha <- (nobs-h)/nobs
  oneminusalpha <- 1.-alpha
  q_alpha <- qchisq(oneminusalpha,nvar)
  c_0 <- pchisq(q_alpha,nvar+2) 
  c_alpha <- oneminusalpha/c_0
  c_2 <- -c_0/2
  c_3 <- -pchisq(q_alpha,nvar+4)/2
  c_4 <- 3*c_3
  b_1 <- c_alpha*(c_3-c_4)/oneminusalpha
  b_2 <- .5 + (c_alpha/oneminusalpha) * (c_3-(q_alpha/nvar)*(c_2+oneminusalpha/2))
  v_1 <- oneminusalpha*b_1^2*(alpha*(c_alpha*q_alpha/nvar-1)^2-1) - 
    2*c_3*c_alpha^2 * (3*(b_1-nvar*b_2)^2+(nvar+2)*b_2*(2*b_1-nvar*b_2))
  v_2 <- nobs * (b_1*(b_1-nvar*b_2)*oneminusalpha)^2 * c_alpha^2
  v <- v_1/v_2
  
  2./(v*c_alpha^2)
}

HDmpred <- function(nobs,nvar,h=floor((nobs+nvar+1)/2),masy=NULL)
{
  if (is.null(masy)) masy <- CHmasy(nobs,nvar,h=h)
  masy * exp(.725-0.00663*nvar-0.078*log(nobs))
}

qHardRoqF <- function(p,nobs,nvar,h=floor((nobs+nvar+1)/2),adj=TRUE,lower.tail=TRUE,log.p=FALSE)
{
  c <- casy(nobs,nvar,h=h)
  if (adj) m <- HDmpred(nobs,nvar,h=h)
  else m <- CHmasy(nobs,nvar,h=h)
  df2 <- m-nvar+1

  (nvar*m/(c*df2)) * qf(p,nvar,df2,lower.tail=lower.tail,log.p=log.p)
}  
