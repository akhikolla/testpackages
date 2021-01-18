#' Uncorrected testing procedure for correlations.
#'
#' @param data         matrix of observations
#' @param alpha        level of multiple testing (used if logical=TRUE)
#' @param stat_test  
#' \describe{
#'   \item{'empirical'}{\eqn{\sqrt{n}*abs(corr)}}
#'   \item{'fisher'}{   \eqn{\sqrt{n-3}*1/2*\log( (1+corr)/(1-corr) )}}
#'   \item{'student'}{  \eqn{\sqrt{n-2}*abs(corr)/\sqrt(1-corr^2)}}
#'   \item{'2nd.order'}{ \eqn{\sqrt{n}*mean(Y)/sd(Y)} with \eqn{Y=(X_i-mean(X_i))(X_j-mean(X_j))}}
#' }
#' @param vect         if TRUE returns a vector of adjusted p-values, corresponding to \code{vectorize(cor(data))};
#'                     if FALSE, returns an array containing the adjusted p-values for each entry of the correlation matrix 
#' @param logical      if TRUE, returns either a vector or a matrix where each element is equal to TRUE if the corresponding null hypothesis is rejected, and to FALSE if it is not rejected
#' @param arr.ind      if TRUE, returns the indexes of the significant correlations, with respect to level alpha
#'
#' @return Returns  \itemize{\item{the non-adjusted p-values, as a vector or a matrix depending of the value of \code{vect},} \item{an array containing indexes \eqn{\lbrace(i,j),\,i<j\rbrace} for which correlation between variables \eqn{i} and \eqn{j} is significant, if \code{arr.ind=TRUE}.}}
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats pnorm
#' @export
#'
#' @examples 
#' n <- 100
#' p <- 10
#' corr_theo <- diag(1,p)
#' corr_theo[1,3] <- 0.5
#' corr_theo[3,1] <- 0.5
#' data <- MASS::mvrnorm(n,rep(0,p),corr_theo)
#' # p-values
#' res <- UncorrectedCor(data,stat_test='empirical')
#' round(res,2)
#' # significant correlations with level alpha:
#' alpha <- 0.05
#' whichCor(res<alpha)
#' # directly
#' UncorrectedCor(data,alpha,stat_test='empirical',arr.ind=TRUE)
UncorrectedCor <- function(data,alpha=0.05,stat_test='empirical',vect=FALSE,logical=FALSE,arr.ind=FALSE){

  if(sum(stat_test==c('empirical','fisher','student','2nd.order'))==0){ stop('Wrong value for stat_test.')}
  
	stat <- abs(eval_stat(data,stat_test))
  m <- length(stat)
  pval <- 2*(1-pnorm(stat))

    if(arr.ind){ 
      vect <- FALSE 
      logical <- TRUE
    }
    
  res <- pval  
  if(!vect){ res <- unvectorize(pval) }
  if(logical){ res <- (res < alpha) }
  if(arr.ind){ res <- whichCor(res) }
  
  return(res)
  
  
}



#--------------------------------method 1 : Bonferroni--------------------------------------
#--------------------------------------------------------------------------------------------

#' Bonferroni multiple testing procedure for correlations.
#'
#' @param data         matrix of observations
#' @param alpha        level of multiple testing (used if logical=TRUE)
#' @param stat_test  
#' \describe{
#'   \item{'empirical'}{\eqn{\sqrt{n}*abs(corr)}}
#'   \item{'fisher'}{   \eqn{\sqrt{n-3}*1/2*\log( (1+corr)/(1-corr) )}}
#'   \item{'student'}{  \eqn{\sqrt{n-2}*abs(corr)/\sqrt(1-corr^2)}}
#'   \item{'2nd.order'}{ \eqn{\sqrt{n}*mean(Y)/sd(Y)} with \eqn{Y=(X_i-mean(X_i))(X_j-mean(X_j))}}
#' }
#' @param vect         if TRUE returns a vector of adjusted p-values, corresponding to \code{vectorize(cor(data))};
#'                     if FALSE, returns an array containing the adjusted p-values for each entry of the correlation matrix 
#' @param logical      if TRUE, returns either a vector or a matrix where each element is equal to TRUE if the corresponding null hypothesis is rejected, and to FALSE if it is not rejected
#' @param arr.ind      if TRUE, returns the indexes of the significant correlations, with respect to level alpha
#'
#' @return Returns  \itemize{\item{the adjusted p-values, as a vector or a matrix depending of the value of \code{vect},} \item{an array containing indexes \eqn{\lbrace(i,j),\,i<j\rbrace} for which correlation between variables \eqn{i} and \eqn{j} is significant, if \code{arr.ind=TRUE}.}}
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats pnorm
#' @export
#'
#' @references  Bonferroni, C. E. (1935). Il calcolo delle assicurazioni su gruppi di teste. Studi in onore del professore salvatore ortu carboni, 13-60.
#' @references  Roux, M. (2018). Graph inference by multiple testing with application to Neuroimaging, Ph.D., Université Grenoble Alpes, France, https://tel.archives-ouvertes.fr/tel-01971574v1.
#' @seealso ApplyFwerCor, BonferroniCor_SD
#'
#' @examples 
#' n <- 100
#' p <- 10
#' corr_theo <- diag(1,p)
#' corr_theo[1,3] <- 0.5
#' corr_theo[3,1] <- 0.5
#' corr_theo <- diag(1,p)
#' data <- MASS::mvrnorm(n,rep(0,p),corr_theo)
#' # adjusted p-values
#' res <- BonferroniCor(data,stat_test='empirical')
#' round(res,2)
#' # significant correlations with level alpha:
#' alpha <- 0.05
#' whichCor(res<alpha)
#' # directly
#' BonferroniCor(data,alpha,stat_test='empirical',arr.ind=TRUE)
BonferroniCor <- function(data,alpha=0.05,stat_test='empirical',vect=FALSE,logical=FALSE,arr.ind=FALSE){

  if(sum(stat_test==c('empirical','fisher','student','2nd.order'))==0){ stop('Wrong value for stat_test.')}
  
	stat <- abs(eval_stat(data,stat_test))
    m <- length(stat)
    pval <- pmin(2*m*(1-pnorm(stat)), 1)

    if(arr.ind){ 
      vect <- FALSE 
      logical <- TRUE
    }
    
    res <- pval  
    if(!vect){ res <- unvectorize(pval) }
    if(logical){ res <- (res < alpha) }
    if(arr.ind){ res <- whichCor(res) }
    
    return(res)
  
}


#--------------------------------method 2 : Sidak-------------------------------------------
#--------------------------------------------------------------------------------------------

#' Sidak multiple testing procedure for correlations.
#'
#' @param data         matrix of observations
#' @param alpha        level of multiple testing (used if logical=TRUE)
#' @param stat_test  
#' \describe{
#'   \item{'empirical'}{\eqn{\sqrt{n}*abs(corr)}}
#'   \item{'fisher'}{   \eqn{\sqrt{n-3}*1/2*\log( (1+corr)/(1-corr) )}}
#'   \item{'student'}{  \eqn{\sqrt{n-2}*abs(corr)/\sqrt(1-corr^2)}}
#'   \item{'2nd.order'}{ \eqn{\sqrt{n}*mean(Y)/sd(Y)} with \eqn{Y=(X_i-mean(X_i))(X_j-mean(X_j))}}
#' }
#' @param vect         if TRUE returns a vector of adjusted p-values, corresponding to \code{vectorize(cor(data))};
#'                     if FALSE, returns an array containing the adjusted p-values for each entry of the correlation matrix 
#' @param logical      if TRUE, returns either a vector or a matrix where each element is equal to TRUE if the corresponding null hypothesis is rejected, and to FALSE if it is not rejected
#' @param arr.ind      if TRUE, returns the indexes of the significant correlations, with respect to level alpha
#'
#' @return Returns  \itemize{\item{the adjusted p-values, as a vector or a matrix depending of the value of \code{vect},} \item{an array containing indexes \eqn{\lbrace(i,j),\,i<j\rbrace} for which correlation between variables \eqn{i} and \eqn{j} is significant, if \code{arr.ind=TRUE}.}}
#'
#' @importFrom stats pnorm
#' @importFrom MASS mvrnorm
#' @export
#'
#' @references  Roux, M. (2018). Graph inference by multiple testing with application to Neuroimaging, Ph.D., Université Grenoble Alpes, France, https://tel.archives-ouvertes.fr/tel-01971574v1.
#' @references  Šidák, Z. (1967). Rectangular confidence regions for the means of multivariate normal distributions. Journal of the American Statistical Association, 62(318), 626-633.
#' @seealso ApplyFwerCor, SidakCor_SD
#'
#' @examples  
#' n <- 100
#' p <- 10
#' corr_theo <- diag(1,p)
#' corr_theo[1,3] <- 0.5
#' corr_theo[3,1] <- 0.5
#' data <- MASS::mvrnorm(n,rep(0,p),corr_theo)
#' # adjusted p-values
#' res <- SidakCor(data,stat_test='empirical')
#' round(res,2)
#' # significant correlations with level alpha:
#' alpha <- 0.05
#' whichCor(res<alpha)
#' # directly
#' SidakCor(data,alpha,stat_test='empirical',arr.ind=TRUE)
SidakCor <- function(data,alpha=0.05,stat_test='empirical',vect=FALSE,logical=FALSE,arr.ind=FALSE){

  if(sum(stat_test==c('empirical','fisher','student','2nd.order'))==0){ stop('Wrong value for stat_test.')}
  
	stat <- abs(eval_stat(data,stat_test))
    m <- length(stat)
    pval <- pmin(1-(2*pnorm(stat)-1)^m, 1)

    if(arr.ind){ 
      vect <- FALSE 
      logical <- TRUE
    }
    
    res <- pval  
    if(!vect){ res <- unvectorize(pval) }
    if(logical){ res <- (res < alpha) }
    if(arr.ind){ res <- whichCor(res) }
    
    return(res)
}


#--------------------------------method 3 : Bootstrap Romano/Wolf---------------------------
#--------------------------------------------------------------------------------------------

#' Bootstrap multiple testing method of Romano & Wolf (2005) for correlations.
#' @description Multiple testing method based on the evaluation of quantile by bootstrap 
#' in the initial dataset (Romano & Wolf (2005)).
#'
#' @param data         matrix of observations
#' @param alpha        level of multiple testing (used if logical=TRUE)
#' @param stat_test  
#' \describe{
#'   \item{'empirical'}{\eqn{\sqrt{n}*abs(corr)}}
#'   \item{'fisher'}{   \eqn{\sqrt{n-3}*1/2*\log( (1+corr)/(1-corr) )}}
#'   \item{'student'}{  \eqn{\sqrt{n-2}*abs(corr)/\sqrt(1-corr^2)}}
#'   \item{'2nd.order'}{ \eqn{\sqrt{n}*mean(Y)/sd(Y)} with \eqn{Y=(X_i-mean(X_i))(X_j-mean(X_j))}}
#' }
#' @param Nboot        number of iterations for Monte-Carlo quantile evaluation
#' @param vect         if TRUE returns a vector of adjusted p-values, corresponding to \code{vectorize(cor(data))};
#'                     if FALSE, returns an array containing the adjusted p-values for each entry of the correlation matrix 
#' @param logical      if TRUE, returns either a vector or a matrix where each element is equal to TRUE if the corresponding null hypothesis is rejected, and to FALSE if it is not rejected
#' @param arr.ind      if TRUE, returns the indexes of the significant correlations, with respect to level alpha
#'
#' @return Returns  \itemize{\item{the adjusted p-values, as a vector or a matrix depending of the value of \code{vect},} \item{an array containing indexes \eqn{\lbrace(i,j),\,i<j\rbrace} for which correlation between variables \eqn{i} and \eqn{j} is significant, if \code{arr.ind=TRUE}.}}
#'
#' @importFrom MASS mvrnorm
#' @export
#'
#' @references  Romano, J. P., & Wolf, M. (2005). Exact and approximate stepdown methods for multiple hypothesis testing. Journal of the American Statistical Association, 100(469), 94-108.
#' @references  Roux, M. (2018). Graph inference by multiple testing with application to Neuroimaging, Ph.D., Université Grenoble Alpes, France, https://tel.archives-ouvertes.fr/tel-01971574v1.
#' @seealso ApplyFwerCor, BootRWCor_SD
#'
#' @examples  
#' n <- 100
#' p <- 10
#' corr_theo <- diag(1,p)
#' corr_theo[1,3] <- 0.5
#' corr_theo[3,1] <- 0.5
#' data <- MASS::mvrnorm(n,rep(0,p),corr_theo)
#' # adjusted p-values
#' res <- BootRWCor(data,stat_test='empirical',Nboot=1000)
#' round(res,2)
#' # significant correlations with level alpha:
#' alpha <- 0.05
#' whichCor(res<alpha)
#' # directly
#' BootRWCor(data,alpha,stat_test='empirical',Nboot=1000,arr.ind=TRUE)
BootRWCor <- function(data,alpha=0.05,stat_test='empirical',Nboot=1000,vect=FALSE,logical=FALSE,arr.ind=FALSE){

  if(sum(stat_test==c('empirical','fisher','student','2nd.order'))==0){ stop('Wrong value for stat_test.')}
  
    # test statistic	
    n <- nrow(data)
    p <- ncol(data)
    stat <- eval_stat(data,stat_test)

    # evaluation of quantile by bootstrap
    max_boot <- rep(0,Nboot)	
 	for (nboot in 1:Nboot){
		
		indb <- sample(seq(1,n,1),replace=TRUE)
		datab <- data[indb,] 		
		stat_boot <- eval_stat(datab,stat_test)
        stat_boot <- abs(stat_boot-stat)
	
	  max_boot[nboot] <- max(stat_boot)
  }
	
  pval <- sapply(abs(stat),function(x){return(mean(max_boot>x))})

  if(arr.ind){ 
    vect <- FALSE 
    logical <- TRUE
  }
  
  res <- pval  
  if(!vect){ res <- unvectorize(pval) }
  if(logical){ res <- (res < alpha) }
  if(arr.ind){ res <- whichCor(res) }
  
  return(res)
}



#--------------------------------method 4 : maxT_infty--------------------------------------
#--------------------------------------------------------------------------------------------

#' Multiple testing method of Drton & Perlman (2007) for correlations.
#' @description Multiple testing method based on the evaluation of quantile by simulation of observations 
#' from the asymptotic distribution (Drton & Perlman (2007)).
#'
#' @param data         matrix of observations
#' @param alpha        level of multiple testing (used if logical=TRUE)
#' @param stat_test  
#' \describe{
#'   \item{'empirical'}{\eqn{\sqrt{n}*abs(corr)}}
#'   \item{'fisher'}{   \eqn{\sqrt{n-3}*1/2*\log( (1+corr)/(1-corr) )}}
#'   \item{'student'}{  \eqn{\sqrt{n-2}*abs(corr)/\sqrt(1-corr^2)}}
#'   \item{'2nd.order'}{ \eqn{\sqrt{n}*mean(Y)/sd(Y)} with \eqn{Y=(X_i-mean(X_i))(X_j-mean(X_j))}}
#' }
#' @param Nboot        number of iterations for Monte-Carlo quantile evaluation
#' @param OmegaChap    matrix of covariance of empirical correlations used for quantile evaluation;
#'                     optional, useful for oracle estimation and step-down
#' @param vect         if TRUE returns a vector of adjusted p-values, corresponding to \code{vectorize(cor(data))};
#'                     if FALSE, returns an array containing the adjusted p-values for each entry of the correlation matrix 
#' @param logical      if TRUE, returns either a vector or a matrix where each element is equal to TRUE if the corresponding null hypothesis is rejected, and to FALSE if it is not rejected
#' @param arr.ind      if TRUE, returns the indexes of the significant correlations, with respect to level alpha
#'
#' @return Returns  \itemize{\item{the adjusted p-values, as a vector or a matrix depending of the value of \code{vect},} \item{an array containing indexes \eqn{\lbrace(i,j),\,i<j\rbrace} for which correlation between variables \eqn{i} and \eqn{j} is significant, if \code{arr.ind=TRUE}.}}
#'
#' @importFrom MASS mvrnorm
#' @export
#'
#' @references  Drton, M., & Perlman, M. D. (2007). Multiple testing and error control in Gaussian graphical model selection. Statistical Science, 22(3), 430-449.
#' @references  Roux, M. (2018). Graph inference by multiple testing with application to Neuroimaging, Ph.D., Université Grenoble Alpes, France, https://tel.archives-ouvertes.fr/tel-01971574v1.
#' @seealso ApplyFwerCor, maxTinftyCor_SD
#'
#' @examples  
#' n <- 100
#' p <- 10
#' corr_theo <- diag(1,p)
#' corr_theo[1,3] <- 0.5
#' corr_theo[3,1] <- 0.5
#' data <- MASS::mvrnorm(n,rep(0,p),corr_theo)
#' # adjusted p-values
#' res <- maxTinftyCor(data,stat_test='empirical',Nboot=1000)
#' round(res,2)
#' # significant correlations with level alpha:
#' alpha <- 0.05
#' whichCor(res<alpha)
#' # directly
#' res <- maxTinftyCor(data,alpha,stat_test='empirical',Nboot=1000,arr.ind=TRUE)
maxTinftyCor <- function(data,alpha=0.05,stat_test='empirical',Nboot=1000,OmegaChap = covDcorNorm(cor(data),stat_test),vect=FALSE,logical=FALSE,arr.ind=FALSE){

  if(sum(stat_test==c('empirical','fisher','student','2nd.order'))==0){ stop('Wrong value for stat_test.')}
  
    stat <- abs(eval_stat(data,stat_test))
    OmegaChap <- as.matrix(OmegaChap)
 
    # evaluation of the (1-alpha/2)-quantile of a N(0,OmegaChap) by simulation
    dataq <- abs(mvrnorm(Nboot,rep(0,nrow(OmegaChap)),OmegaChap))
    maxq <- dataq[cbind(1:nrow(dataq),max.col(dataq))]
	
  pval <- sapply(abs(stat),function(x){return(mean(maxq>x))})

  if(arr.ind){ 
    vect <- FALSE 
    logical <- TRUE
  }
  
  res <- pval  
  if(!vect){ res <- unvectorize(pval) }
  if(logical){ res <- (res < alpha) }
  if(arr.ind){ res <- whichCor(res) }
  
  return(res)
  
}



