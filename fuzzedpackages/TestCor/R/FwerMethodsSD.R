# Multiple testing procedures controlling (asymptotically) the FWER
# for tests on a correlation matrix ** with Stepdown **
###############################################################################

#--------------------------------method 1 : Bonferroni--------------------------------------
#--------------------------------------------------------------------------------------------

#' Bonferroni multiple testing method for correlations 
#' with stepdown procedure.
#'
#' @param data       matrix of observations
#' @param alpha      level of multiple testing
#' @param stat_test  4 test statistics are available:
#' \describe{
#'   \item{'empirical'}{ \eqn{\sqrt{n}*abs(corr)}}
#'   \item{'fisher'}{   \eqn{\sqrt{n-3}*1/2*\log( (1+corr)/(1-corr) )}}
#'   \item{'student'}{  \eqn{\sqrt{n-2}*abs(corr)/\sqrt(1-corr^2)}}
#'   \item{'2nd.order'}{ \eqn{\sqrt{n}*mean(Y)/sd(Y)} with \eqn{Y=(X_i-mean(X_i))(X_j-mean(X_j))}}
#' }
#' @param vect      if TRUE returns a vector of TRUE/FALSE values, corresponding to \code{vectorize(cor(data))};
#'                  if FALSE, returns an array containing TRUE/FALSE values for each entry of the correlation matrix
#' @param logical   if TRUE, returns either a vector or a matrix where each element is equal to TRUE if the corresponding null hypothesis is rejected, and to FALSE if it is not rejected
#'                  if FALSE, returns a list of successive p-values : element [[i+1]] of the list giving the p-values evaluated on the non-rejected hypothesis at step [[i]]; p-values are either as a vector or a list depending on \code{vect}
#' @param arr.ind   if TRUE, returns the indexes of the significant correlations, with respect to level alpha
#'
#' @return Returns  \itemize{\item{logicals, equal to TRUE if the corresponding element of the statistic vector is rejected, as a vector or a matrix depending of the value of \code{vect},} \item{an array containing indexes \eqn{\lbrace(i,j),\,i<j\rbrace} for which correlation between variables \eqn{i} and \eqn{j} is significant, if \code{arr.ind=TRUE}.}}
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats cor
#' @importFrom stats qnorm
#' @export
#'
#' @references Bonferroni, C. E. (1935). Il calcolo delle assicurazioni su gruppi di teste. Studi in onore del professore salvatore ortu carboni, 13-60.
#' @references  Roux, M. (2018). Graph inference by multiple testing with application to Neuroimaging, Ph.D., Université Grenoble Alpes, France, https://tel.archives-ouvertes.fr/tel-01971574v1.
#' @seealso ApplyFwerCor, BonferroniCor
#'
#' @examples  
#' n <- 100
#' p <- 10
#' corr_theo <- diag(1,p)
#' corr_theo[1,3] <- 0.5
#' corr_theo[3,1] <- 0.5
#' data <- MASS::mvrnorm(n,rep(0,p),corr_theo)
#' alpha <- 0.05
#' # significant correlations:
#' BonferroniCor_SD(data,alpha,stat_test='empirical', arr.ind=TRUE)
#' # successive p-values
#' res <- BonferroniCor_SD(data,stat_test='empirical', logical=FALSE)
#' lapply(res,FUN=function(x){round(x,2)})
#' # succesive rejections
#' lapply(res,FUN=function(x){whichCor(x<alpha)})  
BonferroniCor_SD <- function(data,alpha=0.05,stat_test='empirical',vect=FALSE,logical=TRUE,arr.ind=FALSE){

  if(sum(stat_test==c('empirical','fisher','student','2nd.order'))==0){ stop('Wrong value for stat_test.')}
  
	stat <- abs(eval_stat(data, stat_test))
  m <- length(stat)
  pval <- pmin(2*m*(1-pnorm(stat)), 1)
	vBonf <- (pval < alpha)

	  pval <- list(pval)
	  res_SD <- rep(0, length(vBonf))
    indNR <- which(vBonf == FALSE)
    indR <- which(vBonf != FALSE)
    res_SD[indR] <- 1
	
  while((sum(vBonf)!=0)&&(sum(indNR)!=0)){
    stat_SD <- stat[indNR]
    m <- length(stat_SD)
    pval_SD <- pmin(2*m*(1-pnorm(stat_SD)), 1)
    pv <- vector(mode='numeric',length=m)
    pv[indNR] <- pval_SD
    pval <- c(pval, list(pv))

    vBonf <- (pval_SD < alpha)
	  indR <- which(vBonf !=0)            
	  res_SD[indNR[indR]] <- 1
	  indNR <- indNR[-indR]	
  }
	
    if(arr.ind){ 
      logical <- TRUE
      vect <- FALSE 
    }
    
    res <- (res_SD>0)
    if(!vect){ res <- (unvectorize(res_SD)>0) }
    if(arr.ind){ res <- whichCor(res) }
    if(!logical){ 
      res <- pval
      if(!vect){ res <- lapply(pval,unvectorize) }
    }
    
    return(res)
}


#--------------------------------method 2 : Sidak-------------------------------------------
#--------------------------------------------------------------------------------------------

#' Sidak multiple testing method for correlations
#' with stepdown procedure.
#'
#' @param data       matrix of observations
#' @param alpha      level of multiple testing
#' @param stat_test  4 test statistics are available:
#' \describe{
#'   \item{'empirical'}{ \eqn{\sqrt{n}*abs(corr)}}
#'   \item{'fisher'}{ \eqn{\sqrt{n-3}*1/2*\log( (1+corr)/(1-corr) )}}
#'   \item{'student'}{ \eqn{\sqrt{n-2}*abs(corr)/\sqrt(1-corr^2)}}
#'   \item{'2nd.order'}{ \eqn{\sqrt{n}*mean(Y)/sd(Y)} with \eqn{Y=(X_i-mean(X_i))(X_j-mean(X_j))}}
#' }
#' @param vect      if TRUE returns a vector of TRUE/FALSE values, corresponding to \code{vectorize(cor(data))};
#'                  if FALSE, returns an array containing TRUE/FALSE values for each entry of the correlation matrix
#' @param logical   if TRUE, returns either a vector or a matrix where each element is equal to TRUE if the corresponding null hypothesis is rejected, and to FALSE if it is not rejected
#'                  if FALSE, returns a list of successive p-values : element [[i+1]] of the list giving the p-values evaluated on the non-rejected hypothesis at step [[i]]; p-values are either as a vector or a list depending on \code{vect}
#' @param arr.ind   if TRUE, returns the indexes of the significant correlations, with respect to level alpha
#'
#' @return Returns  \itemize{\item{logicals, equal to TRUE if the corresponding element of the statistic vector is rejected, as a vector or a matrix depending of the value of \code{vect},} \item{an array containing indexes \eqn{\lbrace(i,j),\,i<j\rbrace} for which correlation between variables \eqn{i} and \eqn{j} is significant, if \code{arr.ind=TRUE}.}}
#'
#' @export
#' @importFrom MASS mvrnorm
#' @importFrom stats cor
#' @importFrom stats qnorm
#'
#' @references  Roux, M. (2018). Graph inference by multiple testing with application to Neuroimaging, Ph.D., Université Grenoble Alpes, France, https://tel.archives-ouvertes.fr/tel-01971574v1.
#' @references Šidák, Z. (1967). Rectangular confidence regions for the means of multivariate normal distributions. Journal of the American Statistical Association, 62(318), 626-633.
#' @seealso ApplyFwerCor, SidakCor
#'
#' @examples   
#' n <- 100
#' p <- 10
#' corr_theo <- diag(1,p)
#' corr_theo[1,3] <- 0.5
#' corr_theo[3,1] <- 0.5
#' data <- MASS::mvrnorm(n,rep(0,p),corr_theo)
#' alpha <- 0.05
#' # significant correlations:
#' SidakCor_SD(data,alpha,stat_test='empirical', arr.ind=TRUE)
#' # successive p-values
#' res <- SidakCor_SD(data,stat_test='empirical', logical=FALSE)
#' lapply(res,FUN=function(x){round(x,2)})
#' # succesive rejections
#' lapply(res,FUN=function(x){whichCor(x<alpha)})  
SidakCor_SD <- function(data,alpha=0.05,stat_test='empirical',vect=FALSE,logical=TRUE,arr.ind=FALSE){

  if(sum(stat_test==c('empirical','fisher','student','2nd.order'))==0){ stop('Wrong value for stat_test.')}
  
	stat <- abs(eval_stat(data,stat_test))
  m <- length(stat)
  pval <- pmin(1-(2*pnorm(stat)-1)^m, 1)
  vSidak  <- (pval < alpha)

  pval <- list(pval)
  res_SD <- rep(0,length(vSidak))
  indNR <- which(vSidak == FALSE)   # indexes of non rejected hypothesis 
  indR <- which(vSidak != FALSE)               
  res_SD[indR] <- 1
	
  while((sum(vSidak)!=0)&&(sum(indNR)!=0)){
      stat_SD <- stat[indNR]
      m <- length(stat_SD)
      pval_SD <- pmin(1-(2*pnorm(stat_SD)-1)^m, 1)
      pv <- vector(mode='numeric',length=m)
      pv[indNR] <- pval_SD
      pval <- c(pval, list(pv))
      
    vSidak  <- (pval_SD < alpha)
	  indR <- which(vSidak !=0)            
	  res_SD[indNR[indR]] <- 1
	  indNR =indNR[-indR]
  }
	
  if(arr.ind){ 
    logical <- TRUE
    vect <- FALSE 
  }
  
  res <- (res_SD>0)
  if(!vect){ res <- (unvectorize(res_SD)>0) }
  if(arr.ind){ res <- whichCor(res) }
  if(!logical){ 
    res <- pval
    if(!vect){ res <- lapply(pval,unvectorize) }
  }
  
  return(res)
}



#--------------------------------method 3 : Bootstrap Romano/Wolf---------------------------
#--------------------------------------------------------------------------------------------

#' Boootstrap multiple testing method of Romano & Wolf (2005) for correlations, with stepdown procedure.
#'
#' @description  Multiple testing method based on the evaluation of quantile by bootstrap 
#' in the initial dataset (Romano & Wolf (2005)),
#' with stepdown procedure.
#'
#' @param data      matrix of observations
#' @param alpha     level of multiple testing
#' @param stat_test  4 test statistics are available:
#' \describe{
#'   \item{'empirical'}{\eqn{\sqrt{n}*abs(corr)}}
#'   \item{'fisher'}{   \eqn{\sqrt{n-3}*1/2*\log( (1+corr)/(1-corr) )}}
#'   \item{'student'}{  \eqn{\sqrt{n-2}*abs(corr)/\sqrt(1-corr^2)}}
#'   \item{'2nd.order'}{ \eqn{\sqrt{n}*mean(Y)/sd(Y)} with \eqn{Y=(X_i-mean(X_i))(X_j-mean(X_j))}}
#' }
#' @param Nboot      number of iterations for Bootstrap quantile evaluation
#' @param vect      if TRUE returns a vector of TRUE/FALSE values, corresponding to \code{vectorize(cor(data))};
#'                  if FALSE, returns an array containing TRUE/FALSE values for each entry of the correlation matrix
#' @param logical   if TRUE, returns either a vector or a matrix where each element is equal to TRUE if the corresponding null hypothesis is rejected, and to FALSE if it is not rejected
#'                  if FALSE, returns a list of successive p-values : element [[i+1]] of the list giving the p-values evaluated on the non-rejected hypothesis at step [[i]]; p-values are either as a vector or a list depending on \code{vect}
#' @param arr.ind   if TRUE, returns the indexes of the significant correlations, with respect to level alpha
#'
#' @return Returns  \itemize{\item{logicals, equal to TRUE if the corresponding element of the statistic vector is rejected, as a vector or a matrix depending of the value of \code{vect},} \item{an array containing indexes \eqn{\lbrace(i,j),\,i<j\rbrace} for which correlation between variables \eqn{i} and \eqn{j} is significant, if \code{arr.ind=TRUE}.}}
#' 
#' @importFrom stats cor
#' @importFrom MASS mvrnorm
#' @export
#'
#' @references Romano, J. P., & Wolf, M. (2005). Exact and approximate stepdown methods for multiple hypothesis testing. Journal of the American Statistical Association, 100(469), 94-108.
#' @references  Roux, M. (2018). Graph inference by multiple testing with application to Neuroimaging, Ph.D., Université Grenoble Alpes, France, https://tel.archives-ouvertes.fr/tel-01971574v1.
#' @seealso ApplyFwerCor, BootRWCor
#' 
#' @examples  
#' n <- 100
#' p <- 10
#' corr_theo <- diag(1,p)
#' corr_theo[1,3] <- 0.5
#' corr_theo[3,1] <- 0.5
#' data <- MASS::mvrnorm(n,rep(0,p),corr_theo)
#' alpha <- 0.05
#' # significant correlations:
#' BootRWCor_SD(data,alpha,stat_test='empirical', arr.ind=TRUE)
#' # successive p-values
#' res <- BootRWCor_SD(data,stat_test='empirical', logical=FALSE)
#' lapply(res,FUN=function(x){round(x,2)})
#' # succesive rejections
#' lapply(res,FUN=function(x){whichCor(x<alpha)})  
BootRWCor_SD <- function(data,alpha=0.05,stat_test='empirical',Nboot=1000,vect=FALSE,logical=TRUE,arr.ind=FALSE){

  if(sum(stat_test==c('empirical','fisher','student','2nd.order'))==0){ stop('Wrong value for stat_test.')}
  
  pval <-  BootRWCor(data,stat_test=stat_test,Nboot=Nboot,vect=TRUE,logical=FALSE)
  m <- length(pval)
  vBootRW <- (pval < alpha)

  pval <- list(pval)
  res_SD <- rep(0,length(vBootRW))
  indNR <- which(vBootRW == FALSE)   # indexes of non rejected hypothesis 
  indR <- which(vBootRW != FALSE)               
  res_SD[indR] <- 1
	
  n <- ncol(data)
  stat <- eval_stat(data,stat_test)
   
  while((sum(vBootRW)!=0)&&(sum(indNR)!=0)){
		
    # test statistic
    stat_SD <- stat[indNR]

    # evaluation of quantile by bootstrap
    max_boot <- rep(0,Nboot)

    for (nboot in 1:Nboot){
     indb <- sample(seq(1,n,1),replace=TRUE)
	   data_boot <- data[indb,]
	   stat_boot <- eval_stat(data_boot,stat_test)

     stat_boot <- abs(stat_boot[indNR]-stat_SD)
     max_boot[nboot] <- max(stat_boot)
    }
   pval_SD <- sapply(abs(stat_SD),function(x){return(mean(max_boot>x))})
   pv <- vector(mode='numeric',length=m)
   pv[indNR] <- pval_SD
   pval <- c(pval, list(pv))
  
   vBootRW <- (pval_SD < alpha)
   indR <- which(vBootRW !=0)            
   res_SD[indNR[indR]] <- 1
   indNR <- indNR[-indR]
	
  }
	
  if(arr.ind){ 
    logical <- TRUE
    vect <- FALSE 
  }
  
  res <- (res_SD>0)
  if(!vect){ res <- (unvectorize(res_SD)>0) }
  if(arr.ind){ res <- whichCor(res) }
  if(!logical){ 
    res <- pval
    if(!vect){ res <- lapply(pval,unvectorize) }
  }
  
  return(res)
	
}



#--------------------------------method 4 : maxT_infty--------------------------------------
#--------------------------------------------------------------------------------------------
#' Multiple testing method of Drton & Perlman (2007) for correlations, with stepdown procedure.
#'
#' @description Multiple testing method based on the evaluation of quantile by simulation of observations 
#' from the asymptotic distribution (Drton & Perlman (2007)), 
#' with stepdown procedure.
#'
#' @param data         matrix of observations
#' @param alpha        level of multiple testing
#' @param stat_test    4 test statistics are available:
#' \describe{
#'   \item{'empirical'}{ \eqn{\sqrt{n}*abs(corr)}}
#'   \item{'fisher'}{   \eqn{\sqrt{n-3}*1/2*\log( (1+corr)/(1-corr) )}}
#'   \item{'student'}{  \eqn{\sqrt{n-2}*abs(corr)/\sqrt(1-corr^2)}}
#'   \item{'2nd.order'}{ \eqn{\sqrt{n}*mean(Y)/sd(Y)} with \eqn{Y=(X_i-mean(X_i))(X_j-mean(X_j))}}
#' }
#' @param Nboot        number of iterations for Monte-Carlo quantile evaluation
#' @param OmegaChap    matrix of covariance of test statistics;
#'                     optional, useful for oracle estimation and step-down
#' @param vect         if TRUE returns a vector of TRUE/FALSE values, corresponding to \code{vectorize(cor(data))};
#'                     if FALSE, returns an array containing TRUE/FALSE values for each entry of the correlation matrix
#' @param logical   if TRUE, returns either a vector or a matrix where each element is equal to TRUE if the corresponding null hypothesis is rejected, and to FALSE if it is not rejected
#'                  if FALSE, returns a list of successive p-values : element [[i+1]] of the list giving the p-values evaluated on the non-rejected hypothesis at step [[i]]; p-values are either as a vector or a list depending on \code{vect}
#' @param arr.ind      if TRUE, returns the indexes of the significant correlations, with respect to level alpha
#' 
#' @return Returns  \itemize{\item{logicals, equal to TRUE if the corresponding element of the statistic vector is rejected, as a vector or a matrix depending of the value of \code{vect},} \item{an array containing indexes \eqn{\lbrace(i,j),\,i<j\rbrace} for which correlation between variables \eqn{i} and \eqn{j} is significant, if \code{arr.ind=TRUE}.}}
#'
#' @export                     
#' @importFrom MASS mvrnorm
#' @importFrom stats cor
#'
#' @references Drton, M., & Perlman, M. D. (2007). Multiple testing and error control in Gaussian graphical model selection. Statistical Science, 22(3), 430-449.
#' @references  Roux, M. (2018). Graph inference by multiple testing with application to Neuroimaging, Ph.D., Université Grenoble Alpes, France, https://tel.archives-ouvertes.fr/tel-01971574v1.
#' @seealso ApplyFwerCor, maxTinftyCor
#'
#' @examples  
#' n <- 100
#' p <- 10
#' corr_theo <- diag(1,p)
#' corr_theo[1,3] <- 0.5
#' corr_theo[3,1] <- 0.5
#' data <- MASS::mvrnorm(n,rep(0,p),corr_theo)
#' alpha <- 0.05
#' # significant correlations:
#' maxTinftyCor_SD(data,alpha,stat_test='empirical', arr.ind=TRUE)
#' # successive p-values
#' res <- maxTinftyCor_SD(data,stat_test='empirical', logical=FALSE)
#' lapply(res,FUN=function(x){round(x,2)})
#' # succesive rejections
#' lapply(res,FUN=function(x){whichCor(x<alpha)})  
maxTinftyCor_SD <- function(data,alpha=0.05,stat_test='empirical',Nboot=1000,OmegaChap=covDcorNorm(cor(data),stat_test),vect=FALSE,logical=TRUE,arr.ind=FALSE){ 

  if(sum(stat_test==c('empirical','fisher','student','2nd.order'))==0){ stop('Wrong value for stat_test.')}
  
  pval <- maxTinftyCor(data,stat_test=stat_test,Nboot=Nboot,OmegaChap=OmegaChap,vect=TRUE,logical=FALSE)
  m <- length(pval)
  vmaxTinf <- (pval < alpha)
  stat <- abs(eval_stat(data,stat_test))

  pval <- list(pval)
  res_SD <- rep(0,length(vmaxTinf))
  indNR <- which(vmaxTinf == 0)   # indexes of non rejected hypothesis 
  indR <- which(vmaxTinf != 0)               
  res_SD[indR] <- 1
	
  while((sum(vmaxTinf)!=0)&&(sum(indNR)!=0)){
    stat_SD <- stat[indNR]
    OmegaChap_SD <- OmegaChap[indNR,indNR]
 
    # evaluation of the (1-alpha/2)-quantile of a N(0,OmegaChap) by simulation
    dataq <- abs(mvrnorm(Nboot,rep(0,nrow(OmegaChap_SD)),OmegaChap_SD))
    maxq <- dataq[cbind(1:nrow(dataq),max.col(dataq))]
    pval_SD <- sapply(abs(stat_SD),function(x){return(mean(maxq>x))})
    pv <- vector(mode='numeric',length=m)
    pv[indNR] <- pval_SD
    pval <- c(pval, list(pv))
    
    vmaxTinf <- (pval_SD < alpha)
    indR <- which(vmaxTinf !=0)            
    res_SD[indNR[indR]] <- 1
    indNR <- indNR[-indR]
  }

  if(arr.ind){ 
    logical <- TRUE
    vect <- FALSE 
  }
  
  res <- (res_SD>0)
  if(!vect){ res <- (unvectorize(res_SD)>0) }
  if(arr.ind){ res <- whichCor(res) }
  if(!logical){ 
    res <- pval
    if(!vect){ res <- lapply(pval,unvectorize) }
  }
  
  return(res)
	
}



