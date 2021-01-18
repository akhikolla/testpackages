#--------------------------------method 1 : Cai et Liu -------------------------------------------
#-------------------------------------------------------------------------------------------------

#' Procedure LCT-N proposed by Cai & Liu (2016) for correlation testing.
#'
#' @param data         matrix of observations
#' @param alpha        level of multiple testing
#' @param stat_test     
#' \describe{
#'   \item{'empirical'}{\eqn{\sqrt{n}*abs(corr)}}
#'   \item{'fisher'}{   \eqn{\sqrt{n-3}*1/2*\log( (1+corr)/(1-corr) )}}
#'   \item{'student'}{  \eqn{\sqrt{n-2}*abs(corr)/\sqrt(1-corr^2)}}
#'   \item{'2nd.order'}{ \eqn{\sqrt{n}*mean(Y)/sd(Y)} with \eqn{Y=(X_i-mean(X_i))(X_j-mean(X_j))}}
#' }
#' @param vect      if TRUE returns a vector of TRUE/FALSE values, corresponding to \code{vectorize(cor(data))};
#'                  if FALSE, returns an array containing TRUE/FALSE values for each entry of the correlation matrix
#' @param arr.ind   if TRUE, returns the indexes of the significant correlations, with respect to level alpha
#'
#' @return Returns  \itemize{\item{logicals, equal to TRUE if the corresponding element of the statistic vector is rejected, as a vector or a matrix depending of the value of \code{vect},} \item{an array containing indexes \eqn{\lbrace(i,j),\,i<j\rbrace} for which correlation between variables \eqn{i} and \eqn{j} is significant, if \code{arr.ind=TRUE}.}}
#'
#' @importFrom stats cor quantile var
#' 
#' @export
#' @seealso ApplyFdrCor, LCTboot
#' @references  Cai, T. T., & Liu, W. (2016). Large-scale multiple testing of correlations. Journal of the American Statistical Association, 111(513), 229-240.
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
#' LCTnorm(data,alpha,stat_test='empirical',arr.ind=TRUE)
LCTnorm <- function(data,alpha=0.05,stat_test='2nd.order',vect=FALSE,arr.ind=FALSE){

  if(sum(stat_test==c('empirical','fisher','student','2nd.order'))==0){ stop('Wrong value for stat_test.')}
  
  p <- ncol(data)	

  # test statistic	
 	stat <- abs(eval_stat(data,stat_test))
	m <- length(stat)
  stat_sort <- sort(stat,index.return=TRUE,decreasing=TRUE)
  
	# evaluation of pval by normal approximation
	pval <- 2*(1-pnorm(stat_sort$x))
     
    # truncated BH threshold   
    ind <- which( ( pval < alpha*seq(1,m,1)/m )&( stat_sort$x <= sqrt(4*log(p)-2*log(log(p))) ) )
    t <- ifelse(length(ind)==0, 2*sqrt(log(p)), stat_sort$x[max(ind)] )

	res <- (stat >= t) #################### 

	if(arr.ind){ vect <- FALSE }
	if(!vect){ res <- (unvectorize(res)>0) }
	if(arr.ind){ res <- whichCor(res) }
	
	return(res)

}


#--------------------------------method 2 : Cai et Liu with bootstrap ----------------------------
#-------------------------------------------------------------------------------------------------

#' Bootstrap procedure LCT-B proposed by Cai & Liu (2016) for correlation testing.
#'
#' @param data          matrix of observations
#' @param alpha         level of multiple testing
#' @param stat_test     
#' \describe{
#'   \item{'empirical'}{\eqn{\sqrt{n}*abs(corr)}}
#'   \item{'fisher'}{\eqn{\sqrt{n-3}*1/2*\log( (1+corr)/(1-corr) )}}
#'   \item{'student'}{\eqn{\sqrt{n-2}*abs(corr)/\sqrt(1-corr^2)}}
#'   \item{'2nd.order'}{\eqn{\sqrt{n}*mean(Y)/sd(Y)} with \eqn{Y=(X_i-mean(X_i))(X_j-mean(X_j))}}
#' }
#' @param Nboot         number of iterations for bootstrap quantile evaluation
#' @param vect          if TRUE returns a vector of TRUE/FALSE values, corresponding to \code{vectorize(cor(data))};
#'                      if FALSE, returns an array containing TRUE/FALSE values for each entry of the correlation matrix
#' @param arr.ind   if TRUE, returns the indexes of the significant correlations, with respect to level alpha
#'
#' @return Returns  \itemize{\item{logicals, equal to TRUE if the corresponding element of the statistic vector is rejected, as a vector or a matrix depending of the value of \code{vect},} \item{an array containing indexes \eqn{\lbrace(i,j),\,i<j\rbrace} for which correlation between variables \eqn{i} and \eqn{j} is significant, if \code{arr.ind=TRUE}.}}
#'
#' @importFrom stats cor quantile var
#' 
#' @export
#' @seealso ApplyFdrCor, LCTNorm
#' @references  Cai, T. T., & Liu, W. (2016). Large-scale multiple testing of correlations. Journal of the American Statistical Association, 111(513), 229-240.
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
#' LCTboot(data,alpha,stat_test='empirical',Nboot=100,arr.ind=TRUE)
LCTboot <- function(data,alpha=0.05,stat_test='2nd.order',Nboot=100,vect=FALSE,arr.ind=FALSE){
 
  if(sum(stat_test==c('empirical','fisher','student','2nd.order'))==0){ stop('Wrong value for stat_test.')}
  
    n <- nrow(data)
    p <- ncol(data)
	
  # test statistic	
 	stat <- abs(eval_stat(data,stat_test))
	m <- length(stat)
  stat_sort <- sort(stat,index.return=TRUE,decreasing=TRUE)

	# evaluation of pval by bootstrap
	prop <- matrix(0,nrow=Nboot,ncol=m)	
	for (nboot in 1:Nboot){
		
		indb <- sample(seq(1,n,1),replace=TRUE)
		datab <- data[indb,]
		statb <- abs(eval_stat(datab,stat_test))
		
		# comparison with stat test
        for(i in 1:m){
		  prop[nboot,i] <- mean(statb > stat_sort$x[i]) } 
	}
	pval_boot <- colMeans(prop)

    
    # truncated BH threshold   
    ind <- which( ( pval_boot <= alpha*seq(1,m,1)/m )&( stat_sort$x <= sqrt(4*log(p)-2*log(log(p))) ) )
    t_boot <- ifelse(length(ind)==0, 2*sqrt(log(p)), stat_sort$x[max(ind)])

	res <- (stat >= t_boot)
	
	if(arr.ind){ vect <- FALSE }
	if(!vect){ res <- (unvectorize(res)>0) }
	if(arr.ind){ res <- whichCor(res) }
	
	return(res)
	
}


#--------------------------------method 3 : BH-----------------------------------------------
#----------------No theoretical proof for correlation tests----------------------------------
#--------------------------------------------------------------------------------------------

#' Benjamini & Hochberg (1995)'s procedure for correlation testing.
#' @description Benjamini & Hochberg (1995)'s procedure on the correlation matrix entries (no theoretical proof of control).
#'
#' @param data         matrix of observations
#' @param alpha        level of multiple testing
#' @param stat_test     
#' \describe{
#'   \item{'empirical'}{\eqn{\sqrt{n}*abs(corr)}}
#'   \item{'fisher'}{   \eqn{\sqrt{n-3}*1/2*\log( (1+corr)/(1-corr) )}}
#'   \item{'student'}{  \eqn{\sqrt{n-2}*abs(corr)/\sqrt(1-corr^2)}}
#'   \item{'2nd.order'}{ \eqn{\sqrt{n}*mean(Y)/sd(Y)} with \eqn{Y=(X_i-mean(X_i))(X_j-mean(X_j))}}
#' }
#' @param vect      if TRUE returns a vector of TRUE/FALSE values, corresponding to vectorize(cor(data))
#'                  if FALSE, returns an array containing TRUE/FALSE values for each entry of the correlation matrix
#' @param arr.ind   if TRUE, returns the indexes of the significant correlations, with respect to level alpha
#'
#' @return Returns  \itemize{\item{logicals, equal to TRUE if the corresponding element of the statistic vector is rejected, as a vector or a matrix depending of the value of \code{vect},} \item{an array containing indexes \eqn{\lbrace(i,j),\,i<j\rbrace} for which correlation between variables \eqn{i} and \eqn{j} is significant, if \code{arr.ind=TRUE}.}}
#'
#' @importFrom stats cor quantile pnorm
#' 
#' @export
#' @seealso ApplyFdrCor, BHBootCor
#' @references  Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the royal statistical society. Series B (Methodological), 289-300.
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
#' BHCor(data,alpha,stat_test='empirical',arr.ind=TRUE)
BHCor <- function(data,alpha=0.05,stat_test='2nd.order',vect=FALSE,arr.ind=FALSE){

  if(sum(stat_test==c('empirical','fisher','student','2nd.order'))==0){ stop('Wrong value for stat_test.')}
  
    p <- ncol(data)	

  # test statistic	
 	stat <- abs(eval_stat(data,stat_test))
  m <- length(stat)
  pval <- 2*(1-pnorm(stat))
  thresholds_BH <- alpha*seq(1,m,1)/m

    pval_sort <- sort(pval,index.return=TRUE)
    result <- (pval_sort$x <= thresholds_BH)
    if(sum(result) != 0){
	  ind_max <- max(which(result==1))
	  result <- rep(0,m)
	  result[pval_sort$ix[1:ind_max]] <- 1
    }
    
    res <- (result>0)
    if(arr.ind){ vect <- FALSE }
    if(!vect){ res <- (unvectorize(result)>0) }
    if(arr.ind){ res <- whichCor(res) }
    
    return(res)
    
}



#--------------------------------method 3 : BH-boot------------------------------------------
#----------------No theoretical proof for correlation tests----------------------------------
#--------------------------------------------------------------------------------------------

#' Benjamini & Hochberg (1995)'s procedure for correlation testing with bootstrap evaluation of p-values.
#' @description Benjamini & Hochberg (1995)'s procedure on the correlation matrix entries with bootstrap evaluation of p-values (no theoretical proof of control).
#'
#' @param data      matrix of observations
#' @param alpha      level of multiple testing
#' @param stat_test     
#' \describe{
#'   \item{'empirical'}{\eqn{\sqrt{n}*abs(corr)}}
#'   \item{'fisher'}{   \eqn{\sqrt{n-3}*1/2*\log( (1+corr)/(1-corr) )}}
#'   \item{'student'}{  \eqn{\sqrt{n-2}*abs(corr)/\sqrt(1-corr^2)}}
#'   \item{'2nd.order'}{ \eqn{\sqrt{n}*mean(Y)/sd(Y)} with \eqn{Y=(X_i-mean(X_i))(X_j-mean(X_j))}}
#' }
#' @param Nboot     number of iterations for bootstrap quantile evaluation
#' @param vect      if TRUE returns a vector of TRUE/FALSE values, corresponding to \code{vectorize(cor(data))};
#'                  if FALSE, returns an array containing TRUE/FALSE values for each entry of the correlation matrix
#' @param arr.ind   if TRUE, returns the indexes of the significant correlations, with respect to level alpha
#'
#' @return Returns \itemize{\item{a vector or a matrix of logicals, equal to TRUE if the corresponding element of the statistic vector is rejected, if \code{arr.ind=FALSE},} \item{an array containing indexes \eqn{\lbrace(i,j),\,i<j\rbrace} for which correlation between variables \eqn{i} and \eqn{j} is significant, if \code{arr.ind=TRUE}.}}
#'
#' @importFrom stats cor quantile pnorm
#' 
#' @export
#' @seealso ApplyFdrCor, BHCor
#' @references  Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the royal statistical society. Series B (Methodological), 289-300.
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
#' BHBootCor(data,alpha,stat_test='empirical',arr.ind=TRUE)
BHBootCor <- function(data,alpha=0.05,stat_test='2nd.order',Nboot=100,vect=FALSE,arr.ind=FALSE){
 	
  if(sum(stat_test==c('empirical','fisher','student','2nd.order'))==0){ stop('Wrong value for stat_test.')}
  
    n <- nrow(data)
    p <- ncol(data)
	
  # test statistic	
 	stat <- abs(eval_stat(data,stat_test))
	m <- length(stat)
  stat_sort <- sort(stat,index.return=TRUE,decreasing=TRUE)

	# evaluation of pval by bootstrap
	prop <- matrix(0,nrow=Nboot,ncol=m)	
	for (nboot in 1:Nboot){
		
		indb <- sample(seq(1,n,1),replace=TRUE)
		datab <- data[indb,]
		statb <- abs(eval_stat(datab,stat_test))
		
		# comparison with stat test
    for(i in 1:m){
		  prop[nboot,i] <- mean(statb > stat_sort$x[i]) } 
	  }
	  pval_boot <- colMeans(prop)

    thresholds_BH <- alpha*seq(1,m,1)/m
    pval_sort <- sort(pval_boot,index.return=TRUE)
    result <- (pval_sort$x <= thresholds_BH)
    if(sum(result) != 0){
	  ind_max <- max(which(result==1))
	  result <- rep(0,m)
	  result[pval_sort$ix[1:ind_max]] <- 1
    }
    
    res <- (result>0)
    if(arr.ind){ vect <- FALSE }
    if(!vect){ res <- (unvectorize(res)>0) }
    if(arr.ind){ res <- whichCor(res) }
    
    return(res)

}


