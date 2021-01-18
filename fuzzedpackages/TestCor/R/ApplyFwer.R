#' Applies multiple testing procedures controlling (asymptotically) the FWER
#' for tests on a correlation matrix.
#' @description Applies multiple testing procedures controlling (asymptotically) the FWER
#' for tests on a correlation matrix.
#' Methods are described in Chapter 5 of \cite{Roux (2018)}.
#'
#'
#' @return  Returns either \itemize{\item{the adjusted p-values, as a vector or a matrix, depending on \code{vect}} \item{logicals indicating if the corresponding correlation is significant if \code{logical=TRUE}, as a vector or a matrix depending on \code{vect},} \item{an array containing indexes \eqn{\lbrace(i,j),\,i<j\rbrace} for which correlation between variables \eqn{i} and \eqn{j} is significant, if \code{arr.ind=TRUE}.}}
#
#' @param data        matrix of observations
#' @param alpha       level of multiple testing (used if logical=TRUE)
#' @param stat_test  
#' \describe{
#'   \item{'empirical'}{\eqn{\sqrt{n}*abs(corr)}}
#'   \item{'fisher'}{   \eqn{\sqrt{n-3}*1/2*\log( (1+corr)/(1-corr) )}}
#'   \item{'student'}{  \eqn{\sqrt{n-2}*abs(corr)/\sqrt(1-corr^2)}}
#'   \item{'2nd.order'}{ \eqn{\sqrt{n}*mean(Y)/sd(Y)} with \eqn{Y=(X_i-mean(X_i))(X_j-mean(X_j))}}
#' }
#' @param method      choice between 'Bonferroni', 'Sidak', 'BootRW', 'MaxTinfty'
#' @param Nboot       number of iterations for Monte-Carlo of bootstrap quantile evaluation
#' @param stepdown    logical, if TRUE a stepdown procedure is applied
#' @param vect         if TRUE returns a vector of adjusted p-values, corresponding to \code{vectorize(cor(data))};
#'                     if FALSE, returns an array containing the adjusted p-values for each entry of the correlation matrix 
#' @param logical      if TRUE, returns either a vector or a matrix where each element is equal to TRUE if the corresponding null hypothesis is rejected, and to FALSE if it is not rejected
#'                     if \code{stepdown=TRUE} and \code{logical=FALSE}, returns a list of successive p-values.
#' @param arr.ind      if TRUE, returns the indexes of the significant correlations, with repspect to level alpha
#'
#' @importFrom stats cor
#' @importFrom MASS mvrnorm
#' @export
#'
#' @references Bonferroni, C. E. (1935). Il calcolo delle assicurazioni su gruppi di teste. Studi in onore del professore salvatore ortu carboni, 13-60.
#' @references  Drton, M., & Perlman, M. D. (2007). Multiple testing and error control in Gaussian graphical model selection. Statistical Science, 22(3), 430-449.
#' @references  Romano, J. P., & Wolf, M. (2005). Exact and approximate stepdown methods for multiple hypothesis testing. Journal of the American Statistical Association, 100(469), 94-108.
#' @references  Roux, M. (2018). Graph inference by multiple testing with application to Neuroimaging, Ph.D., Université Grenoble Alpes, France, https://tel.archives-ouvertes.fr/tel-01971574v1.
#' @references  Šidák, Z. (1967). Rectangular confidence regions for the means of multivariate normal distributions. Journal of the American Statistical Association, 62(318), 626-633.
#' @seealso ApplyFwerCor_SD, ApplyFdrCor
#' @seealso BonferroniCor, SidakCor, BootRWCor, maxTinftyCor
#' @seealso BonferroniCor_SD, SidakCor_SD, BootRWCor_SD, maxTinftyCor_SD 
#'
#' @examples
#' n <- 100
#' p <- 10
#' corr_theo <- diag(1,p)
#' corr_theo[1,3] <- 0.5
#' corr_theo[3,1] <- 0.5
#' data <- MASS::mvrnorm(n,rep(0,p),corr_theo)
#' # adjusted p-values
#' (res <- ApplyFwerCor(data,stat_test='empirical',method='Bonferroni',stepdown=FALSE))
#' # significant correlations, level alpha:
#' alpha <- 0.05
#' whichCor(res<alpha)
ApplyFwerCor <- function(data,alpha=NULL,stat_test='empirical',method='Sidak',Nboot=1000,stepdown=TRUE,vect=FALSE,logical=stepdown,arr.ind=FALSE){
 
 if(sum(stat_test==c('empirical','fisher','student','2nd.order'))==0){ stop('Wrong value for stat_test.')}
 if(sum(method==c('Bonferroni','Sidak','BootRW','MaxTinfty'))==0){ stop('Wrong value for method.')}
 # if( (!logical) & stepdown){
 #    warning('Stepdown procedures do not return p-values, logical is changed in TRUE.\n')
 #    logical = TRUE
 # }
 if(is.null(alpha) & (logical | arr.ind)){
   warning('For logical or indexes returns, a value of alpha is needed. alpha=0.05 is taken.\n')
   alpha = 0.05
 }
 if(arr.ind){ 
   vect <- FALSE
   logical <- TRUE 
 }
  
 if(method=='Bonferroni'){
	if(stepdown==FALSE){ res <- BonferroniCor(data,alpha,stat_test,vect,logical) }
	else{ res <- BonferroniCor_SD(data,alpha,stat_test,vect,logical=logical) }
 } 
 if(method=='Sidak'){
	if(stepdown==FALSE){ res <- SidakCor(data,alpha,stat_test,vect,logical) }
	else{ res <- SidakCor_SD(data,alpha,stat_test,vect,logical=logical) }
 }
 if(method=='BootRW'){
	if(stepdown==FALSE){ res <- BootRWCor(data,alpha,stat_test,Nboot,vect,logical) }
	else{ res <- BootRWCor_SD(data,alpha,stat_test,Nboot,vect=vect,logical=logical) }
 }
 if(method=='MaxTinfty'){
    if(stepdown==FALSE){ res <- maxTinftyCor(data,alpha,stat_test,Nboot=Nboot,vect=vect,logical=logical) }
    else{ res <- maxTinftyCor_SD(data,alpha,stat_test,Nboot,vect=vect,logical=logical) }
 }
  
 if(arr.ind){ res <- whichCor(res) }
  
 return(res)

}
 
#' Applies an oracle version of MaxTinfty procedure described in Drton & Perlman (2007) for correlation testing.             
#' @description Applies oracle MaxTinfty procedure described in Drton & Perlman (2007) which controls asymptotically the FWER
#' for tests on a correlation matrix. It needs the true correlation matrix. 
#'
#' @return  Returns either \itemize{\item{the adjusted p-values, as a vector or a matrix, depending on \code{vect} (unavailable with stepdown)} \item{logicals indicating if the corresponding correlation is significant if \code{logical=TRUE}, as a vector or a matrix depending on \code{vect},} \item{an array containing indexes \eqn{\lbrace(i,j),\,i<j\rbrace} for which correlation between variables \eqn{i} and \eqn{j} is significant, if \code{arr.ind=TRUE}.}}
#' Oracle estimation of the quantile is used, based on the true correlation matrix
#'
#' @param data          matrix of observations
#' @param corr_theo     true matrix of correlations
#' @param alpha        level of multiple testing (used if logical=TRUE)
#' @param stat_test     
#' \describe{
#'   \item{'empirical'}{\eqn{\sqrt{n}*abs(corr)}}
#'   \item{'fisher'}{\eqn{\sqrt{n-3}*1/2*\log( (1+corr)/(1-corr) )}}
#'   \item{'student'}{\eqn{\sqrt{n-2}*abs(corr)/\sqrt(1-corr^2)}}
#'   \item{'2nd.order'}{\eqn{\sqrt{n}*mean(Y)/sd(Y)} with \eqn{Y=(X_i-mean(X_i))(X_j-mean(X_j))}}
#' }
#' @param method        only 'MaxTinfty' implemented
#' @param Nboot         number of iterations for Monte-Carlo of bootstrap quantile evaluation
#' @param stepdown      logical, if TRUE a stepdown procedure is applied
#' @param vect         if TRUE returns a vector of adjusted p-values, corresponding to \code{vectorize(cor(data))};
#'                     if FALSE, returns an array containing the adjusted p-values for each entry of the correlation matrix 
#' @param logical      if TRUE, returns either a vector or a matrix where each element is equal to TRUE if the corresponding null hypothesis is rejected, and to FALSE if it is not rejected
#'                     if \code{stepdown=TRUE} and \code{logical=FALSE}, returns a list of successive p-values.
#' @param arr.ind      if TRUE, returns the indexes of the significant correlations, with repspect to level alpha
#'                      
#' @importFrom stats cor
#' @importFrom MASS mvrnorm
#' 
#' @references  Drton, M., & Perlman, M. D. (2007). Multiple testing and error control in Gaussian graphical model selection. Statistical Science, 22(3), 430-449.
#' @references  Roux, M. (2018). Graph inference by multiple testing with application to Neuroimaging, Ph.D., Université Grenoble Alpes, France, https://tel.archives-ouvertes.fr/tel-01971574v1.
#' @seealso ApplyFwerCor 
#' @seealso maxTinftyCor, maxTinftyCor_SD 
#'
#' @export
#' @examples
#' n <- 100
#' p <- 10
#' corr_theo <- diag(1,p)
#' corr_theo[1,3] <- 0.5
#' corr_theo[3,1] <- 0.5
#' data <- MASS::mvrnorm(n,rep(0,p),corr_theo)
#' # adjusted p-values:
#' (res <- ApplyFwerCor_oracle(data,corr_theo,stat_test='empirical',Nboot=1000,stepdown=FALSE))
#' # significant correlations, level alpha:
#' alpha <- 0.05
#' whichCor(res<alpha)
ApplyFwerCor_oracle <- function(data,corr_theo,alpha=c(),stat_test='empirical',method='MaxTinfty',Nboot=1000,stepdown=TRUE,vect=FALSE,logical=stepdown,arr.ind=FALSE){

 if(sum(stat_test==c('empirical','fisher','student','2nd.order'))==0){ stop('Wrong value for stat_test.')}

 # if( (!logical) & stepdown){
 #    warning('Stepdown procedures do not return p-values, logical is changed in TRUE.\n')
 #    logical = TRUE
 # }
 if(is.null(alpha) & (logical | arr.ind)){
    warning('For logical or indexes returns, a value of alpha is needed. alpha=0.05 is taken.\n')
    alpha = 0.05
 }
 if(arr.ind){ 
    vect <- FALSE
    logical <- TRUE 
 }
  
 Gtheo <- covDcorNorm(corr_theo,stat_test)
 if(method=='MaxTinfty'){
   	if(stepdown==FALSE){ res <- maxTinftyCor(data,alpha,stat_test,Nboot,Gtheo,vect,logical) }
    else{ res <- maxTinftyCor_SD(data,alpha,stat_test,Nboot,Gtheo,vect) }
 }
 else{ stop('The method is not implemented with oracle version.\n') }

 if(arr.ind){ res <- whichCor(res) }
 
}

#' Applies multiple testing procedures built to control (asymptotically) the FDR for correlation testing.
#' @description Applies multiple testing procedures built to control (asymptotically) the FDR for correlation testing.
#' Some have no theoretical proofs for tests on a correlation matrix.
#'
#'
#' @return  Returns either \itemize{\item{logicals indicating if the corresponding correlation is significant, as a vector or a matrix depending on \code{vect},} \item{an array containing indexes \eqn{\lbrace(i,j),\,i<j\rbrace} for which correlation between variables \eqn{i} and \eqn{j} is significant, if \code{arr.ind=TRUE}.}}
#
#' @param data      matrix of observations
#' @param alpha     level of multiple testing
#' @param stat_test  
#' \describe{
#'   \item{'empirical'}{\eqn{\sqrt{n}*abs(corr)}}
#'   \item{'fisher'}{\eqn{\sqrt{n-3}*1/2*\log( (1+corr)/(1-corr) )}}
#'   \item{'student'}{\eqn{\sqrt{n-2}*abs(corr)/\sqrt(1-corr^2)}}
#'   \item{'2nd.order'}{\eqn{\sqrt{n}*mean(Y)/sd(Y)} with \eqn{Y=(X_i-mean(X_i))(X_j-mean(X_j))}}
#' }
#' @param method    choice between 'LCTnorm' and 'LCTboot' developped by Cai & Liu (2016), 
#'                  'BH', traditional Benjamini-Hochberg's procedure Benjamini & Hochberg (1995)'s
#'                  and 'BHboot', Benjamini-Hochberg (1995)'s procedure with bootstrap evaluation of p-values
#' @param Nboot     number of iterations for bootstrap p-values evaluation
#' @param vect      if TRUE returns a vector of TRUE/FALSE values, corresponding to \code{vectorize(cor(data))};
#'                  if FALSE, returns an array containing rows and columns of significant correlations                  
#' @param arr.ind   if TRUE, returns the indexes of the significant correlations, with repspect to level alpha
#'
#' @importFrom stats cor
#' @importFrom MASS mvrnorm
#' @export
#'
#' @references Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the royal statistical society. Series B (Methodological), 289-300.
#' @references  Cai, T. T., & Liu, W. (2016). Large-scale multiple testing of correlations. Journal of the American Statistical Association, 111(513), 229-240.
#' @references  Roux, M. (2018). Graph inference by multiple testing with application to Neuroimaging, Ph.D., Université Grenoble Alpes, France, https://tel.archives-ouvertes.fr/tel-01971574v1.
#' @seealso ApplyFwerCor 
#' @seealso LCTnorm, LCTboot, BHCor, BHBootCor
#'
#' @examples
#' n <- 100
#' p <- 10
#' corr_theo <- diag(1,p)
#' corr_theo[1,3] <- 0.5
#' corr_theo[3,1] <- 0.5
#' data <- MASS::mvrnorm(n,rep(0,p),corr_theo)
#' res <- ApplyFdrCor(data,stat_test='empirical',method='LCTnorm')
#' # significant correlations, level alpha:
#' alpha <- 0.05
#' whichCor(res<alpha)
ApplyFdrCor <- function(data,alpha=0.05,stat_test='empirical',method='LCTnorm',Nboot=1000,vect=FALSE,arr.ind=FALSE){

  if(sum(stat_test==c('empirical','fisher','student','2nd.order'))==0){ stop('Wrong value for stat_test.')}
  if(sum(method==c('LCTnorm','LCTboot','BH','BHboot'))==0){ stop('Wrong value for method.')}
  if(arr.ind){ vect <- FALSE }
  
 if(method=='LCTnorm'){
        res <- LCTnorm(data,alpha=alpha,stat_test=stat_test,vect)
 } 
 if(method=='LCTboot'){
        res <- LCTboot(data,alpha=alpha,stat_test=stat_test,Nboot=Nboot,vect)
 } 
 if(method=='BH'){
       res <- BHCor(data,alpha=alpha,stat_test=stat_test,vect)
 }
 if(method=='BHboot'){
       res <- BHBootCor(data,alpha=alpha,stat_test=stat_test,Nboot=Nboot,vect)
 }

 if(arr.ind){ res <- whichCor(res) }
  
 return(res)


}


