#' Evaluates the test statistics for tests on correlation matrix entries.
#'
#' @return Returns the test statistics for correlation testing.
#'
#' @param data   matrix of observations
#' @param type     
#' \describe{
#'   \item{'empirical'}{\eqn{\sqrt{n}*abs(corr)}}
#'   \item{'fisher'}{\eqn{\sqrt{n-3}*1/2*\log( (1+corr)/(1-corr) )}}
#'   \item{'student'}{\eqn{\sqrt{n-2}*abs(corr)/\sqrt(1-corr^2)}}
#'   \item{'2nd.order'}{\eqn{\sqrt{n}*mean(Y)/sd(Y)} with \eqn{Y=(X_i-mean(X_i))(X_j-mean(X_j))}}
#' }
#'
#' @export
#'
#' @examples
#' n <- 100
#' p <- 10
#' corr_theo <- diag(1,p)
#' data <- MASS::mvrnorm(n,rep(0,p),corr_theo)
#' stat <- eval_stat(data,'fisher')
eval_stat <- function(data,type='empirical')
{
  n <- nrow(data)
  p <- ncol(data)
  corr_mat <- cor(data)
  corr_vect <- vectorize(corr_mat)

    stat <- sqrt(n)*(corr_vect)

    if(type=='student'){
        if(n>2){
           stat <- sqrt(n-2)*(corr_vect)/sqrt(1-corr_vect^2) 
        }else{
           stat <- sqrt(n)*(corr_vect)/sqrt(1-corr_vect^2)  } 
    }

	if(type=='fisher'){
	    if(n>2){ 
           stat <- sqrt(n-2)*1/2*log( (1+corr_vect)/(1-corr_vect) ) 
        }else{
           stat <- sqrt(n)*1/2*log( (1+corr_vect)/(1-corr_vect) ) }
   }


 	if(type=='2nd.order'){
        meanX <- colMeans(data)
	    x <- (data -meanX)
        theta <- array(0,dim(corr_mat))
        for(ncol in 2:p){
          for(nrow in 1:(ncol-1)){
	        theta[nrow,ncol] <- var(x[,nrow]*x[,ncol])
          }
        } 
        stat <- vectorize(sqrt(n/theta)*corr_mat)
    }


  return(stat)

}





#' Returns the theoretical covariance of test statistics for correlation testing.
#'
#' @param cor_mat    A correlation matrix
#' @param stat_test     
#' \describe{
#'   \item{'empirical'}{\eqn{\sqrt{n}*abs(corr)}}
#'   \item{'fisher'}{\eqn{\sqrt{n-3}*1/2*\log( (1+corr)/(1-corr) )}}
#'   \item{'student'}{\eqn{\sqrt{n-2}*abs(corr)/\sqrt(1-corr^2)}}
#'   \item{'2nd.order'}{\eqn{\sqrt{n}*mean(Y)/sd(Y)} with \eqn{Y=(X_i-mean(X_i))(X_j-mean(X_j))}}
#' }
#' 
#' @return Returns the theoretical covariance of the test statistics.
#'
#' @seealso covDcor, covD2nd, eval_stat
#'
#' @importFrom stats cor
#' @export
#'
#' @examples
#' p <- 10
#' corr_theo <- diag(1,p)
#' corr_theo[2:p,] <- 0.3
#' corr_theo[,2:p] <- 0.3
#' covDcorNorm(corr_theo,stat_test='student')
covDcorNorm <- function(cor_mat,stat_test='empirical')
{
 
    w <- covDcor(cor_mat)

    if(stat_test=='student'){
        cor_vect <- vectorize(cor_mat)
	      fact <- diag(1/sqrt(1-cor_vect^2)^3)
	      w <- fact %*% w %*% fact 
    }
    if(stat_test=='fisher'){
        cor_vect <- vectorize(cor_mat)
	      fact <- diag(1/(1-cor_vect^2))
	      w <- fact %*% w %*% fact 
    }
    if(stat_test=='2nd.order'){
        w <- covD2nd(cor_mat)  
    }


  return(w)

}


#' Returns the indexes of an upper triangular matrix with logical entries.
#'
#' @param mat    A matrix with logical entries in the upper triangular part
#' @return Returns the indexes of the upper triangular part where the entries are TRUE
#' @export
#'
#' @examples
#' n <- 100
#' p <- 10
#' corr_theo <- diag(1,p)
#' corr_theo[1,3] <- 0.5
#' corr_theo[3,1] <- 0.5
#' data <- MASS::mvrnorm(n,rep(0,p),corr_theo)
#' res <- ApplyFwerCor(data,stat_test='empirical',method='Bonferroni',stepdown=FALSE)
#' # significant correlations, level alpha:
#' alpha <- 0.05
#' whichCor(res<alpha)
whichCor <- function(mat){
  ind <- NULL
  if(length(mat)>0){
    ind <- which(mat,arr.ind=TRUE)
    ind <- ind[(ind[,'row']<ind[,'col']),]
  }
  return(ind)
}


