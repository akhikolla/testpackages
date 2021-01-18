#' Cross-valiation for smog 
#' 
#' \code{cv.smog} conducts a greedy-search for optimal lambda's and yields a sparse
#' model given a provided model selection criterion. When type is ''nloglike'', 
#' the method allows the \code{nfolds} to be processed in parallel for speeding up 
#' the cross-validation.
#' 
#' @inheritParams smog.default
#' @param type model selction criterion, should choose from ''nloglike'', 
#'             ''cAIC'', ''AIC'', ''BIC'', and ''GCV'', respectively. 
#' @param lambda.max the maximum value for lambda's. If \code{NULL}, the default \code{lambda.max}
#'                   is \eqn{1/\lambda_{min}(x'x)}.
#' @param nlambda.max the maximum number of lambdas' shrunk down from the maximum lambda \code{lambda.max}.
#'                    Default is 20.
#' @param delta the damping rate for lambda's such that \eqn{\lambda_k = \delta^k\lambda_0}. Default is 0.9.
#' @param nfolds number of folds. One fold of the observations in the data are used
#'               as the testing, and the remaining are fitted for model training. 
#'               Default is 5.
#' @param parallel Whether or not process the \code{nfolds} cross-validations in
#'                 parallel. If \code{TRUE}, use \code{\link[foreach]{foreach}} to do each 
#'                 cross-validation in parallel. Default is \code{FALSE}.
#' @param ncores number of cpu's for parallel computing. See
#'               \code{\link[parallel]{makeCluster}} and \code{\link[doParallel]{registerDoParallel}}.
#'               Default is \code{NULL}. 
#' @param ... other arguments that can be supplied to \code{smog}.
#'
#' @return Includes the profile containing a path of lambda's and the corresponding model 
#'         selectio criterion value, the optimal lambda's, and the optimal model, respectively.
#'         The \code{type} comes from a list of model selection criteria values, includes the 
#'         average of the negative log-likelihood values and the correction AIC for each fold of the data.
#' 
#' \subsection{}{ 
#'  \describe{
#'    \item{cvfit}{the fitted model based on the optimal lambda's.}
#'    \item{lhat}{the optimal lambda's which has the minimum model selection criterion.}
#'    \item{profile}{a data frame contains the path of lambda's and the corresponding model selection
#'                   criterion, which is determined by the \code{type}.}
#'  }
#' }        
#' 
#' @details When the \code{type} is ''nloglike'', it requires doing \code{nfolds} cross-validations. 
#'          For each cross-validation, evenly split the whole data into \code{nfolds}, and one fold of 
#'          the observations are used as the testing data, and the remaining are used for model training. 
#'          After calculating the (deviance) residuals for each fold of testing data, return the average of the 
#'          (deviance) residuals. Note that we keep \code{lambda2}\eqn{=0} during the greedy search for lambda's. 
#' 
#' \subsection{Model selection criteria}{
#' Besides the n-fold cross-validation, `cv.smog` provides several AIC based model selection criteria.  
#' 
#' \describe{
#' \item{cAIC}{ \eqn{\frac{n}{2}}log(|2\*log-likelihood|) + \eqn{\frac{n}{2} (\frac{1+k/n}{1-k+2/n})} }       
#' \item{AIC}{ log(|2\*log-likelihood|/n) + \eqn{2\frac{k}{n}} }                         
#' \item{BIC}{ log(|2\*log-likelihood |/n) + \eqn{2\frac{k}{n}}log(n) }                         
#' \item{GCV}{ |2\*log-likelihood| / \eqn{(n(1-k/n)^2)} }
#' }
#' 
#' Where k is the degrees of freedom \code{DF}, which is related to the penalty parameters \eqn{\lambda}'s. 
#' }
#' 
#' @examples
#' 
#' # generate design matrix x
#' set.seed(2018)
#' n=100;p=20
#' s=10
#' x=matrix(0,n,1+2*p)
#' x[,1]=sample(c(0,1),n,replace = TRUE)
#' x[,seq(2,1+2*p,2)]=matrix(rnorm(n*p),n,p)
#' x[,seq(3,1+2*p,2)]=x[,seq(2,1+2*p,2)]*x[,1]
#' 
#' g=c(p+1,rep(1:p,rep(2,p)))  # groups 
#' v=c(0,rep(1,2*p))           # penalization status
#' label=c("t",rep(c("prog","pred"),p))  # type of predictor variables
#' 
#' # generate beta
#' beta=c(rnorm(13,0,2),rep(0,ncol(x)-13))
#' beta[c(2,4,7,9)]=0
#' 
#' # generate y
#' data=x%*%beta
#' noise=rnorm(n)
#' snr=as.numeric(sqrt(var(data)/(s*var(noise))))
#' y=data+snr*noise
#' 
#' cvfit=cv.smog(x,y,g,v,label,type = "GCV", family="gaussian")
#' 
#' @seealso \code{\link{smog.default}}, \code{\link{smog.formula}}, \code{\link{predict.smog}}, \code{\link{plot.smog}}.
#' @author Chong Ma, \email{chongma8903@@gmail.com}.
#' 
#' @references \insertRef{ma2019structural}{smog}
#' 
#' @importFrom foreach foreach %do% %dopar%
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' 
#' @export
cv.smog <- function(x, y, g, v, label, type = "nloglike", family = "gaussian", lambda.max = NULL, 
                    nlambda.max = 20, delta = 0.9, nfolds = 10, parallel = FALSE, ncores = NULL, ...){
  
  # define the maximum lambda value 
  typelist=c("nloglike","cAIC","AIC","BIC","GCV")
  n = nrow(x)
  p = ncol(x)
  
  if(!(type %in% typelist)){
    stop(paste(c(type,"is not one of",typelist,"!"),collapse = " "))
  }
  
  if(family == "gaussian"){
    wx <- cbind(rep(1,nrow(x)),x)
  }
  
  if(type %in% c("nloglike")){
    setlist = as.integer(seq(0,nrow(x),length.out = nfolds+1))
    ncores = ifelse(parallel,ncores,1)
    
    if(parallel){
      `%mypar%` =  `%dopar%`
      ncores = ifelse(is.null(ncores),1,ncores)
      cl <- parallel::makeCluster(ncores)
      doParallel::registerDoParallel(cl)
    }else{
      `%mypar%` = `%do%`
    }
  }
  
  if(family == "gaussian"){
    lambda.max = ifelse(is.null(lambda.max),10*max(t(x)%*%y)/n*max(log(max(p,n)/n),1),lambda.max)
  }
  
  if(family == "binomial"){
    lambda.max = ifelse(is.null(lambda.max),1/min(svd(x)$d)*max(log(max(p,n)/n),1),lambda.max) 
  }
  
  if(family == "coxph"){
    lambda.max = ifelse(is.null(lambda.max),1/min(svd(x)$d)*p/n,lambda.max) 
  }
  
  temp.fit0 = smog(x,y,g,v,label,lambda1=lambda.max,lambda2=0,lambda3=lambda.max,family,...)
  
  if(nrow(coef(temp.fit0))<4){
    repeat{
      lambda.max = lambda.max*delta
      temp.fit1 = smog(x,y,g,v,label,lambda1=lambda.max,lambda2=0,lambda3=lambda.max,family,...)
      if(nrow(coef(temp.fit1))>=4)  break
    }
  }else{
    repeat{
      lambda.max = lambda.max/delta
      temp.fit1 = smog(x,y,g,v,label,lambda1=lambda.max,lambda2=0,lambda3=lambda.max,family,...)
      if(nrow(coef(temp.fit1))<=4)  break
    }
  }
  
  l1 = l2 = lambda.max
  profile = list(NULL)
  lhat = cvfit = NULL
  
  it = 1    # initialize the search step
  temp=rep(NA,3)
  repeat{
    temp.lambda=list(c(l1*delta,0,l2),c(l1,0,l2*delta),c(l1*delta,0,l2*delta))
    
    if(!(type %in% c("nloglike"))){
      # non-cross-validation type model selection
      for(k in 1:3){
        temp[k] = smog(x,y,g,v,label,lambda1=temp.lambda[[k]][1],lambda2=temp.lambda[[k]][2],lambda3=temp.lambda[[k]][3],family,...)$criteria[[type]]
      }
    }else{
      # cross-validation to calculate the maximum log-likelihood 
      for(k in 1:3){
        i = 1
        cv_res <- foreach::foreach(i=1:nfolds,.combine = c,
                                   .packages = c("foreach"))%mypar%{
                                     tset <- (setlist[i]+1):setlist[i+1]
                                     lset <- c(1:nrow(x))[-tset]
                                     cvmodel <- smog::smog(x,y,g,v,label,lambda1=temp.lambda[[k]][1],lambda2=temp.lambda[[k]][2],lambda3=temp.lambda[[k]][3],family,subset=lset,...)
                                     nt <- length(tset)
                                     kt = cvmodel$DF
                                     
                                     if(family == "gaussian"){
                                       tx = wx[tset,coef(cvmodel)$Id+1]
                                       ty = as.vector(as.matrix(tx) %*% as.matrix(coef(cvmodel)$Estimate))
                                       tloglike = -1/2*sum((y[tset] - ty)^2)
                                     }
                                     
                                     if(family == "binomial"){
                                       probTAB = exp(as.matrix(x[tset,coef(cvmodel)$Id])%*%as.matrix(coef(cvmodel)$Estimate))/
                                         (1+exp(as.matrix(x[tset,coef(cvmodel)$Id])%*%as.matrix(coef(cvmodel)$Estimate)))
                                       probTAB = cbind(1-rowSums(as.matrix(probTAB)),probTAB)
                                       ty = match(y[tset],cvmodel$levels)
                                       tloglike = sum(log(apply(cbind(ty,probTAB),1,function(t) t[-1][t[1]]))) 
                                     }
                                     
                                     if(family == "coxph"){
                                       tx=x[tset,coef(cvmodel)$Id]
                                       ty=y[tset,]
                                       
                                       ttheta = exp(as.matrix(tx)%*%as.matrix(coef(cvmodel)$Estimate))
                                       tloglike = as.numeric(colSums(as.matrix(tx))%*%as.matrix(coef(cvmodel)$Estimate))
                                       
                                       t=NULL
                                       for(t in ty$time[ty$status>0]){
                                         tid0 = which(ty$time == t & ty$status>0)
                                         tid1 = which(ty$time >= t)
                                         ti = NULL
                                         for(ti in 0:(length(tid0)-1)){
                                           tloglike = tloglike - log(sum(ttheta[tid1])-ti/length(tid0)*sum(ttheta[tid0]))
                                         }
                                       }
                                     }
                                     
                                     c(-tloglike)
                                   }
        temp[k] = mean(cv_res, na.rm = TRUE)
      }
    }
    
    if(1 %in% which(temp == min(temp)) | 3 %in% which(temp == min(temp)) 
       | min(temp[c(1,3)]) - temp[2] < 0.001*abs(min(temp[c(1,3)]))
    ){
      profile[[it]] = c(temp.lambda[[3]][c(1,3)],temp[3])
    }else{
      profile[[it]] = c(temp.lambda[[2]][c(1,3)],temp[2])
    }
    
    l1 = profile[[it]][1]
    l2 = profile[[it]][2]
    
    # determine the stopping criterion
    if(it >= 5){
      if((abs(profile[[it-1]][3] - profile[[it]][3]) > 0.01*abs(profile[[it-1]][3]) & 
          abs(profile[[it-2]][3] - profile[[it-1]][3]) > 0.01*abs(profile[[it-2]][3]) &
          (profile[[it-2]][3] - profile[[it-1]][3] > 10*(profile[[it-1]][3] - profile[[it]][3]) |
           (profile[[it-2]][3] - profile[[it-1]][3] > 0.01*abs(profile[[it-2]][3]) & 
            profile[[it-1]][3] - profile[[it]][3] > 10*(profile[[it-2]][3] - profile[[it-1]][3])))
      ) | it >= nlambda.max){
        break
      } 
    }
    it = it + 1
  }
  
  if(parallel)  parallel::stopCluster(cl)
  
  profile = data.frame(do.call(rbind,profile))
  colnames(profile) = c("lambda1","lambda2",type)
  lhat = as.numeric(profile[which.min(profile[,3]),1:2])
  cvfit = smog(x,y,g,v,label,lambda1=lhat[1],lambda2=0,lambda3=lhat[2],family,...)$model
  
  cv_result = list(cvfit = cvfit, lhat = lhat, profile = profile)
  cv_result$call = match.call()
  class(cv_result) = "cv.smog"
  
  cv_result
}



#' plot method for objects of `cv.smog` class
#' 
#' Yield a search path for optimal group penalty \eqn{G-\lambda} and \eqn{I-\lambda} using
#' the mean-squared errors from the cross-validations. 
#' 
#' @param x An fitted object in "cv.smog" class.
#' @param ... Other graphical parameters to ggplot2.
#' 
#' @details 
#' x-axis represents the group tuning parameter \eqn{\lambda_G} and y-axis for the interaction 
#' tuning parameter \eqn{\lambda_I}, respectively. The point size demonstrates the maganitude
#' of MSE or negative log-likelihood. 
#' 
#' @seealso [smog], [cv.smog], [cv.cglasso].
#' 
#' @author Chong Ma, \email{chongma8903@@gmail.com}.
#' @references \insertRef{ma2019structural}{smog}
#' 
#' @import ggplot2
#' @examples 
#' 
#' # generate design matrix x
#' set.seed(2018)
#' n=100;p=20
#' s=10
#' x=matrix(0,n,1+2*p)
#' x[,1]=sample(c(0,1),n,replace = TRUE)
#' x[,seq(2,1+2*p,2)]=matrix(rnorm(n*p),n,p)
#' x[,seq(3,1+2*p,2)]=x[,seq(2,1+2*p,2)]*x[,1]
#' 
#' g=c(p+1,rep(1:p,rep(2,p)))  # groups 
#' v=c(0,rep(1,2*p))           # penalization status
#' label=c("t",rep(c("prog","pred"),p))  # type of predictor variables
#' 
#' # generate beta
#' beta=c(rnorm(13,0,2),rep(0,ncol(x)-13))
#' beta[c(2,4,7,9)]=0
#' 
#' # generate y
#' data=x%*%beta
#' noise=rnorm(n)
#' snr=as.numeric(sqrt(var(data)/(s*var(noise))))
#' y=data+snr*noise
#' 
#' cvfit=cv.smog(x,y,g,v,label,type = "AIC", family="gaussian")
#' plot(cvfit)
#' 
#' @export
plot.cv.smog <- function(x, ...){
  
  profile = x$profile
  profile_name = colnames(profile)
  
  p = ggplot(data = profile, aes(x = log(lambda1), y = log(lambda2))) +
    geom_line(color = alpha("grey50", 0.8)) +
    geom_point(aes(size = get(profile_name[3])), color = alpha("red",0.8)) +
    labs(x = expression(log(lambda[G])), y = expression(log(lambda[I]))) +
    geom_vline(xintercept = log(x$lhat[1]), linetype = "dashed", color = alpha("grey30",0.8)) +
    geom_hline(yintercept = log(x$lhat[2]), linetype = "dashed", color = alpha("grey30",0.8)) +
    theme(panel.grid.minor = element_blank(),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)) +
    guides(size = guide_legend(profile_name[3]))
  
  suppressWarnings(print(p))
}







