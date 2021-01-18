#' cross-validation for conditional group lasso method
#' 
#' This function uses the cross-validation approach to search the optimal group tuning 
#' paramete \eqn{\lambda_1}, conditional on fixing \eqn{\lambda_2} and \eqn{\lambda_3}
#' at a small value. 
#' 
#' @inheritParams cv.smog
#' @param ratio The ratio of \eqn{\lambda_1} and \eqn{\lambda_2} to `lambda.max`. Smaller
#'              value means less penalty on the coefficients of interactions.  
#' 
#' @importFrom stats sd
#' @details 
#' The idea of this conditional group lasso function is to reduce the computing time,
#' by merely searching the optimal group penalty rather than searching a grid of two-dimensional
#' penalties. By controling the ridge and interaction penalties at a small value, it 
#' still honors the hierarchy structure, but also leverage the multicolliearity problems.   
#' 
#' @examples
#' 
#' # generate design matrix x
#' set.seed(2018)
#' n=50;p=20
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
#' cvfit=cv.cglasso(x,y,g,v,label,family="gaussian",nlambda.max = 20)
#' 
#' @seealso \code{\link{smog.default}}, \code{\link{smog.formula}}, [cv.smog]
#' @author Chong Ma, \email{chongma8903@@gmail.com}.
#' 
#' @references \insertRef{ma2019structural}{smog}
#' 
#' @export
cv.cglasso <- function(x, y, g, v, label, family = "gaussian", lambda.max = NULL, 
                       nlambda.max = 50, delta = 0.9, nfolds = 10, ratio = 0.01, parallel = FALSE, 
                       ncores = NULL, ...) {
  n = nrow(x)
  p = ncol(x)
  cl = NULL
  
  if (family == "gaussian") {
    wx <- cbind(rep(1, nrow(x)), x)
    lambda.max = ifelse(is.null(lambda.max), 10*max(t(x)%*%y)/n*max(log(max(p, n)/n), 1), lambda.max)
  }
  
  if (family == "binomial") {
    lambda.max = ifelse(is.null(lambda.max), 1/min(svd(x)$d)*max(log(max(p, n)/n), 1), lambda.max)
  }
  if (family == "coxph") {
    lambda.max = ifelse(is.null(lambda.max), 1/min(svd(x)$d)*p/n, lambda.max)
  }
  
  temp.fit0 = smog(x, y, g, v, label, lambda1 = lambda.max, lambda2 = lambda.max*ratio, 
                   lambda3 = lambda.max*ratio, family, ...)
  
  if (nrow(coef(temp.fit0)) < 4) {
    repeat {
      lambda.max = lambda.max*delta
      temp.fit1 = smog(x, y, g, v, label, lambda1 = lambda.max, lambda2 = lambda.max*ratio, 
                       lambda3 = lambda.max*ratio, family, ...)
      if (nrow(coef(temp.fit1)) >= 4) 
        break
    }
  } else {
    repeat {
      lambda.max = lambda.max/delta
      temp.fit1 = smog(x, y, g, v, label, lambda1 = lambda.max, lambda2 = lambda.max*ratio, 
                       lambda3 = lambda.max*ratio, family, ...)
      if (nrow(coef(temp.fit1)) <= 4) 
        break
    }
  }
  
  setlist = as.integer(seq(0, nrow(x), length.out = nfolds + 1))
  lambda = exp(seq(log(lambda.max), log(lambda.max*delta^nlambda.max), length.out = nlambda.max))
  cv_res = matrix(NA, nlambda.max, nfolds)
  
  if (parallel) {
    `%mypar%` = `%dopar%`
    ncores = ifelse(is.null(ncores), 1, ncores)
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
  } else {
    `%mypar%` = `%do%`
  }
  
  cv_fit = cv_mean = cv_sd = wh_min = lambda.min = lambda.1se = NULL
  
  j = 1
  for(j in 1:nlambda.max){
    i = 1
    cv_res[j,] <- foreach::foreach(i = 1:nfolds, .combine = c, .packages = c("foreach")) %mypar% {
      tset <- (setlist[i] + 1):setlist[i + 1]
      lset <- c(1:nrow(x))[-tset]
      cvmodel <- smog::smog(x, y, g, v, label, lambda1 = lambda[j], lambda2 = lambda[j]*ratio, 
                      lambda3 = lambda[j]*ratio, family, subset = lset, ...)
      nt <- length(tset)
      kt = cvmodel$DF
      
      if (family == "gaussian") {
        tx = wx[tset, coef(cvmodel)$Id + 1]
        ty = as.vector(as.matrix(tx) %*% as.matrix(coef(cvmodel)$Estimate))
        tloglike = -1/2 * mean((y[tset] - ty)^2)
      }
      
      if (family == "binomial") {
        probTAB = exp(as.matrix(x[tset, coef(cvmodel)$Id]) %*% as.matrix(coef(cvmodel)$Estimate))/
          (1+exp(as.matrix(x[tset, coef(cvmodel)$Id]) %*% as.matrix(coef(cvmodel)$Estimate)))
        probTAB = cbind(1 - rowSums(as.matrix(probTAB)), probTAB)
        ty = match(y[tset], cvmodel$levels)
        tloglike = mean(log(apply(cbind(ty, probTAB), 1, function(t) t[-1][t[1]])))
      }
      
      if (family == "coxph") {
        tx = x[tset, coef(cvmodel)$Id]
        ty = y[tset, ]
        ttheta = exp(as.matrix(tx) %*% as.matrix(coef(cvmodel)$Estimate))
        tloglike = as.numeric(colSums(as.matrix(tx)) %*% as.matrix(coef(cvmodel)$Estimate))
        t = NULL
        for (t in ty$time[ty$status > 0]) {
          tid0 = which(ty$time == t & ty$status > 0)
          tid1 = which(ty$time >= t)
          ti = NULL
          for (ti in 0:(length(tid0) - 1)) {
            tloglike = tloglike - log(sum(ttheta[tid1]) - ti/length(tid0) * sum(ttheta[tid0]))
          }
        }
      }
      c(-tloglike)
    }
  }
  
  if (parallel) 
    parallel::stopCluster(cl)
  
  cv_mean = apply(cv_res, 1, mean)
  cv_sd = apply(cv_res, 1, sd)
  wh_min = which.min(cv_mean)
  lambda.min = lambda[wh_min]
  lambda.1se = tryCatch(min(lambda[cv_mean >= cv_mean[wh_min] + cv_sd[wh_min]]),
                        warning = function(w) NA,
                        error = function(e) NA)
  lambda.opt = max(lambda.min, lambda.1se, na.rm = TRUE)
  
  profile = data.frame(cv_mean, cv_sd, lambda)
  cvfit = smog(x, y, g, v, label, lambda1 = lambda.opt, lambda2 = lambda.opt*ratio, 
               lambda3 = lambda.opt*ratio, family, ...)$model
  
  cv_obj = list(cvfit = cvfit, profile = profile, lambda.min = lambda.min, lambda.1se = lambda.1se)
  cv_obj$call = match.call()
  class(cv_obj) = "cv.cglasso"
  
  cv_obj 
}


#' plot method for objects of `cv.cglasso` class
#' 
#' Yields a cross-validation curve, and error bars within one standard deviation of the curve, 
#' as a function of the group penalty \eqn{\lambda_1}.
#' 
#' @param x An fitted object in "cv.cglasso" class.
#' @param ... Other graphical parameters to ggplot2.
#' 
#' @seealso [cv.cglasso], [cv.smog], [smog].
#' 
#' @author Chong Ma, \email{chongma8903@@gmail.com}.
#' @references \insertRef{ma2019structural}{smog}
#' @examples 
#' 
#' # generate design matrix x
#' set.seed(2018)
#' n=50;p=20
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
#' cvfit=cv.cglasso(x,y,g,v,label,family="gaussian", nlambda.max = 20)
#' plot(cvfit)
#' 
#' @export
plot.cv.cglasso <- function(x, ...){
  
  p = ggplot(data = x$profile, aes(x = log(lambda), y = cv_mean), ...) +
    geom_line() +
    geom_point(color = alpha("red",0.8), size = 3) +
    geom_errorbar(aes(ymin = cv_mean - cv_sd, ymax = cv_mean + cv_sd), width = 0.05,
                  color = alpha("grey50",0.8)) +
    labs(x = expression(log(lambda)), y = "MSE") +
    geom_vline(xintercept = log(x$lambda.min), linetype = "dashed", color = alpha("grey30",0.8)) +
    geom_vline(xintercept = log(x$lambda.1se), linetype = "dashed", color = alpha("grey30",0.8)) +
    theme(panel.grid.minor = element_blank(),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12))
  
  suppressWarnings(print(p))
}



