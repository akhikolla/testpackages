

ypbpSurv <- function(time, z, par, tau, degree,
                     baseline=c("hazard", "odds")){
  baseline <- match.arg(baseline)
  q <- length(z)
  base <- bp(time, degree, tau)

  psi <- par[1:q]
  phi <- par[(q+1):(2*q)]
  gamma <- par[(2*q+1):(2*q+degree)]
  theta_S <- exp( as.numeric(z%*%psi) )
  theta_L <- exp( as.numeric(z%*%phi) )


  if(baseline=="hazard"){
    Ht0 <- as.numeric(base$B%*%gamma) + as.numeric(time>tau)*(time-tau)*degree*gamma[degree]/tau
    St0 <- exp(-Ht0)
    Ft0 <- 1-St0
    St <- exp( -theta_L*(log(theta_L*St0 + theta_S*Ft0)-(as.numeric(z%*%phi) - Ht0)  ) )
  }else{
    ratio <- exp( as.numeric(z%*%psi) - as.numeric(z%*%phi) )
    Rt0 <- as.numeric(base$B%*%gamma) + as.numeric(time>tau)*(time-tau)*degree*gamma[degree]/tau
    St <- exp(-theta_L*log(1 + ratio*Rt0))
  }
  return(St)
}


ypbpSurv2 <- function(time, z, x, par, tau, degree,
                      baseline=c("hazard", "odds")){
  baseline <- match.arg(baseline)
  q <- length(z)
  p <- length(x)
  base <- bp(time, degree, tau)

  psi <- par[1:q]
  phi <- par[(q+1):(2*q)]
  beta <- par[(2*q+1):(2*q+p)]
  gamma <- par[(2*q+p+1):(2*q+p+degree)]
  theta_S <- exp( as.numeric(z%*%psi) )
  theta_L <- exp( as.numeric(z%*%phi) )
  theta_C <- exp( as.numeric(x%*%beta) )
  ratio <- exp( as.numeric(z%*%psi) - as.numeric(z%*%phi) )


  if(baseline=="hazard"){
    Ht0 <- as.numeric(base$B%*%gamma) + as.numeric(time>tau)*(time-tau)*degree*gamma[degree]/tau
    St0 <- exp(-Ht0)
    Ft0 <- 1-St0
    St <- exp( -theta_L*theta_C*(log(theta_L*St0 + theta_S*Ft0)-(as.numeric(z%*%phi) - Ht0)  ) )
  }else{
    St <- exp(-theta_L*theta_C*log(1 + ratio*Rt0))
    Rt0 <- as.numeric(base$B%*%gamma) + as.numeric(time>tau)*(time-tau)*degree*gamma[degree]/tau
  }
  return(St)
}


#---------------------------------------------
#' survfit method for ypbp models
#'
#' @aliases survfit.ypbp
#' @description Computes the predicted survivor function for a ypbp model.
#' @importFrom survival survfit
#' @export
#' @param formula an object of the class ypbp
#' @param newdata a data frame containing the set of explanatory variables.
#' @param ... further arguments passed to or from other methods.
#' @return  a list containing the estimated survival probabilities.
#' @examples
#' \donttest{
#' # ML approach:
#' library(YPBP)
#' mle <- ypbp(Surv(time, status)~arm, data=ipass, approach="mle")
#' summary(mle)
#' ekm <- survival::survfit(Surv(time, status)~arm, data=ipass)
#' newdata <- data.frame(arm=0:1)
#' St <- survfit(mle, newdata)
#' plot(ekm, col=1:2)
#' with(St, lines(time, surv[[1]]))
#' with(St, lines(time, surv[[2]], col=2))
#'
#' # Bayesian approach:
#' bayes <- ypbp(Surv(time, status) ~ arm, data = ipass,
#'               approach = "bayes", chains = 2, iter = 100)
#' summary(bayes)
#' ekm <- survival::survfit(Surv(time, status)~arm, data=ipass)
#' newdata <- data.frame(arm=0:1)
#' St <- survfit(bayes, newdata)
#' plot(ekm, col=1:2)
#' with(St, lines(time, surv[[1]]))
#' with(St, lines(time, surv[[2]], col=2))
#' }
#'
survfit.ypbp <- function(formula, newdata, ...){
  object <- formula
  mf <- object$mf
  labels <- names(mf)[-1]
  baseline <- object$baseline
  time <- sort( stats::model.response(mf)[,1])
  status <- sort( stats::model.response(mf)[,2])
  data <- data.frame(cbind(time, status, mf[,-1]))
  names(data) <- c("time", "status", names(mf)[-1])
  degree <- object$degree
  tau <- object$tau
  labels <- match.arg(names(newdata), labels, several.ok = TRUE)
  formula <- object$formula
  Z <- as.matrix(stats::model.matrix(formula, data = newdata, rhs = 1)[,-1])
  X <- suppressWarnings(try( as.matrix(stats::model.matrix(formula, data = newdata, rhs = 2)[,-1]), TRUE))
  St <- list()


  if(object$approach=="mle"){
    par <- object$fit$par
    if(object$p==0){
      for(i in 1:nrow(newdata)){
        St[[i]] <- ypbpSurv(time, Z[i, ], par, tau, degree, baseline)
      }
    }else{
      for(i in 1:nrow(newdata)){
        St[[i]] <- ypbpSurv2(time, Z[i, ], X[i, ], par, tau, degree, baseline)
      }
    }
  }else{ # Bayesian approach
    samp <- rstan::extract(object$fit)
    if(object$p==0){
      par <- cbind(samp$psi, samp$phi, samp$gamma)
      for(i in 1:nrow(newdata)){
        aux <- apply(par, 1, ypbpSurv, time=time, z=Z[i,], tau=tau, degree=degree, baseline=baseline)
        St[[i]] <- apply(aux, 1, mean)
      }
    }else{
      par <- cbind(samp$psi, samp$phi, samp$beta, samp$gamma)
      for(i in 1:nrow(newdata)){
        aux <- apply(par, 1, ypbpSurv, time=time, z=Z[i,], x=X[i, ], tau=tau, degree=degree, baseline=baseline)
        St[[i]] <- apply(aux, 1, mean)
      }
    }
  }

  out <- list(time = time, surv = St)
  class(out) <- "survfit.ypbp"
  return(out)
}


ypbpCrossSurv <- function(z1, z2, par, tau0, tau, degree, baseline){
  diff_St <- function(time, z1, z2, par, tau, degree, baseline){
    St1 <- ypbpSurv(time=time, z=z1, par=par, tau=tau, degree=degree, baseline=baseline)
    St2 <- ypbpSurv(time=time, z=z2, par=par, tau=tau, degree=degree, baseline=baseline)
    return(St1-St2)
  }
  I <- c(tau0, 1.5*tau)
  t <- try(stats::uniroot(diff_St, interval=I, z1=z1, z2=z2, par=par, tau=tau, degree=degree, baseline=baseline)$root, TRUE)

  if(class(t)=="try-error")
  {
    return(NA)
  }else{
    return(t)
  }
}

ypbpCrossSurv2 <- function(z1, z2, x, par, tau0, tau, degree, baseline){
  diff_St <- function(time, z1, z2, x, par, tau, degree, baseline){
    St1 <- ypbpSurv2(time=time, z=z1, x=x, par=par, tau=tau, degree=degree, baseline=baseline)
    St2 <- ypbpSurv2(time=time, z=z2, x=x, par=par, tau=tau, degree=degree, baseline=baseline)
    return(St1-St2)
  }
  I <- c(tau0, 1.5*tau)
  #I <- c(tau0, tau)
  t <- try(stats::uniroot(diff_St, interval=I, z1=z1, z2=z2, x=x, par=par, tau=tau, degree=degree, baseline=baseline)$root, TRUE)
  if(class(t)=="try-error")
  {
    return(NA)
  }else{
    return(t)
  }
}



#---------------------------------------------
#' Generic S3 method crossTime
#' @aliases crossTime
#' @export
#' @param object a fitted model object
#' @param ... further arguments passed to or from other methods.
#' @return the crossing survival time
#'
crossTime <- function(object, ...) UseMethod("crossTime")


#' Computes the crossing survival times
#'
#' @aliases crossTime.ypbp
#' @description Computes the crossing survival times along with their corresponding confidence/credible intervals.
#' @rdname crossTime-methods
#' @method crossTime ypbp
#' @export
#' @export crossTime
#' @param object an object of class ypbp
#' @param newdata1 a data frame containing the first set of explanatory variables
#' @param newdata2 a data frame containing the second set of explanatory variables
#' @param conf.level level of the confidence/credible intervals; default is conf.level = 0.95
#' @param nboot number of bootstrap samples (default nboot=4000); ignored if approach="bayes".
#' @param ... further arguments passed to or from other methods.
#' @return  the crossing survival time
#' @examples
#' \donttest{
#' # ML approach:
#' library(YPBP)
#' mle <- ypbp(Surv(time, status)~arm, data=ipass, approach="mle")
#' summary(mle)
#' newdata1 <- data.frame(arm=0)
#' newdata2 <- data.frame(arm=1)
#' tcross <- crossTime(mle, newdata1, newdata2, nboot = 100)
#' tcross
#' ekm <- survival::survfit(Surv(time, status)~arm, data=ipass)
#' newdata <- data.frame(arm=0:1)
#' St <- survfit(mle, newdata)
#' plot(ekm, col=1:2)
#' with(St, lines(time, surv[[1]]))
#' with(St, lines(time, surv[[2]], col=2))
#' abline(v=tcross, col="blue")
#'
#' # Bayesian approach:
#' bayes<-ypbp(Surv(time,status)~arm,data=ipass,approach="bayes",chains=2,iter=100)
#' summary(bayes)
#' newdata1 <- data.frame(arm=0)
#' newdata2 <- data.frame(arm=1)
#' tcross <- crossTime(bayes, newdata1, newdata2)
#' tcross
#' ekm <- survival::survfit(Surv(time, status)~arm, data=ipass)
#' newdata <- data.frame(arm=0:1)
#' St <- survfit(bayes, newdata)
#' plot(ekm, col=1:2)
#' with(St, lines(time, surv[[1]]))
#' with(St, lines(time, surv[[2]], col=2))
#' abline(v=tcross, col="blue")
#' }
#'
crossTime.ypbp <- function(object, newdata1, newdata2,
                           conf.level=0.95, nboot=4000, ...){
  q <-object$q
  p <-object$p
  mf <- object$mf
  labels <- names(mf)[-1]
  baseline <- object$baseline
  time <- stats::model.response(mf)[,1]
  status <- stats::model.response(mf)[,2]
  o <- order(time)
  data <- data.frame(cbind(time, status, mf[,-1]))[o,]
  names(data) <- c("time", "status", labels)
  tau0 <- min(time)
  degree <- object$degree
  tau <- object$tau
  labels <- match.arg(names(newdata1), names(newdata2), several.ok=TRUE)
  labels <- match.arg(names(mf)[-1], names(newdata1), several.ok=TRUE)
  z1 <- matrix(stats::model.matrix(object$formula, data = newdata1, rhs = 1)[,-1], ncol=q)
  z2 <- matrix(stats::model.matrix(object$formula, data = newdata2, rhs = 1)[,-1], ncol=q)
  if(p>0){
    x <- matrix(stats::model.matrix(object$formula, data = newdata2, rhs = 2)[,-1], ncol=p)
  }

  I <- c(tau0, 1.5*tau)
  alpha <- 1 - conf.level
  prob <- c(alpha/2, 1-alpha/2)

  if(object$approach=="mle"){
    t <- c()
    par <- object$fit$par
    for(i in 1:nrow(newdata1)){
      if(p==0){
        t[i] <- ypbpCrossSurv(z1=z1[i,], z2=z2[i,], par=par, tau0, tau, degree=degree, baseline=baseline)
      }else{
        t[i] <- ypbpCrossSurv2(z1=z1[i,], z2=z2[i,], x=x[i,], par=par, tau0, tau, degree=degree, baseline=baseline)
      }
    }
    par <- with(object, ypbpBoot(formula=formula, data=data, n_int=n_int,
                                rho=rho, tau=tau, nboot=nboot, prob=prob))

    ci <- matrix(nrow=nrow(newdata1), ncol=2)
    if(object$p==0){
      for(i in 1:nrow(newdata1)){
        aux <- apply(par, 1, ypbpCrossSurv, z1=z1[i,], z2=z2[i,], tau0, tau, degree=degree, baseline=baseline)
        ci[i,] <- stats::quantile(aux, probs=prob, na.rm=TRUE)
      }
    }else{
      for(i in 1:nrow(newdata1)){
        aux <- apply(par, 1, ypbpCrossSurv2, z1=z1[i,], z2=z2[i,], x=x[i,], tau0, tau, degree=degree, baseline=baseline)
        ci[i,] <- stats::quantile(aux, probs=prob, na.rm=TRUE)
      }
    }

    t <- data.frame(cbind(t, ci))
    names(t) <- c("Est.", paste(100*prob, "%", sep=""))
  }else{ # Bayesian approach
    t <- matrix(nrow=nrow(newdata1), ncol=3)
    samp <- rstan::extract(object$fit)
    if(object$p==0){
      par <- cbind(samp$psi, samp$phi, samp$gamma)
      for(i in 1:nrow(newdata1)){
        aux <- apply(par, 1, ypbpCrossSurv, z1=z1[i,], z2=z2[i,], tau0, tau, degree=degree, baseline=baseline)
        ci <- stats::quantile(aux, probs=prob, na.rm=TRUE)
        t[i,] <- c(mean(aux, na.rm=TRUE), ci)
      }
    }else{
      par <- cbind(samp$psi, samp$phi, samp$beta, samp$gamma)
      for(i in 1:nrow(newdata1)){
        aux <- apply(par, 1, ypbpCrossSurv2, z1=z1[i,], z2=z2[i,], x=x[i,], tau0, tau, degree=degree, baseline=baseline)
        ci <- stats::quantile(aux, probs=prob, na.rm=TRUE)
        t[i,] <- c(mean(aux, na.rm=TRUE), ci)
      }
    }
    t <- as.data.frame(t)
    names(t) <- c("Est.", names(ci) )
  }
  #class(t) <- "crossTime.ypbp"
  return(t)
}

