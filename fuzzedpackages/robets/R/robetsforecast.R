#' Forecasting using ROBETS models
#'
#' Returns forecasts and other information for univariate ROBETS models.
#' 
#' @param object An object of class "\code{robets}". Usually the result of a call to \code{\link{robets}}.
#' @param h Number of periods for forecasting
#' @param level Confidence level for prediction intervals.
#' @param PI If \code{TRUE}, prediction intervals are calculated.
#' @param lambda Box-Cox transformation parameter. Ignored if NULL. Otherwise, forecasts back-transformed via an inverse Box-Cox transformation.
#' @param ... Other arguments.
#' @return An object of class "\code{forecast}". The function \code{summary} is used to obtain and print a summary of the results, while the function \code{plot} produces a plot of the forecasts. The generic accessor functions \code{fitted.values} and \code{residuals} extract useful features of the value returned by \code{forecast.robets}. An object of class \code{"forecast"} is a list containing at least the following elements:
#' \itemize{
#' \item{model: }{A list containing information about the fitted model}
#' \item{method: }{The name of the forecasting method as a character string}
#' \item{mean: }{Point forecasts as a time series}
#' \item{x: }{The original time series (either \code{object} itself or the time series used to create the model stored as \code{object}).}
#' \item{residuals: }{Residuals from the fitted model. For models with additive errors, the residuals are x - fitted values. For models with multiplicative errors, the residuals are equal to x /(fitted values) - 1.}
#' \item{fitted: }{Fitted values (one-step ahead forecasts)}
#' }
#' @details The code of this function is based on the function \code{forecast.ets} of the package \code{forecast} of Hyndman and Khandakar (2008).
#' @examples 
#' library(forecast)
#' model <- robets(nottem)
#' plot(forecast(model))
#' @references Crevits, R., and Croux, C (2016) "Forecasting with Robust Exponential Smoothing with Damped Trend and Seasonal Components".\emph{Working paper}. \url{https://doi.org/10.13140/RG.2.2.11791.18080}
#' @references Hyndman, R. J., and Khandakar, Y (2008) "Automatic time series forecasting: The forecasting package for R".\emph{Journal of Statistical Software} \strong{27}(3). \url{https://doi.org/10.18637/jss.v027.i03}
#' @author Ruben Crevits, \email{ruben.crevits@@kuleuven.be}, \url{https://rcrevits.wordpress.com/research}
#' @seealso \code{\link{robets}}
#'
#' @export
forecast.robets <- function(object, h=ifelse(object$m>1, 2*object$m, 10),
  level=c(80,95), PI=TRUE, lambda=object$lambda, ...)
{
  # Check inputs
  #if(h>2000 | h<=0)
  if(h <= 0)
    stop("Forecast horizon out of bounds")
  if(!PI)
  {
    level <- 90
  }

  n <- length(object$x)
  damped <- as.logical(object$components[4])

  if(object$components[1]=="A" & is.element(object$components[2],c("A","N")) & is.element(object$components[3],c("N","A")))
    f <- class1(h,object$states[n+1,(2:ncol(object$states))],object$components[2],object$components[3],damped,object$m,object$states[n+1,1]^2,object$par)
  else if(object$components[1]=="M" & is.element(object$components[2],c("A","N")) & is.element(object$components[3],c("N","A")))
    f <- class2(h,object$states[n+1,(2:ncol(object$states))],object$components[2],object$components[3],damped,object$m,object$states[n+1,1]^2,object$par)
  else if(object$components[1]=="M" & object$components[3]=="M" & object$components[2]!="M")
    f <- class3(h,object$states[n+1,(2:ncol(object$states))],object$components[2],object$components[3],damped,object$m,object$states[n+1,1]^2,object$par)
  else
    stop("Cannot forecast this model")
  if(!PI)
    f$var <- f$lower <- f$upper <- NULL

  tsp.x <- tsp(object$x)
  if(!is.null(tsp.x))
    start.f <- tsp(object$x)[2] + 1/object$m
  else
    start.f <- length(object$x)+1
  out <- list(model=object,mean=ts(f$mu,frequency=object$m,start=start.f),level=level,x=object$x)
  if(!is.null(f$var))
  {
    out$lower <- out$upper <- ts(matrix(NA,ncol=length(level),nrow=h))
    for(i in 1:length(level))
    {
      marg.error <- sqrt(f$var) * abs(qnorm((100-level[i])/200))
      out$lower[,i] <- out$mean - marg.error
      out$upper[,i] <- out$mean + marg.error
    }
    tsp(out$lower) <- tsp(out$upper) <- tsp(out$mean)
  }
  else if(!is.null(f$lower))
  {
    out$lower <- ts(f$lower)
    out$upper <- ts(f$upper)
    tsp(out$lower) <- tsp(out$upper) <- tsp(out$mean)
  }
  else if(PI)
    warning("No prediction intervals for this model")
		
  out$fitted <- fitted(object)
  out$method <- object$method
  out$residuals <- residuals(object)

  if(!is.null(lambda))
  {
	  #out$x <- InvBoxCox(object$x,lambda)
	  #out$fitted <- InvBoxCox(out$fitted,lambda)
    out$mean <- InvBoxCox(out$mean,lambda)
	  if(PI)  # PI = TRUE
	  {
		out$lower <- InvBoxCox(out$lower,lambda)
		out$upper <- InvBoxCox(out$upper,lambda)
	  }
  }
	
  return(structure(out,class="forecast"))
}

class1 <- function(h,last.state,trendtype,seasontype,damped,m,sigma2,par)
{
  p <- length(last.state)
  H <- matrix(c(1,rep(0,p-1)),nrow=1)
  if(seasontype=="A")
    H[1,p] <- 1
  if(trendtype=="A")
  {
    if(damped)
      H[1,2] <- par["phi"]
    else
      H[1,2] <- 1
  }
  F <- matrix(0,p,p)
  F[1,1] <- 1
  if(trendtype=="A")
  {
    if(damped)
      F[1,2] <- F[2,2] <- par["phi"]
    else
      F[1,2] <- F[2,2] <- 1
  }
  if(seasontype=="A")
  {
    F[p-m+1,p] <- 1
    F[(p-m+2):p,(p-m+1):(p-1)] <- diag(m-1)
  }
  G <- matrix(0,nrow=p,ncol=1)
  G[1,1] <- par["alpha"]
  if(trendtype=="A")
    G[2,1] <- par["beta"]
  if(seasontype=="A")
    G[3,1] <- par["gamma"]
  mu <- numeric(h)
  Fj <- diag(p)
  cj <- numeric(h-1)
  if(h>1)
  {
    for(i in 1:(h-1))
    {
      mu[i] <- H %*% Fj %*% last.state
      Fj <- Fj %*% F
      cj[i] <- H %*% Fj %*% G
    }
    cj2 <- cumsum(cj^2)
    var <- sigma2 * c(1,1+cj2)
  }
  else
    var <- sigma2
  mu[h] <- H %*% Fj %*% last.state

  return(list(mu=mu,var=var,cj=cj))
}

class2 <- function(h,last.state,trendtype,seasontype,damped,m,sigma2,par)
{
  tmp <- class1(h,last.state,trendtype,seasontype,damped,m,sigma2,par)
  theta <- numeric(h)
  theta[1] <- tmp$mu[1]^2
  if(h>1)
  {
    for(j in 2:h)
      theta[j] <- tmp$mu[j]^2 + sigma2 * sum(tmp$cj[1:(j-1)]^2*theta[(j-1):1])
  }
  var <- (1+sigma2)*theta - tmp$mu^2
  return(list(mu=tmp$mu,var=var))
}

class3 <- function(h,last.state,trendtype,seasontype,damped,m,sigma2,par)
{
  p <- length(last.state)
  H1 <- matrix(rep(1,1+(trendtype!="N")),nrow=1)
  H2 <- matrix(c(rep(0,m-1),1),nrow=1)
  if(trendtype=="N")
  {
    F1 <- 1
    G1 <- par["alpha"]
  }
  else
  {
    F1 <- rbind(c(1,1),c(0,ifelse(damped,par["phi"],1)))
    G1 <- rbind(c(par["alpha"],par["alpha"]),c(par["beta"],par["beta"]))
  }
  F2 <- rbind(c(rep(0,m-1),1),cbind(diag(m-1),rep(0,m-1)))

  G2 <- matrix(0,m,m)
  G2[1,m] <- par["gamma"]
  Mh <- matrix(last.state[1:(p-m)]) %*% matrix(last.state[(p-m+1):p],nrow=1)
  Vh <- matrix(0,length(Mh),length(Mh))
  H21 <- H2 %x% H1
  F21 <- F2 %x% F1
  G21 <- G2 %x% G1
  K <- (G2 %x% F1) + (F2 %x% G1)
  mu <- var <- numeric(h)
  for(i in 1:h)
  {
    mu[i] <- H1 %*% Mh %*% t(H2)
    var[i] <- (1+sigma2) * H21 %*% Vh %*% t(H21) + sigma2*mu[i]^2
    vecMh <- c(Mh)
    Vh <- F21 %*% Vh %*% t(F21) + sigma2 * (F21 %*% Vh %*% t(G21) + G21 %*% Vh %*% t(F21) +
      K %*% (Vh + vecMh %*% t(vecMh)) %*% t(K) + sigma2 * G21 %*% (3*Vh + 2*vecMh%*%t(vecMh))%*%t(G21))
    Mh <- F1 %*% Mh %*% t(F2) + G1 %*% Mh %*% t(G2) * sigma2
  }
  return(list(mu=mu,var=var))
}

