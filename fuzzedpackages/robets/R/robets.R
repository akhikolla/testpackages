
#' Robust exponential smoothing model
#'
#' Returns robets model applied to \code{y}.
#' 
#' @param y a numeric vector or time series
#' @param model A three-letter string indicating the method using the framework terminology of Hyndman et al. (2008). The first letter denotes the error type ("A", "M" or "Z"); the second letter denotes the trend type ("N","A" or "Z"); and the third letter denotes the season type ("N","A","M" or "Z"). In all cases, "N"=none, "A"=additive, "M"=multiplicative and "Z"=automatically selected. So, for example, "ANN" is simple exponential  smoothing with additive errors, "MAM" is multiplicative Holt-Winters' method with multiplicative errors, and so on. It is also possible for the model to be of class "\code{robets}", and equal to the output from a previous call to \code{robets}. In this case, the same model is fitted to \code{y} without re-estimating any smoothing parameters. See also the \code{use.initial.values} argument.
#' @param damped If TRUE, use a damped trend. If NULL, both damped and non-damped trends will be tried and the best model (according to the information criterion \code{ic}) will be returned.
#' @param alpha Value of alpha. If NULL, it is estimated.
#' @param beta Value of beta. If NULL, it is estimated.
#' @param gamma Value of gamma. If NULL, it is estimated.
#' @param phi Value of phi. If NULL, it is estimated.
#' @param additive.only If TRUE, will only consider additive models. Default is FALSE.
#' @param lambda Box-Cox transformation parameter. Ignored if NULL. Otherwise, data transformed before model is estimated. When \code{lambda=TRUE}, \code{additive.only} is set to FALSE.
#' @param lower Lower bounds for the parameters (alpha, beta, gamma, phi)
#' @param upper Upper bounds for the parameters (alpha, beta, gamma, phi)
#' @param opt.crit Optimization criterion. One of "roblik" (Robust Log-likelihood, default), "\link{tau2}" (Tau squared error of the residuals), "mse" (Mean Square Error), "amse" (Average MSE over first \code{nmse} forecast horizons), "sigma" (Standard deviation of residuals), "mae" (Mean of absolute residuals), or "lik" (Log-likelihood).
#' @param bounds Type of parameter space to impose: \code{"usual"} indicates all parameters must lie between specified lower and upper bounds; \code{"admissible"} indicates parameters must lie in the admissible space; \code{"both"} (default) takes the intersection of these regions.
#' @param ic Information criterion to be used in model selection.
#' @param use.initial.values If \code{TRUE} (default) and \code{model} is of class \code{"robets"}, then the initial values in the model are also not re-estimated.
#' @param opt.initial.values If \code{FALSE} (default) a robust heuristic is used for chosing the initial values. If \code{TRUE}  the initial values are part of the problem to optimize \code{opt.crit}.  Neglected if \code{use.initial.values} is \code{TRUE} and \code{model} is of class \code{"robets"}.
#' @param rob.start.initial.values If \code{TRUE} (default) the initial values are computed via the robust heuristic described in Crevits and Croux (2016). If \code{FALSE} the initial values are computed via the same heuristic as in Hyndman et al. (2008). The initial values computed with these methods are further optimized if \code{opt.initial.values} is \code{TRUE}.
#' @param opt.sigma0 If \code{FALSE} (default) sigma0 is equal to the value computed together with the other initial values via a heuristic. If \code{TRUE} sigma0 is included as a variable in the optimization problem. It is not recommended to set \code{opt.sigma0 = TRUE}.
#' @param k Value of k in forecasting equations. \code{k=3} is default. If NULL, \code{k} is included as a variable in the optimization problem. It is not recommended to set \code{k = NULL}.
#' @param nmse Number of steps for AMSE (1<=\code{nmse}<=30), \code{nmse=1} is default.
#' @param ... Other undocumented arguments.
#' @return An object of class "\code{robets}".
#' @details The code is an extended version of the code of the function \code{ets} of the package \code{forecast} of Hyndman and Khandakar (2008). The methodology is an extended version of Gelper et al. (2008). In Crevits and Croux (2016) the methodology of \code{robets} is described in full.
#' @examples
#' library(forecast)
#' model <- robets(nottem)
#' plot(forecast(model))
#' @references Crevits, R., and Croux, C (2016) "Forecasting with Robust Exponential Smoothing with Damped Trend and Seasonal Components".\emph{Working paper}. \url{https://doi.org/10.13140/RG.2.2.11791.18080}
#' @references Gelper S., Fried R. and Croux C. (2010) "Robust Forecasting with Exponential and Holt-Winters Smoothing".\emph{Journal of Forecasting}, \strong{29}, 285-300. \url{https://doi.org/10.1002/for.1125}
#' @references Hyndman, R. J., and Khandakar, Y (2008) "Automatic time series forecasting: The forecasting package for R".\emph{Journal of Statistical Software} \strong{27}(3). \url{https://doi.org/10.18637/jss.v027.i03}
#' @author Ruben Crevits, \email{ruben.crevits@@kuleuven.be}, \url{https://rcrevits.wordpress.com/research}
#' @seealso \code{\link{forecast.robets}, \link{plot.robets}, \link{plotOutliers}, \link{tau2}, \link{ets}}
#' @export
robets <- function(y, model="ZZZ", damped=NULL,
                alpha=NULL, beta=NULL, gamma=NULL, phi=NULL, additive.only=FALSE, lambda=NULL,
                lower=c(rep(0.0001,3), 0.8), upper=c(rep(0.9999,3),0.98),
                opt.crit=c("roblik","tau2","lik","mse","amse","sigma","mae"), bounds=c("both","usual","admissible"),
                ic=c("robaicc","robaic","robbic","aicc","bic","aic"), 
                use.initial.values=TRUE,opt.initial.values=FALSE,rob.start.initial.values=TRUE,opt.sigma0=FALSE,k=3,nmse=1,...)
{
  #dataname <- substitute(y)
  opt.crit <- match.arg(opt.crit)
  bounds <- match.arg(bounds)
  ic <- match.arg(ic)
  
  if(any(class(y) %in% c("data.frame","list","matrix","mts")))
    stop("y should be a univariate time series")
  y <- as.ts(y)
  
  # Remove missing values near ends
  ny <- length(y)
  y <- na.contiguous(y)
  if(ny != length(y))
    warning("Missing values encountered. Using longest contiguous portion of time series")
  
  orig.y <- y
  if(class(model)=="robets" & is.null(lambda))
    lambda <- model$lambda
  if(!is.null(lambda))
  {
    y <- BoxCox(y,lambda)
    additive.only=TRUE
  }
  
  m <- frequency(y)
  
  if(sum((upper-lower)>0)<4)
    stop("Lower limits must be less than upper limits")
  
  # If model is an robets object, re-fit model to new data
  if(class(model)=="robets")
  {
    if(use.initial.values == FALSE && opt.initial.values == FALSE)
      stop("Set either use.initial.values=TRUE or opt.initial.values=TRUE.")
    alpha <- model$par["alpha"]
    beta <- model$par["beta"]
    if(is.na(beta))
      beta <- NULL
    gamma <- model$par["gamma"]
    if(is.na(gamma))
      gamma <- NULL
    phi <- model$par["phi"]
    if(is.na(phi))
      phi <- NULL
    k <- model$par["k"]
    if(is.na(k))
      k <- 3
    modelcomponents <- paste(model$components[1],model$components[2],model$components[3],sep="")
    damped <- (model$components[4]=="TRUE")
    if(use.initial.values)
    {
      errortype  <- substr(modelcomponents,1,1)
      trendtype  <- substr(modelcomponents,2,2)
      seasontype <- substr(modelcomponents,3,3)
      
      # Recompute errors from pegelsresid.C
      e <- robpegelsresid.C(y, m, model$initstate, errortype, trendtype, seasontype, damped, alpha, beta, gamma, phi,nmse,k)
      
      # Compute error measures
      np <- length(model$par)
      model$loglik <- -0.5*e$lik
      model$aic <- e$lik + 2*np
      model$bic <- e$lik + log(ny)*np
      model$aicc <- e$lik + 2*ny*np/(ny-np-1)
      model$roblik <- e$roblik
      model$robaic <- e$roblik + 2*np
      model$robbic <- e$roblik + log(ny)*np
      model$robaicc <- e$roblik + 2*ny*np/(ny-np-1)
      model$mse <- e$amse[1]
      model$amse <- mean(e$amse)
      model$tau2 <- e$tau2
      
      # Compute states, fitted values and residuals
      tsp.y <- tsp(y)
      model$states=ts(e$states,frequency=tsp.y[3],start=tsp.y[1]-1/tsp.y[3])
      colnames(model$states)[1] <- "sigma"
      colnames(model$states)[2] <- "l"
      if(trendtype!="N")
        colnames(model$states)[3] <- "b"
      if(seasontype!="N")
        colnames(model$states)[(3+(trendtype!="N")):ncol(model$states)] <- paste("s",1:m,sep="")
      if(errortype=="A")
        model$fitted <- ts(y-e$e,frequency=tsp.y[3],start=tsp.y[1])
      else
        model$fitted <- ts(y/(1+e$e),frequency=tsp.y[3],start=tsp.y[1])
      model$residuals <- ts(e$e,frequency=tsp.y[3],start=tsp.y[1])
      model$sigma2 <- mean(model$residuals^2,na.rm=TRUE)
      model$x <- orig.y
      if(!is.null(lambda))
      {
        model$fitted <- InvBoxCox(model$fitted,lambda)
      }
      model$lambda <- lambda
      
      # Return model object
      return(model)
    }
    else
    {
      model <- modelcomponents
    }
  }
  
  errortype  <- substr(model,1,1)
  trendtype  <- substr(model,2,2)
  seasontype <- substr(model,3,3)
  
  if(!is.element(errortype,c("M","A","Z")))
    stop("Invalid error type")
  if(is.element(trendtype,c("M"))) #NO MULTIPLICATIVE TREND
    stop("Invalid trend type: no multiplicative trend allowed")
  if(!is.element(trendtype,c("N","A","Z"))) 
    stop("Invalid trend type")
  if(!is.element(seasontype,c("N","A","M","Z")))
    stop("Invalid season type")
  
  if(m < 1 | length(y) <= m)
  {
    #warning("I can't handle data with frequency less than 1. Seasonality will be ignored.")
    seasontype <- "N"
  }
  if(m == 1)
  {
    if(seasontype=="A" | seasontype=="M")
      stop("Nonseasonal data")
    else
      substr(model,3,3) <- seasontype <- "N"
  }
  if(m > 24)
  {
    if(is.element(seasontype,c("A","M")))
      stop("Frequency too high")
    else if(seasontype=="Z")
    {
      warning("I can't handle data with frequency greater than 24. Seasonality will be ignored.")
      substr(model,3,3) <- seasontype <- "N"
      #m <- 1
    }
  }
  
  data.positive <- (min(y) > 0)
  
  if(!data.positive & errortype=="M")
    stop("Inappropriate model for data with negative or zero values")
  
  if(!is.null(damped))
  {
    if(damped & trendtype=="N")
      stop("Forbidden model combination")
  }
  
  # Check we have enough data to fit a model
  n <- length(y)
  npars <- 2L # alpha + l0
  if(trendtype=="A") 
    npars <- npars + 2L # beta + b0
  if(seasontype=="A" | seasontype=="M")
    npars <- npars + m # gamma + s
  if(!is.null(damped))
    npars <- npars + as.numeric(damped)
  if(n <= npars + 1)
    stop("You've got to be joking. I need more data!")
  
  # Fit model (assuming only one nonseasonal model)
  if(errortype=="Z")
    errortype <- c("A","M")
  if(trendtype=="Z")
  {
    trendtype <- c("N","A") # NO MULTIPLICATIVE TREND
  }
  if(seasontype=="Z")
    seasontype <- c("N","A","M")
  if(is.null(damped))
    damped <- c(TRUE,FALSE)
  best.ic <- Inf
  for(i1 in 1:length(errortype))
  {
    for(i2 in 1:length(trendtype))
    {
      for(i3 in 1:length(seasontype))
      {
        for(i4 in 1:length(damped))
        {
          if(trendtype[i2]=="N" & damped[i4])
            next
          if(!data.positive & errortype[i1]=="M")
            next
          if(errortype[i1]=="A" & seasontype[i3]=="M")
            next
          if(additive.only & (errortype[i1]=="M" | seasontype[i3]=="M"))
            next
          fit <- robetsmodel(y,errortype[i1],trendtype[i2],seasontype[i3],damped[i4],alpha,beta,gamma,phi,lower=lower,upper=upper,opt.crit=opt.crit,bounds=bounds,opt.initial.values=opt.initial.values,rob.start.initial.values=rob.start.initial.values,opt.sigma0=opt.sigma0,k=k,nmse=nmse,...)
          fit.ic <- switch(ic,aic=fit$aic,bic=fit$bic,aicc=fit$aicc,robaic=fit$robaic,robbic=fit$robbic,robaicc=fit$robaicc)
          if(!is.na(fit.ic))
          {
            if(fit.ic < best.ic)
            {
              model <- fit
              best.ic <- fit.ic
              best.e <- errortype[i1]
              best.t <- trendtype[i2]
              best.s <- seasontype[i3]
              best.d <- damped[i4]
            }
          }
        }
      }
    }
  }
  if(best.ic == Inf)
    stop("No model able to be fitted")
  
  model$m <- m
  model$method <- paste("ROBETS(",best.e,",",best.t,ifelse(best.d,"d",""),",",best.s,")",sep="")
  model$components <- c(best.e,best.t,best.s,best.d)
  model$call <- match.call()
  model$initstate <- model$states[1,]
  model$sigma2 <- mean(model$residuals^2,na.rm=TRUE)
  model$x <- orig.y
  model$lambda <- lambda
  if(!is.null(lambda))
  {
    model$fitted <- InvBoxCox(model$fitted,lambda)
  }
  
  if(is.null(k)) k <- model$par["k"]
  sigmas <- model$states[2:nrow(model$states),1]
  model$outlyingness <- model$residuals/sigmas 
  model$outliers <-  abs(model$outlyingness) > k
  model$outliers <- ts(model$outlier, frequency = tsp(y)[3], start = tsp(y)[1])
  #model$call$data <- dataname
  
  return(structure(model,class="robets"))
}

robpegelsresid.C <- function(y,m,initstate,errortype,trendtype,seasontype,damped,alpha,beta,gamma,phi,nmse,k)
{
  if(!damped){
    phi <- 1;
  }
  if(trendtype == "N"){
    beta <- 0;
  }
  if(seasontype == "N"){
    gamma <- 0;
  }
  
  res <- .Call("robets_calc_out",y,m,initstate,switch(errortype,"A"=1,"M"=2), switch(trendtype,"N"=0,"A"=1,"M"=2), switch(seasontype,"N"=0,"A"=1,"M"=2),damped,alpha,beta,gamma,phi,nmse,k,PACKAGE ="robets")
  if(!is.na(res$lik)){
    if(abs(res$lik+99999) < 1e-7)
      res$lik <- NA
  }
  
  tsp.y <- tsp(y)
  e <- ts(res$e)
  tsp(e) <- tsp.y  
  
  return(list(lik=res$lik,amse=res$amse[1:nmse],roblik=res$roblik,tau2=res$tau2, e=res$e, states=matrix(res$states, nrow=length(y)+1, ncol=(2+(trendtype!="N")+m*(seasontype!="N")), byrow=TRUE)))
}

robetsmodel <- function(y, errortype, trendtype, seasontype, damped,
                     alpha=NULL, beta=NULL, gamma=NULL, phi=NULL,
                     lower, upper, opt.crit, bounds, maxit=2000, control=NULL, seed=NULL, trace=FALSE, 
                     opt.initial.values = TRUE,rob.start.initial.values=TRUE,opt.sigma0=FALSE, k=NULL,nmse)
{
  tsp.y <- tsp(y)
  if(is.null(tsp.y))
    tsp.y <- c(1,length(y),1)
  
  if(seasontype != "N"){
    m <- tsp.y[3]
  }else{
    m <- 1
  }
  
  # Initialize smoothing parameters
  par <- initparam(alpha,beta,gamma,phi,trendtype,seasontype,damped,lower,upper,m)
  names(alpha) <- names(beta) <- names(gamma) <- names(phi) <- NULL
  par.noopt <- c(alpha=alpha,beta=beta,gamma=gamma,phi=phi)
  if(!is.null(par.noopt))
    par.noopt <- c(na.omit(par.noopt))
  if(!is.na(par["alpha"]))
    alpha <- par["alpha"]
  if(!is.na(par["beta"]))
    beta <- par["beta"]
  if(!is.na(par["gamma"]))
    gamma <- par["gamma"]
  if(!is.na(par["phi"]))
    phi <- par["phi"]
  
  #    if(errortype=="M" | trendtype=="M" | seasontype=="M")
  #        bounds="usual"
  if(!check.param(alpha,beta,gamma,phi,lower,upper,bounds,m))
  {
    print(paste("Model: ROBETS(",errortype,",",trendtype,ifelse(damped,"d",""),",",seasontype,")",sep=""))
    stop("Parameters out of range")
  }
  
  # Initialize state
  if(rob.start.initial.values){
    initval <- robinitstate(y,errortype,trendtype,seasontype) 
  }else{
    initval <- initstate2(y,errortype,trendtype,seasontype) 
  }
  
  sigma0 <- initval[1]
  names(sigma0) <- "sigma0"
  initstate <- initval[2:length(initval)]
  if(opt.sigma0){
    par <- c(par,sigma0)
    lower <- c(lower,sigma0 = 0)
    upper <- c(upper,sigma0 = Inf)   
  }else{
    par.noopt <- c(par.noopt,sigma0)
  }
  if(opt.initial.values){
    par <- c(par,initstate = initstate)
    bo <- initstatebounds(y,errortype,trendtype,seasontype)
    lower <- c(lower,initstate = bo$lower)
    upper <- c(upper,initstate = bo$upper)    
  }else{
    par.noopt <- c(par.noopt,initstate = initstate)
  }
  if(is.null(k)){
    par <- c(par,k=100)
    lower <- c(lower,k=1)
    upper <- c(upper,k=200)  
    optK <- TRUE
  }else if(!(class(k)=="numeric")){
    par.noopt <- c(par.noopt,k=3)
    optK <- FALSE
  }else if(k<=0){
    par.noopt <- c(par.noopt,k=3)
    optK <- FALSE
  }else{
    par.noopt <- c(par.noopt,k=unname(k))
    optK <- FALSE
  }
  
  np <- length(par)
  if(np >= length(y)-1) # Not enough data to continue
    return(list(roblik=Inf,robaic=Inf,robbic=Inf,robaicc=Inf,
                loglik=-Inf,aic=Inf,bic=Inf,aicc=Inf,mse=Inf,
                amse=Inf,tau2=Inf,fit=NULL,states=initstate,par=par))  

  env <- robetsTargetFunctionInit(par=par, y=y, errortype=errortype, trendtype=trendtype,
                             seasontype=seasontype, damped=damped, par.noopt=par.noopt, lowerb=lower, upperb=upper,
                             opt.crit=opt.crit, nmse=as.integer(nmse), bounds=bounds, m=m,pnames=names(par),pnames2=names(par.noopt),
                             opt.sigma0=opt.sigma0, opt.initial.values = opt.initial.values, optK = optK)

  fred <- .Call("robetsNelderMead", par, env, -Inf, sqrt(.Machine$double.eps), 1.0, 0.5, 2.0, trace, maxit, PACKAGE = "robets")
  
    fit.par <- fred$par
    names(fit.par) <- names(par)
  if(opt.sigma0 && !is.na(fit.par["sigma0"])){
    sigma0 <- fit.par["sigma0"]
  }
  if(opt.initial.values){
    initstate <- fit.par[grep('^initstate', names(fit.par))]
  }
  if(optK && !is.na(fit.par["k"])){
    k <- fit.par["k"]
  }
  
  if(!is.na(fit.par["alpha"]))
    alpha <- fit.par["alpha"]
  if(!is.na(fit.par["beta"]))
    beta <- fit.par["beta"]
  if(!is.na(fit.par["gamma"]))
    gamma <- fit.par["gamma"]
  if(!is.na(fit.par["phi"]))
    phi <- fit.par["phi"]
  e <- robpegelsresid.C(y,m,c(sigma0,initstate),errortype,trendtype,seasontype,damped,alpha,beta,gamma,phi,nmse,k)
  
  n <- length(y)  
  loglik <- -0.5*e$lik
  aic <- e$lik + 2*np
  bic <- e$lik + log(n)*np
  aicc <- e$lik + 2*n*np/(n-np-1)
  roblik <- e$roblik
  robaic <- e$roblik + 2*np
  robbic <- e$roblik + log(n)*np
  robaicc <- e$roblik + 2*n*np/(n-np-1)
  mse <- e$amse[1]
  amse <- mean(e$amse)
  tau2 <- e$tau2  
  
  states=ts(e$states,frequency=tsp.y[3],start=tsp.y[1]-1/tsp.y[3])
  colnames(states)[1] <- "sigma"
  colnames(states)[2] <- "l"
  if(trendtype!="N")
    colnames(states)[3] <- "b"
  if(seasontype!="N")
    colnames(states)[(3+(trendtype!="N")):ncol(states)] <- paste("s",1:m,sep="")
  
  
  fit.par <- c(fit.par,par.noopt)
  #    fit.par <- fit.par[order(names(fit.par))]
  if(errortype=="A"){
    fits <- y-e$e
  }else{
    fits <- y/(1+e$e)
  }
  
  return(list(roblik=e$roblik,robaic=robaic,robbic=robbic,robaicc=robaicc,loglik=-0.5*e$lik,aic=aic,bic=bic,aicc=aicc,mse=mse,amse=amse,tau2=tau2,fit=fred,residuals=ts(e$e,frequency=tsp.y[3],start=tsp.y[1]),fitted=ts(fits,frequency=tsp.y[3],start=tsp.y[1]),states=states,par=fit.par))
}

robetsTargetFunctionInit <- function(par,y,errortype,trendtype,seasontype,damped,par.noopt,lowerb,upperb,
                                  opt.crit,nmse,bounds,m,pnames,pnames2,opt.sigma0,opt.initial.values,optK)
{  
  names(par) <- pnames
  names(par.noopt) <- pnames2
  alpha <- c(par["alpha"],par.noopt["alpha"])["alpha"]
  if(is.na(alpha))
    stop("alpha problem!")
  if(trendtype!="N"){
    beta <- c(par["beta"],par.noopt["beta"])["beta"]
    if(is.na(beta))
      stop("beta Problem!")
  }
  else
    beta <- NULL
  if(seasontype!="N"){
    gamma <- c(par["gamma"],par.noopt["gamma"])["gamma"]
    if(is.na(gamma))
      stop("gamma Problem!")
  }
  else{
    m <- 1
    gamma <- NULL
  }
  if(damped){
    phi <- c(par["phi"],par.noopt["phi"])["phi"]
    if(is.na(phi))
      stop("phi Problem!")
  }else
    phi <- NULL
  
  #determine which values to optimize and which ones are given by the user/not needed
  optAlpha <- !is.null(alpha)
  optBeta <- !is.null(beta)
  optGamma <- !is.null(gamma)
  optPhi <- !is.null(phi)
  
  givenAlpha <- FALSE
  givenBeta <- FALSE
  givenGamma <- FALSE
  givenPhi <- FALSE
  
  if(!is.null(par.noopt["alpha"])) if(!is.na(par.noopt["alpha"])) {
    optAlpha <- FALSE
    givenAlpha <- TRUE
  }
  if(!is.null(par.noopt["beta"])) if(!is.na(par.noopt["beta"])) {
    optBeta <- FALSE
    givenBeta <- TRUE
  }
  if(!is.null(par.noopt["gamma"])) if(!is.na(par.noopt["gamma"])) {
    optGamma <- FALSE
    givenGamma <- TRUE
  }
  if(!is.null(par.noopt["phi"])) if(!is.na(par.noopt["phi"])) {
    optPhi <- FALSE
    givenPhi <- TRUE
  }
  
  givenInit <- !opt.initial.values
  givenK <- !optK
  givenSigma0 <- !opt.sigma0
  
  initstate <- c(par[grep('^initstate', names(par))],par.noopt[grep('^initstate', names(par.noopt))])
  if(is.na(sum(initstate)))
    stop("initstate Problem!")  
  
  sigma0 <- c(par["sigma0"],par.noopt["sigma0"])["sigma0"]
  if(is.na(sigma0))
    stop("sigma0 Problem!")
  
  k <- c(par["k"],par.noopt["k"])["k"]
  if(is.na(k))
    stop("k Problem!")
  
  
  if(!damped)
    phi <- 1;
  if(trendtype == "N")
    beta <- 0;
  if(seasontype == "N")
    gamma <- 0;
  
  env <- new.env()
  
  res <- .Call("robetsTargetFunctionInit", y=y, errortype=switch(errortype,"A"=1,"M"=2),
               trendtype=switch(trendtype,"N"=0,"A"=1,"M"=2), seasontype=switch(seasontype,"N"=0,"A"=1,"M"=2),
               damped=damped, lowerb=lowerb, upperb=upperb,
               opt.crit=opt.crit, nmse=as.integer(nmse), bounds=bounds, m=m,
               optAlpha, optBeta, optGamma, optPhi, opt.sigma0, opt.initial.values, optK,
               givenAlpha, givenBeta, givenGamma, givenPhi, givenSigma0, givenInit, givenK,
               alpha, beta, gamma, phi,sigma0,initstate,k, env, PACKAGE = "robets")
  res
}

initparam <- function(alpha,beta,gamma,phi,trendtype,seasontype,damped,lower,upper,m)
{
  # Set up initial parameters
  par <- numeric(0)
  if(is.null(alpha))
  {
    if(m > 12)
      alpha <- 0.0002
    if(is.null(beta) & is.null(gamma))
      alpha <- lower[1] + .5*(upper[1]-lower[1])
    else if(is.null(gamma))
      alpha <- beta+0.001
    else if(is.null(beta))
      alpha <- 0.999-gamma
    else
      alpha <- 0.5*(beta - gamma + 1)
    if(alpha < lower[1] | alpha > upper[1])
      stop("Inconsistent parameter limits")
    par <- alpha
    names(par) <- "alpha"
  }
  if(is.null(beta))
  {
    if(trendtype !="N")
    {
      if(m > 12)
        beta <- 0.00015
      else
        beta <- lower[2] + .1*(upper[2]-lower[2])
      if(beta > alpha)
        beta <- min(alpha - 0.0001,0.0001)
      if(beta < lower[2] | beta > upper[2])
        stop("Can't find consistent starting parameters")
      par <- c(par,beta)
      names(par)[length(par)] <- "beta"
    }
  }
  if(is.null(gamma))
  {
    if(seasontype !="N")
    {
      if(m > 12)
        gamma <- 0.0002
      else
        gamma <- lower[3] + .01*(upper[3]-lower[3])
      if(gamma > 1-alpha)
        gamma <- min(0.999-alpha,0.001)
      if(gamma < lower[3] | gamma > upper[3])
        stop("Can't find consistent starting parameters")
      par <- c(par,gamma)
      names(par)[length(par)] <- "gamma"
    }
  }
  if(is.null(phi))
  {
    if(damped)
    {
      phi <- lower[4] + .99*(upper[4]-lower[4])
      par <- c(par,phi)
      names(par)[length(par)] <- "phi"
    }
  }
  
  return(par)
}

robinitstate <-  function(y, errortype, trendtype, seasontype)
{
  m <- frequency(y)
  n <- length(y)
  startup <- min(max(ceiling(10 / m) * m, 5 * m), floor(length(y) / m) * m)
  if (n < 2 * m)
    stop("Not enough data: need at least three periods.")

  if (trendtype == "N") {
    l0 <- median(y[1:startup])
    b0 <- NULL
    ydt <- l0
  } else{
    # trendtype = "A"
    fit <- myrm(y[1:startup]) # Repeated median
    l0 <- fit$level
    b0 <- fit$slope
    # If error type is "M", then we don't want l0+b0=0.
    # So perturb just in case.
    if (abs(l0 + b0) < 1e-8) {
      l0 <- l0 * (1 + 1e-3)
      b0 <- b0 * (1 - 1e-3)
    }
    ydt<-l0 + (1:startup) * b0
  }

  if (seasontype != "N")
  {
    # Seasonally adjusted data
    if (seasontype == "A") {
      init.seas <- rowMedians(matrix(y[1:startup] - ydt, nrow = m))
      if (errortype == "A")
        sigma0 <- mad(y[1:startup] - ydt - init.seas)
      else
        # if errortype = "M"
        sigma0 <- mad((y[1:startup] - ydt - init.seas) / (ydt - init.seas))
    } else{
      init.seas <- rowMedians(matrix(y[1:startup] / ydt, nrow = m))
      init.seas <-
        pmax(init.seas, 1e-2) # We do not want negative seasonal indexes
      if (errortype == "A")
        sigma0 <- mad(y[1:startup] - ydt / init.seas)
      else
        # if errortype = "M"
        sigma0 <- mad((y[1:startup] - ydt / init.seas) / (ydt / init.seas))
    }

    init.seas <- rev(init.seas) # reverse to match convention in ETS
    names(init.seas) <- paste("s", 0:(m - 1), sep = "")

    if (seasontype == "A") {
      l0 <- l0 + mean(init.seas)
      init.seas <- init.seas[1:(m - 1)] - mean(init.seas)
    }
    if (seasontype == "M") {
      me <- mean(init.seas)
      l0 <- l0 * me
      init.seas <- init.seas[1:(m - 1)] / me
    }
  } else
    # non-seasonal model
  {
    m <- 1
    init.seas <- NULL
    if (errortype == "A")
      sigma0 <- mad(y[1:startup] - ydt)
    else
      # if errortype = "M"
      sigma0 <- mad((y[1:startup] - ydt) / ydt)
  }

  names(l0) <- "l"
  if(!is.null(b0))
    names(b0) <- "b"
  names(sigma0) <- "sigma0"

  return(c(sigma0,l0,b0,init.seas))
}

initstate2 <- function(y,errortype,trendtype,seasontype)
{
  if(seasontype!="N")
  {
    # Do decomposition
    m <- frequency(y)
    n <- length(y)
    if(n < 4)
      stop("You've got to be joking (not enough data).")
    else if(n < 3*m) # Fit simple Fourier model.
    {
      fouriery <- fourier(y,1)
      fit <- tslm(y ~ trend + fouriery)
      if(seasontype=="A")
        y.d <- list(seasonal=y -fit$coef[1] - fit$coef[2]*(1:n))
      else # seasontype=="M". Biased method, but we only need a starting point
        y.d <- list(seasonal=y / (fit$coef[1] + fit$coef[2]*(1:n)))
    }
    else # n is large enough to do a decomposition
      y.d <- decompose(y,type=switch(seasontype, A="additive", M="multiplicative"))
    
    init.seas <- rev(y.d$seasonal[2:m]) # initial seasonal component
    names(init.seas) <- paste("s",0:(m-2),sep="")
    # Seasonally adjusted data
    if(seasontype=="A")
      y.sa <- y-y.d$seasonal
    else
    {
      init.seas <- pmax(init.seas, 1e-2) # We do not want negative seasonal indexes
      if(sum(init.seas) > m)
        init.seas <- init.seas/sum(init.seas + 1e-2)
      y.sa <- y/pmax(y.d$seasonal, 1e-2)
    }
  }
  else # non-seasonal model
  {
    m <- 1
    init.seas <- NULL
    y.sa <- y
  }
  
  maxn <- min(max(10,2*m),length(y.sa))
  
  if(trendtype=="N")
  {
    l0 <- mean(y.sa[1:maxn])
    b0 <- NULL
    sigma0 <- sqrt(mean((y.sa[1:maxn]-l0)^2))
  }
  else  # Simple linear regression on seasonally adjusted data
  {
    fit <- lsfit(1:maxn,y.sa[1:maxn])
    if(trendtype=="A")
    {
      l0 <- fit$coef[1]
      b0 <- fit$coef[2]
      # If error type is "M", then we don't want l0+b0=0.
      # So perturb just in case.
      if(abs(l0+b0) < 1e-8)
      {
        l0 <- l0*(1+1e-3)
        b0 <- b0*(1-1e-3)
      }
      sigma0 <- sqrt(mean(fit$residuals^2))
    }
    else #if(trendtype=="M")
    {
      l0 <- fit$coef[1]+fit$coef[2] # First fitted value
      if(abs(l0) < 1e-8)
        l0 <- 1e-7
      b0 <- (fit$coef[1] + 2*fit$coef[2])/l0 # Ratio of first two fitted values
      l0 <- l0/b0 # First fitted value divided by b0
      if(abs(b0) > 1e10) # Avoid infinite slopes
        b0 <- sign(b0)*1e10
      if(l0 < 1e-8 | b0 < 1e-8) # Simple linear approximation didn't work.
      {
        l0 <- max(y.sa[1],1e-3)
        b0 <- max(y.sa[2]/y.sa[1],1e-3)
      }
      sigma0 <- sqrt(mean((fit$resid/(y-fit$resid))^2))
    }
  }
  
  names(sigma0) <- "sigma0"
  names(l0) <- "l"
  if(!is.null(b0))
    names(b0) <- "b"
  return(c(sigma0,l0,b0,init.seas))
}

initstatebounds <- function(y,errortype,trendtype,seasontype){
  m <- frequency(y)
  
  up <- max(y)
  lo <- min(y)
  
  j=1
  if(trendtype !="N"){
    d <- diff(y)
    up[2] <- max(d)
    lo[2] <- min(d)
    j=j+1
  }
  if(seasontype == "A"){
    up[(j+1):(j+m-1)] <- rep(max(y)-min(y),m-1)
    lo[(j+1):(j+m-1)] <- -rep(max(y)-min(y),m-1)
  }
  if(seasontype == "M"){
    up[(j+1):(j+m-1)] <- rep(max(y)/min(y),m-1)
    lo[(j+1):(j+m-1)] <- rep(min(y)/max(y),m-1)
  }
  return(list(upper=up,lower=lo))
}

admissible <- function(alpha,beta,gamma,phi,m)
{
  if(is.null(phi))
    phi <- 1
  if(phi < 0 | phi > 1+1e-8)
    return(0)
  if(is.null(gamma))
  {
    if(alpha < 1-1/phi | alpha > 1+1/phi)
      return(0)
    if(!is.null(beta))
    {
      if(beta < alpha * (phi-1) | beta > (1+phi)*(2-alpha))
        return(0)
    }
  }
  else if(m > 1) # Seasonal model
  {
    if(is.null(beta))
      beta <- 0
    if(gamma < max(1-1/phi-alpha,0) | gamma > 1+1/phi-alpha)
      return(0)
    if(alpha < 1-1/phi-gamma*(1-m+phi+phi*m)/(2*phi*m))
      return(0)
    if(beta < -(1-phi)*(gamma/m+alpha))
      return(0)
    
    # End of easy tests. Now use characteristic equation
    P <- c(phi*(1-alpha-gamma),alpha+beta-alpha*phi+gamma-1,rep(alpha+beta-alpha*phi,m-2),(alpha+beta-phi),1)
    roots <- polyroot(P)
    
    #cat("maxpolyroots: ", max(abs(roots)), "\n")
    
    if(max(abs(roots)) > 1+1e-10)
      return(0)
  }
  # Passed all tests
  return(1)
}

getNewBounds <- function(par, lower, upper, nstate) {
  
  myLower <- NULL
  myUpper <- NULL
  
  if("alpha" %in% names(par)) {
    myLower <- c(myLower, lower[1])
    myUpper <- c(myUpper, upper[1])
  }
  if("beta" %in% names(par)) {
    myLower <- c(myLower, lower[2])
    myUpper <- c(myUpper, upper[2])
  }
  if("gamma" %in% names(par)) {
    myLower <- c(myLower, lower[3])
    myUpper <- c(myUpper, upper[3])
  }
  if("phi" %in% names(par)) {
    myLower <- c(myLower, lower[4])
    myUpper <- c(myUpper, upper[4])
  }
  
  myLower <- c(myLower,rep(-1e8,nstate))
  myUpper <- c(myUpper,rep(1e8,nstate))
  
  list(lower=myLower, upper=myUpper)
}

check.param <- function(alpha,beta,gamma,phi,lower,upper,bounds,m)
{
  if(bounds != "admissible")
  {
    if(!is.null(alpha))
    {
      if(alpha < lower[1] | alpha > upper[1])
        return(0)
    }
    if(!is.null(beta))
    {
      if(beta < lower[2] | beta > alpha | beta > upper[2])
        return(0)
    }
    if(!is.null(phi))
    {
      if(phi < lower[4] | phi > upper[4])
        return(0)
    }
    if(!is.null(gamma))
    {
      if(gamma < lower[3] | gamma > 1-alpha | gamma > upper[3])
        return(0)
    }
  }
  if(bounds != "usual")
  {
    if(!admissible(alpha,beta,gamma,phi,m))
      return(0)
  }
  return(1)
}

# Repeated Median
myrm <- function(x){
  l <- 1:(length(x))
  medquotient.t <- c()
  # median of the values in (x[k] - all the values that are diff to x[k]) divided by
  # the index k - all the indices diff than k, if k= 2 and l = 1 3 4 5 then k - l[l!=k]
  # is 1 -1 -2 -3
  # this is eq are in page 6 at the bottom
  for (k in l) medquotient.t <- c(medquotient.t, median((x[k] - x[l[l!=k]])/( k - l[l!=k])) ) 
  beta.RM <- median(medquotient.t)
  # this is alpha in the paper
  mu.RM <- median(x - beta.RM*l)
  res=x-mu.RM-beta.RM*l
  list(level = mu.RM, slope = beta.RM, res = res)
}

# Rowwise Medians
rowMedians <- function(x){
  apply(x,1,median)
}