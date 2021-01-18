
########### to find shape and rate of a Gamma distribution ############
getPrior_mu <- function(x, px, y, py, s1=1, s2=0.001, plot=TRUE){
  # x = quantile for p1
  # y = quantile for p2
  # p1 = the probability between 0 and x
  # p2 = the probability between 0 and y
  if (any(missing(x), missing(y), missing(px), missing(py))) stop("specify all arguments ('x','y','px','py')")
  if ((x<y) != (px<py)) stop("the arguments are not compatible, check their probability 'px' and 'py' arguments again.")
  f <- function(k){c(F1 = qgamma(px, k[1], k[2])-x, 
                     F2 = qgamma(py, k[1], k[2])-y)}
  ss <- multiroot(f = f, start = c(s1, s2))
  param <- rbind(ss$root)
  colnames(param) <- c("shape", "rate")
  if (plot) {
    musx <- seq(0, 1.5*max(x,y), length.out = 1000)
    musy <- dgamma(musx, shape = param[1], rate = param[2])
    plot(musx, musy, ylab = "Density", xlab = bquote(mu), type="l")
    polygon(c( 0, musx[musx<=x]),  c(musy[musx<=x],0 ), col="snow4")
    polygon(c( x, musx[musx<=y & musx>=x],x),  c(musy[musx<=y & musx>=x],0,0 ), col="snow3")
  }
  print(paste0("muPrior = list(priorDist = 'gamma', hyperpars = c(", round(param[1],3), ", ", round(param[2],3), "))"))
  return(invisible(param))
}

########### to find shapes of Beta distribution using confidence intervals #############
getPrior_delta <- function(lower, upper, p = 0.7, mode, conc, plot = TRUE){
  case1 <- case2 <- FALSE
  if (!missing(lower) & !missing(upper) & missing(mode) & missing(conc)) case1 <- TRUE;
  if (!missing(mode) & !missing(conc) & missing(lower) & missing(upper) ) case2 <- TRUE;
  if (case1 == case2) stop("specify arguments ('mode','conc') or ('lower','y') but not both")
  
  if (case1){
    # x = 2.5% quantile
    # y = 97.5% quantile
    # p = probability, area "most likely"
    f <- function(k){c(F1 = qbeta((1-p)/2, k[1], k[2])-lower, 
                       F2 = qbeta((1-p)/2, k[1], k[2], lower.tail = FALSE)-upper)}
    ss <- multiroot(f = f, start = c(1, 1))
    param <- rbind(ss$root)
  } else {
  # omega = assumed true mode
  # k = assumed true concentration
  alpha <- mode*(conc-2)+1
  beta <- (1-mode)*(conc-2)+1
  param <- cbind(alpha, beta)
  }
  colnames(param) <- c("alpha", "beta")
  denbeta <- function(x) dbeta(x, param[1], param[2])
  if (plot) {
    curve(denbeta, ylab = "Density", xlab = bquote(delta))
    polygon(c( lower, seq(lower, upper, length.out  = 101),lower), c(denbeta(seq(lower, upper, length.out  = 101)),0,0 ), col="snow4") 
  }
  print(paste0("deltaPrior = list(priorDist = 'beta', hyperpars = c(", round(param[1],3), ", ", round(param[2],3), "))"))
  return(invisible(param))
}
# shape and rate parameter for beta distribution where I believe 70% of the distribution lies between 0.5 and 0.85
# probability can be changed if practitioners are very sure or less sure about their CIs (default is 70%)
