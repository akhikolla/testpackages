#' This script simulates data using familial frailty model. 
#' We use the following variation: gamma(mu, ssq), where mu is the mean and ssq is sigma square.
#' See: https://www.rocscience.com/help/swedge/webhelp/swedge/Gamma_Distribution.htm
#' @param N Number of individuals.
#' @param f a list of formulas that define age (time) - dependency. Default: list(at="a", f1t="f1", Qt="Q*exp(theta*t)", ft="f", bt="b", mu0t="mu0*exp(theta*t)")
#' @param step An interval between two observations, a random uniformally-distributed value is then added to this step.
#' @param tstart Starting time (age).
#' Can be a number (30 by default) or a vector of two numbers: c(a, b) - in this case, starting value of time 
#' is simulated via uniform(a,b) distribution.
#' @param tend A number, defines final time (105 by default).
#' @param ystart A starting value of covariates.
#' @param sd0 A standard deviation for modelling the next covariate value, sd0 = 1 by default. 
#' @param nobs A number of observations (lines) for individual observations.
#' @param gamma_mu A parameter which is a mean value, default = 1
#' @param gamma_ssq A sigma squared, default = 0.5.
#' @return A table with simulated data.
#'@references Yashin, A. et al (2007), Health decline, aging and mortality: how are they related? 
#'Biogerontology, 8(3), 291-302.<DOI:10.1007/s10522-006-9073-3>.
#'@export
#' @examples
#' library(stpm)
#' dat <- simdata_gamma_frailty(N=10)
#' head(dat)
#'
simdata_gamma_frailty <- function(N=10,f=list(at="-0.05", f1t="80", Qt="2e-8", ft="80", bt="5", mu0t="1e-3"),
                             step=1, tstart=30, tend=105, ystart=80, sd0=1, nobs=NULL, gamma_mu=1, gamma_ssq=0.5) {
  
  simGamma <- function(gamma_mu, gamma_ssq) {
    #Generate a plot of gamma distributions. 
    # https://www.rocscience.com/help/swedge/webhelp/swedge/Gamma_Distribution.htm
    aplha <- (gamma_mu)^2/gamma_ssq # shape
    betta <- gamma_ssq/gamma_mu # scale
    simgamma <- dgamma(1, shape=aplha, scale = betta) #Make probability density function for default gamma distribution
    simgamma
  }
  
  formulas <- f  
  at <- NULL
  f1t <- NULL
  Qt <- NULL
  ft <- NULL
  bt <- NULL
  mu0t <- NULL
  
  comp_func_params <- function(astring, f1string, qstring, fstring, bstring, mu0string) {
    at <<- eval(bquote(function(t) .(parse(text = astring)[[1]])))
    f1t <<- eval(bquote(function(t) .(parse(text = f1string)[[1]]))) 
    Qt <<- eval(bquote(function(t) .(parse(text = qstring)[[1]])))
    ft <<- eval(bquote(function(t) .(parse(text = fstring)[[1]])))
    bt <<- eval(bquote(function(t) .(parse(text = bstring)[[1]])))
    mu0t <<- eval(bquote(function(t) .(parse(text = mu0string)[[1]])))
  }
  
  sigma_sq <- function(t1, t2) {
    # t2 = t_{j}, t1 = t_{j-1}
    ans <- bt(t1)*(t2-t1)
    ans
  }
  
  m <- function(y, t1, t2) {
    # y = y_{j-1}, t1 = t_{j-1}, t2 = t_{j}
    ans <- y + at(t1)*(y - f1t(t1))*(t2 - t1)
    ans
  }
  
  mu <- function(y, t, z) {
    ans <- z*mu0t(t) + (y - ft(t))^2*Qt(t)
  }
  
  
  comp_func_params(formulas$at, formulas$f1t, formulas$Qt, formulas$ft, formulas$bt, formulas$mu0t)
  
  data <- matrix(nrow=1,ncol=6, NA)
  record <- 1
  id <- 0
  famid <- 0
  for(i in 1:N/2) { # N must be even
    # Compute Z:
    Z <- simGamma(gamma_mu, gamma_ssq)
    for(ii in 1:2) {
      if(length(tstart) == 1) {
        t2 <- tstart # Starting time
      } else if(length(tstart) == 2){
        t2 <- runif(1,tstart[1], tstart[2]) # Starting time
      } else {
        stop(paste("Incorrect tstart:", tstart))
      }
    
      # Starting point
      new_person <- FALSE
      y2 <- rnorm(1,mean=ystart, sd=sd0) 
    
      n_observ <- 0
    
      while(new_person == FALSE){
        t1 <- t2
        t2 <- t1 + runif(1,-step/10,step/10) + step
        y1 <- y2
      
        S <- exp(-1*mu(y1, t1, Z)*(t2-t1))
      
        xi <- 0
        if (S > runif(1,0,1)) {
          xi <- 0
          y2 <- rnorm(1,mean=m(y1, t1, t2), sd=sqrt(sigma_sq(t1,t2)))
          new_person <- FALSE
          cov <- c(y1, y2)
          data <- rbind(data, c(id, xi, t1, t2, cov))
        } else {
          xi <- 1
          y2 <- NA
          new_person <- TRUE
          cov <- c(y1, y2)
          data <- rbind(data, c(id, xi, t1, t2, cov))
        }
      
        n_observ <- n_observ + 1
        if(!is.null(nobs)) {
          if(n_observ == nobs) {
            new_person <- TRUE;
          }
        }
      
        if(t2 > tend & new_person == FALSE) {
          new_person <- TRUE
        }
      }
      id <- id + 1
    }
    famid <- famid + 1
  }
  
  # One last step:
  data <- data[2:dim(data)[1],]
  colnames(data) <- c("id","xi","t1","t2", "y", "y.next")
  invisible(data.frame(data))
}