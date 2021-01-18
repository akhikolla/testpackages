#' Multi-dimensional simulation function for continuous-time SPM.
#'@references Yashin, A.I. et al (2007). Stochastic model for analysis of longitudinal data on aging 
#'and mortality. Mathematical Biosciences, 208(2), 538-551.<DOI:10.1016/j.mbs.2006.11.006>.
#' @param N Number of individuals.
#' @param a A k by k matrix, represents the adaptive capacity of the organism
#' @param f1 A trajectory that corresponds to the long-term average value of the stochastic process Y(t), 
#' which describes a trajectory of individual covariate (physiological variable) influenced by different 
#' factors represented by a random Wiener process W(t).
#' This is a vector with length of k.
#' @param Q A matrix k by k, which is a non-negative-definite symmetric matrix, 
#' represents a sensitivity of risk function to deviation from the norm.
#' @param f A vector with length of k, represents the normal (or optimal) state of physiological variable.
#' @param b A diffusion coefficient, k by k matrix, 
#' characterizes a strength of the random disturbances from Wiener process W(t).
#' @param mu0 A baseline mortality.
#' @param theta A displacement coefficient.
#' @param ystart A vector with length equal of k, defines starting values of covariates.
#' @param dt A discrete step size between two observations. A random uniform value is then added to this step size.
#' @param tstart A number that defines starting time (30 by default).
#' @param tend A number, defines final time (105 by default).
#' @param sd0 a standard deviation for modelling the next covariate value.
#' @param gomp A flag (FALSE by default). When it is set, then time-dependent exponential form of mu0 and Q are used:
#' mu0 = mu0*exp(theta*t).
#' @param nobs A number of observations (lines) for individual observations.
#' @param format Data format: "long" (default), "short".
#' @return A table with simulated data.
#' @export
#' @examples
#' library(stpm)
#' dat <- simdata_cont(N=50)
#' head(dat)
#'
simdata_cont <- function(N=10, 
                          a=-0.05, 
                          f1=80, 
                          Q=2e-8, 
                          f=80, 
                          b=5, 
                          mu0=1e-5, 
                          theta=0.08, 
                          ystart=80, 
                          tstart=30, 
                          tend=105, 
                          dt=1, 
                          sd0=1,
                          nobs=NULL, 
                          gomp=TRUE,
                         format="long") {
  
  k <- length(ystart)
  
  if ( (dim(as.data.frame(a))[1] != k) & (dim(as.data.frame(a))[2] != k) &
       (dim(as.data.frame(Q))[1] != k) & (dim(as.data.frame(Q))[2] != k) & 
       (dim(as.data.frame(f1))[1] != k) & (dim(as.data.frame(f))[1] != k) &
       (dim(as.data.frame(b))[1] != k) & 
       (dim(as.data.frame(ystart))[1] != k) ) {
    stop("Dimenstions of provided parameters are not equal.")
  }  
  
  if ( (dim(as.data.frame(a))[1] != k) & (dim(as.data.frame(a))[2] != k) &
       (dim(as.data.frame(Q))[1] != k) & (dim(as.data.frame(Q))[2] != k) & 
       (dim(as.data.frame(f1))[1] != k) & (dim(as.data.frame(f))[1] != k) &
       (dim(as.data.frame(b))[1] != k) & 
       (dim(as.data.frame(ystart))[1] != k) ) {
    stop("Dimenstions of provided parameters are not equal.")
  }  
  
  aH<-matrix(a,nrow=k,ncol=k)
  f1H<-matrix(f1,nrow=k,ncol=1,byrow=FALSE)
  QH<-matrix(Q,nrow=k,ncol=k)
  fH<-matrix(f,nrow=k,ncol=1,byrow=FALSE)
  bH<-matrix(b,nrow=k,ncol=1,byrow=FALSE)
  ystart<-matrix(ystart,nrow=k,ncol=1,byrow=FALSE)
  
  simulated = .Call("simCont", N, aH, f1H, QH, fH, bH, mu0, theta, tstart, ystart, tend, k, dt, sd0, nobs, gomp);
  
  data_names <- c()
  for(n in 1:k) {
    data_names <- c(data_names, paste("y",n, sep=''), paste("y",n, ".next", sep=''))
  }
  colnames(simulated) <- c("id", "xi", "t1", "t2", data_names)
  simulated <- data.frame(simulated)
  
  if(format == "short")
  {
    simulated <- make.short.format(simulated)
  }
  
  invisible(return(simulated))
}

