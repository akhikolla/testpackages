#'Multi-dimension simulation function for data with partially observed covariates (multidimensional GenSPM) with arbitrary intervals
#'@references Arbeev, K.G. et al (2009). Genetic model for longitudinal studies of aging, health, and longevity
# and its potential application to incomplete data. Journal of Theoretical
# Biology 258(1), 103{111 (2009).<doi:10.1016/j.jtbi.2009.01.023>
#'@references Yashin, A.I. et al (2007). Stochastic model for analysis of longitudinal data on aging 
#'and mortality. Mathematical Biosciences, 208(2), 538-551.<DOI:10.1016/j.mbs.2006.11.006>.
#' @param N Number of individuals.
#' @param aH A k by k matrix, which characterize the rate of the adaptive response when Z = 1.
#' @param aL A k by k matrix, which characterize the rate of the adaptive response when Z = 0.
#' @param f1H A particular state, which if a deviation from the normal (or optimal) when Z = 1.
#' This is a vector with length of k.
#' @param f1L A particular state, which if a deviation from the normal (or optimal) when Z = 0. 
#' This is a vector with length of k.
#' @param QH A matrix k by k, which is a non-negative-definite symmetric matrix when Z = 1.
#' @param QL A matrix k by k, which is a non-negative-definite symmetric matrix when Z = 0.
#' @param fH A vector-function (with length k) of the normal (or optimal) state when Z = 1.
#' @param fL A vector-function (with length k) of the normal (or optimal) state when Z = 0.
#' @param bH A diffusion coefficient, k by k matrix when Z = 1.
#' @param bL A diffusion coefficient, k by k matrix when Z = 0.
#' @param mu0H mortality at start period of time when Z = 1.
#' @param mu0L mortality at start period of time when Z = 0.
#' @param thetaH A displacement coefficient of the Gompertz function when Z = 1.
#' @param thetaL A displacement coefficient of the Gompertz function when Z = 0.
#' @param p A proportion of carriers in a sumulated population (default p = 0.25).
#' @param ystart A vector with length equal to number of dimensions used, defines starting values of covariates.
#' @param dt A discrete step size between two observations. A random uniform value is then added to this step size.
#' @param tstart A number that defines starting time (30 by default).
#' @param tend A number, defines final time (105 by default).
#' @param sd0 A standard deviation for modelling the next physiological variable (covariate) value.
#' @param mode Can have the following values: "observed" (default), "unobserved".
#' This represents a type of group to simulate: a group with observed variable Z, or group with unbobserved variable Z.
#' @param gomp A flag (FALSE by default). When it is set, then time-dependent exponential form of mu0 and Q are used:
#' mu0 = mu0*exp(theta*t).
#' @param nobs A number of observations (lines) for individual observations.
#' @return A table with simulated data.
#' @export
#' @examples
#' library(stpm)
#' dat <- sim_pobs(N=50)
#' head(dat)
#'
sim_pobs <- function(N=10, 
                          aH=-0.05, aL=-0.01, 
                          f1H=60, f1L=80, 
                          QH=2e-8, QL=2.5e-8, 
                          fH=60, fL=80, 
                          bH=4, bL=5, 
                          mu0H=0.8e-5, mu0L=1e-5, 
                          thetaH=0.08, thetaL=0.1, 
                          p=0.25,
                          ystart=80, 
                          tstart=30, 
                          tend=105, 
                          dt=1, 
                          sd0=1,
                          mode="observed", 
                          gomp=FALSE, 
                          nobs = NULL) {
  
  k <- length(ystart)
  
  if ( (dim(as.data.frame(aH))[1] != k) & (dim(as.data.frame(aH))[2] != k) &
       (dim(as.data.frame(QH))[1] != k) & (dim(as.data.frame(QH))[2] != k) & 
       (dim(as.data.frame(f1H))[1] != k) & (dim(as.data.frame(fH))[1] != k) &
       (dim(as.data.frame(bH))[1] != k) & 
       (dim(as.data.frame(ystart))[1] != k) ) {
    stop("Dimenstions of provided parameters are not equal.")
  }  
  
  if ( (dim(as.data.frame(aL))[1] != k) & (dim(as.data.frame(aL))[2] != k) &
       (dim(as.data.frame(QL))[1] != k) & (dim(as.data.frame(QL))[2] != k) & 
       (dim(as.data.frame(f1L))[1] != k) & (dim(as.data.frame(fL))[1] != k) &
       (dim(as.data.frame(bL))[1] != k) & 
       (dim(as.data.frame(ystart))[1] != k) ) {
    stop("Dimenstions of provided parameters are not equal.")
  }  
  
  aH<-matrix(aH,nrow=k,ncol=k); aL<-matrix(aL,nrow=k,ncol=k)
  f1H<-matrix(f1H,nrow=k,ncol=1,byrow=FALSE); f1L<-matrix(f1L,nrow=k,ncol=1,byrow=FALSE)
  QH<-matrix(QH,nrow=k,ncol=k); QL<-matrix(QL,nrow=k,ncol=k)
  fH<-matrix(fH,nrow=k,ncol=1,byrow=FALSE); fL<-matrix(fL,nrow=k,ncol=1,byrow=FALSE)
  bH<-matrix(bH,nrow=k,ncol=1,byrow=FALSE); bL<-matrix(bL,nrow=k,ncol=1,byrow=FALSE)
  
  ystart<-matrix(ystart,nrow=k,ncol=1,byrow=FALSE)
  
  genmode <- 1
  if(mode == "observed") {
    genmode <- 1
  } else if(mode == "unobserved") {
    genmode <- 0
  } else {
    stop("Mode can be set to 'observed' or 'unobserved'.")
  }
  
  if(is.null(nobs)) {
    nobs <- 0
  }
  
  if(length(tstart) > 2) {
    stop(paste("Incorrect tstart:", tstart))
  }
  
  simulated = .Call("simGenCont", 
                    N, 
                    aH, aL, 
                    f1H, f1L, 
                    QH, QL, 
                    fH, fL, 
                    bH, bL, 
                    mu0H, mu0L,
                    thetaH, thetaL,
                    p,
                    tstart, ystart, tend, 
                    k, dt, sd0, genmode, gomp, nobs);
  
  data_names <- c()
  for(n in 1:k) {
    data_names <- c(data_names, paste("y",n, sep=''), paste("y",n, ".next", sep=''))
  }
  
  if(genmode == 1) {
    colnames(simulated) <- c("id", "xi", "t1", "t2", "Z", data_names)
  } else {
    colnames(simulated) <- c("id", "xi", "t1", "t2", data_names)
  }
  
  invisible(data.frame(simulated))
}

