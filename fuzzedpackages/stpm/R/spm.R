#'A central function that estimates Stochastic Process Model parameters a from given dataset.
#'@references Yashin, A. et al (2007), Stochastic model for analysis of longitudinal data on aging 
#'and mortality. Mathematical Biosciences, 208(2), 538-551.
#'@references Akushevich I., Kulminski A. and Manton K. (2005). Life tables with covariates: Dynamic model 
#'for Nonlinear Analysis of Longitudinal Data. Mathematical Popu-lation Studies, 12(2), pp.: 51-80.
#'<DOI: 10.1080/08898480590932296>.
#'@references Yashin, A. et al (2007), Health decline, aging and mortality: how are they related? 
#'Biogerontology, 8(3), 291-302.<DOI:10.1007/s10522-006-9073-3>.
#'@param x A dataset: is the output from prepare_data(...) function and consists of two separate data tables:
#'(1) a data table for continuous-time model and (2) a data table for discrete-time model.
#'@param model A model type. Choices are: "discrete", "continuous" or "time-dependent".
#'@param formulas A list of parameter formulas used in the "time-dependent" model.
#'Default: \code{formulas=list(at="a", f1t="f1", Qt="Q", ft="f", bt="b", mu0t="mu0")}.
#'@param start A starting values of coefficients in the "time-dependent" model.
#'@param tol A tolerance threshold for matrix inversion (NULL by default).
#'@param stopifbound A flag (default=FALSE) if it is set then the optimization stops 
#'when any of the parametrs achives lower or upper boundary.
#'@param lb Lower boundary, default \code{NULL}.
#'@param ub Upper boundary, default \code{NULL}.
#'@param pinv.tol A tolerance threshold for matrix pseudo-inverse. Default: 0.01.
#'@param theta.range A user-defined range of the parameter \code{theta} used in 
#'discrete-time optimization and estimating of starting point for continuous-time optimization.
#'@param verbose A verbosing output indicator (FALSE by default).
#'@param gomp A flag (FALSE by default). When it is set, then time-dependent exponential form of mu0 and Q are used:
#' mu0 = mu0*exp(theta*t), Q = Q*exp(theta*t).
#'@param opts A list of options for \code{nloptr}.
#'Default value: \code{opt=list(algorithm="NLOPT_LN_NELDERMEAD", 
#'maxeval=100, ftol_rel=1e-8)}.
#'Please see \code{nloptr} documentation for more information.
#'@return For "discrete" (dmodel) and "continuous" (cmodel) model types: 
#'(1) a list of model parameter estimates for the discrete model type described in 
#'"Life tables with covariates: Dynamic Model for Nonlinear Analysis of Longitudinal Data", 
#'Akushevich et al, 2005.<DOI:10.1080/08898480590932296>,  and  
#'(2) a list of model parameter estimates for the continuous model type described in 
#'"Stochastic model for analysis of longitudinal data on aging and mortality", 
#'Yashin et al, 2007, Math Biosci.<DOI:10.1016/j.mbs.2006.11.006>.
#'
#'For the "time-dependent" model (model parameters depend on time): a set of model parameter estimates.
#'@export
#'@examples \dontrun{ 
#'library(stpm)
#'data.continuous <- simdata_cont(N=1000)
#'data.discrete <- simdata_discr(N=1000)
#'data <- list(data.continuous, data.discrete)
#'p.discr.model <- spm(data)
#'p.discr.model
#'p.cont.model <- spm(data, model="continuous")
#'p.cont.model
#'p.td.model <- spm(data, 
#'model="time-dependent",f=list(at="aa*t+bb", f1t="f1", Qt="Q", ft="f", bt="b", mu0t="mu0"), 
#'start=list(a=-0.001, bb=0.05, f1=80, Q=2e-8, f=80, b=5, mu0=1e-3))
#'p.td.model
#'}
spm <- function(x, model="discrete", 
                formulas = list(at="a", f1t="f1", Qt="Q", ft="f", bt="b", mu0t="mu0"), 
                start=NULL, tol=NULL, 
                stopifbound=FALSE, 
                lb=NULL, ub=NULL,
                pinv.tol = 0.01,
                theta.range=seq(0.01, 0.2, by=0.001),
                verbose=FALSE, gomp=FALSE,
                opts=list(algorithm="NLOPT_LN_NELDERMEAD", 
                          maxeval=100, ftol_rel=1e-8)) {
  
  # List of available models:
  models <- c("discrete", "continuous", "time-dependent")
  
  if(!(model %in% models)) {
    stop(cat(model, " - unknown model type!"))
  }
  
  # Number of variables (dimensions):
  k <- (dim(x[[1]])[2] - 4)/2
  
  
  if(model == "discrete") {
    # Estimation of starting point with discrete optimization:
    #pars <- spm_discrete(dat=x[[2]],verbose = verbose, tol = tol, theta_range=theta.range)
    pars <- spm_discrete(dat=x[[1]],verbose = verbose, tol = tol, theta_range=theta.range)
    res <- list(dmodel=list(u=pars$dmodel$u, 
                            R=pars$dmodel$R, 
                            b=pars$dmodel$b, 
                            Q=pars$dmodel$Q, 
                            Sigma=pars$dmodel$Sigma,
                            mu0=pars$dmodel$mu0,
                            theta=pars$dmodel$theta), 
                cmodel=list(a=pars$cmodel$a, 
                            f1=pars$cmodel$f1,
                            Q=pars$cmodel$Q,
                            f=pars$cmodel$f, 
                            b=pars$cmodel$b, 
                            mu0=pars$cmodel$mu0, 
                            theta=pars$cmodel$theta))
    
  }
  
  
  if(model == "continuous") {
    #pars <- spm_discrete(dat=x[[2]],verbose = verbose, tol = tol)
    pars <- spm_discrete(dat=x[[1]],verbose = verbose, tol = tol)
    data <- data.frame(x[[1]])
  
    if(verbose) {
      cat("Starting parameters:\n")
      print(pars)
    }
    
    if(det(pars$cmodel$Q) < 0) {
      cat("Error: determinant of Q < 0\n")
      cat("Q:\n")
      print(pars$cmodel$Q)
      cat("Det(Q):\n")
      print(det(pars$cmodel$Q))
      
      res <- NA
    
    } else {
      res.t <- spm_continuous(as.matrix(data), 
                    a=pars$cmodel$a, 
                    f1=pars$cmodel$f1, 
                    Q=pars$cmodel$Q, 
                    f=pars$cmodel$f, 
                    b=pars$cmodel$b, 
                    mu0=pars$cmodel$mu0, 
                    theta=pars$cmodel$theta, 
                    stopifbound = stopifbound,
                    lb = lb, ub = ub,
                    pinv.tol = pinv.tol,
                    verbose = verbose, 
                    gomp=gomp, 
                    opts=opts)
      
      Q.c <- res.t$Q
      R.c <- res.t$a + diag(k)
      Sigma.c <- as.matrix(res.t$b)
      u.c <- (-1)*(t(res.t$f1) %*% res.t$a)
      b.c <- -2*t(res.t$f) %*% res.t$Q
      mu0.c <- res.t$mu0 + t(res.t$f) %*% res.t$Q %*% res.t$f
      theta.c <- res.t$theta
      
      res <- list(dmodel=list(u=u.c, 
                              R=R.c, 
                              b=b.c, 
                              Q=Q.c, 
                              Sigma=Sigma.c,
                              mu0=mu0.c,
                              theta=theta.c), 
                  cmodel=list(a=res.t$a, 
                              f1=res.t$f1,
                              Q=res.t$Q,
                              f=res.t$f, 
                              b=res.t$b, 
                              mu0=res.t$mu0, 
                              theta=res.t$theta))
      
      #print(res.t)
      
    }
  }
  
  if(model == "time-dependent") {
    
    if(k > 1) {
        stop("Number of variables > 1. Model with time-dependent parameters can be used only with one variable!")
    }
    
    #if(length(formulas) != 6) {
    #  stop("It must be 6 equations for corresponding coefficients.")
    #}
    
    # Raw parameters estimates
    #pars <- spm_discrete(dat=x[[2]],verbose = verbose, tol = tol, theta_range=theta.range)
    # Parameter optimization for time-dependent model
    if(is.null(start)) {
        warning("Default starting values will be used.")
    }
    
    res <- spm_time_dep(x[[1]], 
                        frm=formulas,
                        start=start,
                        lb=lb, ub=ub,
                        verbose=verbose, 
                        opts = opts)
    
  }
  class(res) <- "spm"
  invisible(res)
}