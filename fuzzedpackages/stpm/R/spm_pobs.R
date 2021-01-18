#'Continuous-time multi-dimensional optimization for SPM with partially observed covariates (multidimensional GenSPM)
#'@references Arbeev, K.G. et al (2009). Genetic model for longitudinal studies of aging, health, and longevity
# and its potential application to incomplete data. Journal of Theoretical
# Biology 258(1), 103{111 (2009).<doi:10.1016/j.jtbi.2009.01.023>
#'@references Yashin, A.I. et al (2007). Stochastic model for analysis of longitudinal data on aging 
#'and mortality. Mathematical Biosciences, 208(2), 538-551.<DOI:10.1016/j.mbs.2006.11.006>.
#'@param x A data table with genetic component.
#'@param y A data table without genetic component.
#'@param aH A k by k matrix. Characterizes the rate of the adaptive response for Z = 1.
#'@param aL A k by k matrix. Characterize the rate of the adaptive response for Z = 0.
#'@param f1H A deviation from the norm (or optimal) state for Z = 1.
#'This is a vector of length k.
#'@param f1L A deviation from the norm (or optimal) for Z = 0. 
#'This is a vector of length k.
#'@param QH A matrix k by k, which is a non-negative-definite symmetric matrix for Z = 1.
#'@param QL A matrix k by k, which is a non-negative-definite symmetric matrix for Z = 0.
#'@param fH A vector with length of k. Represents the normal (or optimal) state for Z = 1.
#'@param fL A vector with length of k. Represents the normal (or optimal) state for Z = 0.
#'@param bH A diffusion coefficient, k by k matrix for Z = 1.
#'@param bL A diffusion coefficient, k by k matrix for Z = 0.
#'@param mu0H A baseline mortality for Z = 1.
#'@param mu0L A baseline mortality for Z = 0.
#'@param thetaH A displacement coefficient for Z = 1.
#'@param thetaL A displacement coefficient for Z = 0.
#'@param p a hyphotetical percentage of presence of partially observed covariate in a population (default p=0.25).
#'@param stopifbound If TRUE then estimation stops if at least one parameter achieves lower or upper boundaries.
#'@param algorithm An optimization algorithm used, can be one of those provided by \code{nloptr}. 
#'#'Check the NLopt website for a description of
#'the algorithms. Default: NLOPT_LN_NELDERMEAD
#'@param lb Lower bound of parameter values.
#'@param ub Upper bound of parameter values.
#'@param maxeval Maximum number of iterations of the algorithm for \code{nloptr} optimization. 
#'The program stops when the number of function evaluations exceeds maxeval. Default: 500.
#'@param verbose An indicator of verbosing output (FALSE by default).
#'@param pinv.tol A tolerance value for pseudo-inverse of matrix gamma (see Yashin, A.I. et al (2007). Stochastic model for analysis of longitudinal data on aging 
#'and mortality. Mathematical Biosciences, 208(2), 538-551.<DOI:10.1016/j.mbs.2006.11.006>.)
#'@param mode Can be one of the following: "observed" (default), "unobserved" or "combined". 
#'mode = "observed" represents analysing only dataset with observed variable Z.
#'mode = "unobserved" represents analysing only dataset of unobserved variable Z.
#'mode = "combined" denoted joint analysis of both observed and unobserved datasets.
#'@param gomp A flag (FALSE by default). When it is set, then time-dependent exponential form of mu0 is used:
#' mu0 = mu0*exp(theta*t).
#'@param ftol_rel Relative tolerance threshold for likelihood function (defalult: 1e-6), 
#'see http://ab-initio.mit.edu/wiki/index.php/NLopt_Reference
#'@return A set of estimated parameters aH, aL, f1H, f1H, QH, QL, fH, fL, bH, bL, mu0H, mu0L, thetaH, thetaL, p and
#'additional variable \code{limit} which indicates if any parameter 
#'achieved lower or upper boundary conditions (FALSE by default).
#'@export
#'@examples \dontrun{
#'library(stpm)
#'#Reading the data:
#'data <- sim_pobs(N=1000)
#'head(data)
#'#Parameters estimation:
#'pars <- spm_pobs(x=data)
#'pars
#'}
spm_pobs <- function(x=NULL, y=NULL,
                    aH=-0.05, aL=-0.01, 
                    f1H=60, f1L=80, 
                    QH=2e-8, QL=2.5e-8, 
                    fH=60, fL=80, 
                    bH=4, bL=5, 
                    mu0H=0.8e-5, mu0L=1e-5, 
                    thetaH=0.08, thetaL=0.1,
                    p=0.25,
                    stopifbound=FALSE, 
                    algorithm="NLOPT_LN_NELDERMEAD",
                    lb=NULL, ub=NULL,
                    maxeval=500,
                    verbose=FALSE,
                    pinv.tol=0.01,
                    mode="observed",
                    gomp=TRUE,
                    ftol_rel=1.0e-6) {
  
  setlb <- function(k, params) {
    # This function sets lower and upper boundaries for optim.
    # - k - number of dimensions
    # - params - initial parameters, a vector
    #
    # Lower boundaries:
    lower_bound <- c()
    # Setting boundaries for coefficients:
    # aH, aL
    start=1; end=k^2
    lower_bound <- c(lower_bound, unlist(lapply(start:end, function(n){ params[n] + ifelse(params[n] > 0, -0.1*params[n], 0.1*params[n]) }))) 
    start=end+1; end=start+k^2-1
    lower_bound <- c(lower_bound, unlist(lapply(start:end, function(n){ params[n] + ifelse(params[n] > 0, -0.1*params[n], 0.1*params[n]) }))) 
    
    # f1H, f1L
    start=end+1; end=start+k-1
    lower_bound <- c(lower_bound, unlist(lapply(start:end, function(n){ params[n] + ifelse(params[n] > 0, -0.1*params[n], 0.1*params[n]) }))) 
    start=end+1; end=start+k-1
    lower_bound <- c(lower_bound, unlist(lapply(start:end, function(n){ params[n] + ifelse(params[n] > 0, -0.1*params[n], 0.1*params[n]) }))) 
    
    # QH, QL
    start=end+1; end=start+k^2-1
    lower_bound <- c(lower_bound, unlist(lapply(start:end, function(n){ ifelse(params[n] > 0, 0, params[n]+0.1*params[n]) })))
    start=end+1; end=start+k^2-1
    lower_bound <- c(lower_bound, unlist(lapply(start:end, function(n){ ifelse(params[n] > 0, 0, params[n]+0.1*params[n]) })))
    
    # fH, fL
    start=end+1; end=start+k-1
    lower_bound <- c(lower_bound, unlist(lapply(start:end, function(n){ params[n] + ifelse(params[n] > 0, -0.1*params[n], 0.1*params[n]) })) )
    start=end+1; end=start+k-1
    lower_bound <- c(lower_bound, unlist(lapply(start:end, function(n){ params[n] + ifelse(params[n] > 0, -0.1*params[n], 0.1*params[n]) })) )
    
    # bH, bL
    start=end+1; end=start+k-1
    lower_bound <- c(lower_bound, unlist(lapply(start:end, function(n){params[n] + ifelse(params[n] > 0, -0.1*params[n], 0.1*params[n]) })) )
    start=end+1; end=start+k-1
    lower_bound <- c(lower_bound, unlist(lapply(start:end, function(n){params[n] + ifelse(params[n] > 0, -0.1*params[n], 0.1*params[n]) })) )
    
    # mu0h, mu0L
    start=end+1; end=start
    lower_bound <- c( lower_bound, params[start:end] - 0.1*params[start:end])
    start=end+1; end=start
    lower_bound <- c( lower_bound, params[start:end] - 0.1*params[start:end])
    
    # thetah, thetal
    start=end+1; end=start
    lower_bound <- c( lower_bound, params[start:end] - 0.1*params[start:end])
    start=end+1; end=start
    lower_bound <- c( lower_bound, params[start:end] - 0.1*params[start:end])
    
    # p
    start=end+1; end=start
    lower_bound <- c( lower_bound, ifelse((params[start:end] - 0.1*params[start:end]) <= 0,
                                          0,
                                          params[start:end] - 0.1*params[start:end]))
    
    lower_bound
  }
  
  setub <- function(k, params) {
    # This function sets lower and upper boundaries for optim.
    # - k - number of dimensions
    # - params - initial parameters, a vector
    #
    #Upper boundaries:
    upper_bound <- c()
    # Setting boundaries for coefficients:
    # aH, aL
    start=1; end=k^2
    upper_bound <- c(upper_bound, unlist(lapply(start:end, function(n){ifelse(params[n] > 0, params[n] + 0.1*params[n], 0*params[n]) })))
    start=end+1; end=start+k^2-1
    upper_bound <- c(upper_bound, unlist(lapply(start:end, function(n){ifelse(params[n] > 0, params[n] + 0.1*params[n], 0*params[n]) })))
    
    # f1H, f1L
    start=end+1; end=start+k-1
    upper_bound <- c(upper_bound, unlist(lapply(start:end, function(n){ifelse(params[n] > 0, params[n] + 0.1*params[n], params[n]- 0.1*params[n]) })))
    start=end+1; end=start+k-1
    upper_bound <- c(upper_bound, unlist(lapply(start:end, function(n){ifelse(params[n] > 0, params[n] + 0.1*params[n], params[n]- 0.1*params[n]) })))
    
    # QH, QL
    start=end+1; end=start+k^2-1
    upper_bound <- c(upper_bound, unlist(lapply(start:end, function(n) {ifelse(params[n] > 0, params[n] + 0.1*params[n], params[n]-0.1*params[n] )})))
    start=end+1; end=start+k^2-1
    upper_bound <- c(upper_bound, unlist(lapply(start:end, function(n) {ifelse(params[n] > 0, params[n] + 0.1*params[n], params[n]-0.1*params[n] )})))
    
    # fH, fL
    start=end+1; end=start+k-1
    upper_bound <- c(upper_bound, unlist(lapply(start:end, function(n){ifelse(params[n] > 0, params[n]+0.1*params[n], params[n]-0.1*params[n]) })))
    start=end+1; end=start+k-1
    upper_bound <- c(upper_bound, unlist(lapply(start:end, function(n){ifelse(params[n] > 0, params[n]+0.1*params[n], params[n]-0.1*params[n]) })))
    
    # bH, bL
    start=end+1; end=start+k-1
    upper_bound <- c(upper_bound, unlist(lapply(start:end, function(n){ifelse(params[n] > 0, params[n]+0.1*params[n], params[n]-0.1*params[n]) })))
    start=end+1; end=start+k-1
    upper_bound <- c(upper_bound, unlist(lapply(start:end, function(n){ifelse(params[n] > 0, params[n]+0.1*params[n], params[n]-0.1*params[n]) })))
    
    # mu0h, mu0L
    start=end+1; end=start
    upper_bound <- c( upper_bound, params[start:end] + 0.1*params[start:end])
    start=end+1; end=start
    upper_bound <- c( upper_bound, params[start:end] + 0.1*params[start:end])
    
    # thetah, thetal
    start=end+1; end=start
    upper_bound <- c( upper_bound, params[start:end] + 0.1*params[start:end])
    start=end+1; end=start
    upper_bound <- c( upper_bound, params[start:end] + 0.1*params[start:end])
    
    # p
    start=end+1; end=start
    upper_bound <- c( upper_bound, 
                      ifelse((params[start:end] + 0.1*params[start:end]) >= 1, 
                             1, 
                             (params[start:end] + 0.1*params[start:end])))
    
    
    upper_bound
  }
  
  
  
  
  ###=======For DEBUG========###
  #dat = dat
  #
  #a=matrix(c(-0.05,  0.001, 0.001, -0.05), nrow = 2, ncol = 2, byrow = T)
  #f1=t(matrix(c(100, 200), nrow = 2, ncol = 1, byrow = F))
  #Q=matrix(c(1e-06, 1e-7, 1e-7,  1e-06), nrow = 2, ncol = 2, byrow = T)
  #f=t(matrix(c(100, 200), nrow = 2, ncol = 1, byrow = F))
  #b=matrix(c(2, 5), nrow = 2, ncol = 1, byrow = F)
  #mu0=1e-4
  #theta=0.08
  #k=2
  #
  #stopifbound=FALSE
  #algorithm="NLOPT_LN_NELDERMEAD"
  #lb=NULL
  #ub=NULL
  #verbose=FALSE
  ###=========================###
  
  
  
  avail_algorithms <- c("NLOPT_GN_DIRECT", "NLOPT_GN_DIRECT_L",
                        "NLOPT_GN_DIRECT_L_RAND", "NLOPT_GN_DIRECT_NOSCAL",
                        "NLOPT_GN_DIRECT_L_NOSCAL",
                        "NLOPT_GN_DIRECT_L_RAND_NOSCAL",
                        "NLOPT_GN_ORIG_DIRECT", "NLOPT_GN_ORIG_DIRECT_L",
                        "NLOPT_GD_STOGO", "NLOPT_GD_STOGO_RAND",
                        "NLOPT_LD_SLSQP", "NLOPT_LD_LBFGS_NOCEDAL",
                        "NLOPT_LD_LBFGS", "NLOPT_LN_PRAXIS", "NLOPT_LD_VAR1",
                        "NLOPT_LD_VAR2", "NLOPT_LD_TNEWTON",
                        "NLOPT_LD_TNEWTON_RESTART",
                        "NLOPT_LD_TNEWTON_PRECOND",
                        "NLOPT_LD_TNEWTON_PRECOND_RESTART",
                        "NLOPT_GN_CRS2_LM", "NLOPT_GN_MLSL", "NLOPT_GD_MLSL",
                        "NLOPT_GN_MLSL_LDS", "NLOPT_GD_MLSL_LDS",
                        "NLOPT_LD_MMA", "NLOPT_LN_COBYLA", "NLOPT_LN_NEWUOA",
                        "NLOPT_LN_NEWUOA_BOUND", "NLOPT_LN_NELDERMEAD",
                        "NLOPT_LN_SBPLX", "NLOPT_LN_AUGLAG", "NLOPT_LD_AUGLAG",
                        "NLOPT_LN_AUGLAG_EQ", "NLOPT_LD_AUGLAG_EQ",
                        "NLOPT_LN_BOBYQA", "NLOPT_GN_ISRES")
  
  if(!(algorithm %in% avail_algorithms)) {
    stop(cat("Provided algorithm", algorithm, "not in the list of available optimization methods."))
  }
  
  mode_avail <- c("observed", "unobserved", "combined")
  if(!(mode %in% mode_avail)){
    stop(cat("Provided mode ", mode, " not found in a set of available modes: ", mode_avail))
  }
  
  if(verbose)
    cat("Provided mode: ", mode, "\n")
  
  if((mode == "observed" | mode == "combined") & !is.null(x)) {
    dat <- as.matrix(x[, 2:dim(x)[2]])
  }
  
  if( !is.null(y) & (mode == "unobserved" | mode == "combined") ) {
    #cat("Provided mode: ", mode, "\n")
    unobsdat <- as.matrix(y[, 2:dim(y)[2]])
  } 
  
  if(!(mode %in% mode_avail)) {
    stop(cat("Provided mode:", mode, " is unknown or provided data is null."))
  }
  
  
  k <- dim(as.matrix(aH))[1]
  final_res <- list()
  
  if(mu0H < 0) {mu0H <- 0}
  if(mu0L < 0) {mu0L <- 0}
  
  #######################################################################################################
  if(mode == "observed" | mode == "combined") {
    dd <- cbind(dat[,1], dat[,5])
    colnames(dd) <- c("id", "Z")
    ddd <- data.frame(dd)
    d<-aggregate(Z ~ id, data=ddd, FUN=mean)
    N.c <- length(which(d$Z == 1)) # Carriers
    N.nc <- length(which(d$Z == 0)) # Non-carriers
  }
  #######################################################################################################
  
  parameters <- c(t(aH), t(aL), f1H, f1L, t(QH), t(QL), fH, fL, bH, bL, mu0H, mu0L, thetaH, thetaL, p)
  # Current results:
  results_tmp <- list(aH=NULL, aL=NULL, 
                      f1H=NULL, f1L=NULL,
                      QH=NULL, QL=NULL,
                      fH=NULL, fL=NULL,
                      bH=NULL, bL=NULL,
                      mu0H=NULL, mu0L=NULL,
                      thetaH=NULL, thetaL=NULL,
                      p=NULL)
  iteration <- 0
  
  bounds <- list()
  
  if(is.null(lb)) {
    bounds$lower_bound <- setlb(k, parameters)
  } else {
    bounds$lower_bound <- lb
  }
  
  if(is.null(ub)) {
    bounds$upper_bound <- setub(k, parameters)
  } else {
    bounds$upper_bound <- ub
  }
  
  # Reading parameters:
  #start=1; end=k^2
  #a <- matrix(parameters[start:end],ncol=k, nrow=k, byrow=T)
  #results$a <- a
  ##print(results$a)
  #start=end+1; end=start+k-1
  #f1 <- matrix(parameters[start:end],ncol=1, nrow=k, byrow=T)
  #results$f1 <- f1
  ##print(results$f1)
  #start=end+1; end=start+k^2-1
  #Q <- matrix(parameters[start:end],ncol=k, nrow=k, byrow=T)
  #results$Q <- Q
  ##print(results$Q)
  #start=end+1; end=start+k-1
  #f <- matrix(parameters[start:end],ncol=1, nrow=k, byrow=F)
  #results$f <- f
  ##print(results$f)
  #start=end+1; end=start+k-1
  #b <- matrix(parameters[start:end],nrow=k, ncol=1, byrow=F)
  #results$b <- b
  ##print(results$b)
  #start=end+1; end=start
  #mu0 <- parameters[start:end]
  #results$mu0 <- mu0
  ##print(results$mu0)
  #start=end+1; end=start
  #theta <- parameters[start:end]
  #results$theta <- theta
  ##print(results$theta)
  ## End of reading parameters
  
  maxlik <- function(par) {
    
    stopflag <- F
    # Reading parameters:
    
    # ah, al
    start=1; end=k^2
    ah <- matrix(par[start:end],ncol=k, nrow=k, byrow=TRUE)
    results_tmp$aH <<- ah
    start=end+1; end=start+k^2-1
    al <- matrix(par[start:end],ncol=k, nrow=k, byrow=TRUE)
    results_tmp$aL <<- al
    
    # f1h, f1l
    start=end+1; end=start+k-1
    f1h <- matrix(par[start:end],ncol=1, nrow=k, byrow=FALSE)
    results_tmp$f1H <<- f1h
    start=end+1; end=start+k-1
    f1l <- matrix(par[start:end],ncol=1, nrow=k, byrow=FALSE)
    results_tmp$f1L <<- f1l
    
    # Qh, Ql
    start=end+1; end=start+k^2-1
    Qh <- matrix(par[start:end],ncol=k, nrow=k, byrow=TRUE)
    results_tmp$QH <<- Qh
    start=end+1; end=start+k^2-1
    Ql <- matrix(par[start:end],ncol=k, nrow=k, byrow=TRUE)
    results_tmp$QL <<- Ql
    
    # fh, fl
    start=end+1; end=start+k-1
    fh <- matrix(par[start:end],ncol=1, nrow=k, byrow=FALSE)
    results_tmp$fH <<- fh
    start=end+1; end=start+k-1
    fl <- matrix(par[start:end],ncol=1, nrow=k, byrow=FALSE)
    results_tmp$fL <<- fl
    
    # bh, bl
    start=end+1; end=start+k-1
    bh <- matrix(par[start:end],nrow=k, ncol=1, byrow=FALSE)
    results_tmp$bH <<- bh
    start=end+1; end=start+k-1
    bl <- matrix(par[start:end],nrow=k, ncol=1, byrow=FALSE)
    results_tmp$bL <<- bl
    
    # mu0h, mu0l
    start=end+1; end=start
    mu0h <- par[start:end]
    results_tmp$mu0H <<- mu0h
    start=end+1; end=start
    mu0l <- par[start:end]
    results_tmp$mu0L <<- mu0l
    
    # thetah, thetal
    start=end+1; end=start
    thetah <- par[start:end]
    results_tmp$thetaH <<- thetah
    start=end+1; end=start
    thetal <- par[start:end]
    results_tmp$thetaL <<- thetal
    
    # p
    start=end+1; end=start
    p <- par[start:end]
    results_tmp$p <<- p
    
    # End reading parameters
    
    if(stopifbound) {
      for(i in 1:length(results_tmp)) {
        if(length(intersect(results_tmp[[i]],c(bounds$lower_bound[i], bounds$upper_bound[i]))) >= 1) {
          cat("Parameter", names(results_tmp)[i], "achieved lower/upper bound. Process stopped.\n")
          cat(results_tmp[[i]],"\n")
          stopflag <- T
          break
        }
      }
    }
    
    res <- 0
    if(stopflag == FALSE) {
      if(mode == "observed" | mode == "combined") {
        res1 <- .Call("complik_gen", dat, dim(dat)[1], dim(dat)[2], 
                   ah, al,
                   f1h, f1l,
                   Qh, Ql,
                   bh, bl, 
                   fh, fl, 
                   mu0h, mu0l, 
                   thetah, thetal, 
                   p, 
                   N.c, N.nc, 
                   k, pinv.tol,
                   gomp)
        if(mode == "combined") {
          res2 <- .Call("complikGenNonGenetic", unobsdat, dim(unobsdat)[1], dim(unobsdat)[2], 
                        ah, al,
                        f1h, f1l,
                        Qh, Ql,
                        bh, bl, 
                        fh, fl, 
                        mu0h, mu0l, 
                        thetah, thetal, 
                        p,
                        k, pinv.tol, gomp)
          
          res <- res1 + res2
        } else {
          res <- res1
        }
      } else if(mode == "unobserved") {
        res <- .Call("complikGenNonGenetic", unobsdat, dim(unobsdat)[1], dim(unobsdat)[2],
                     ah, al,
                     f1h, f1l,
                     Qh, Ql,
                     bh, bl, 
                     fh, fl, 
                     mu0h, mu0l, 
                     thetah, thetal, 
                     p,
                     k, pinv.tol, gomp)
      }
      
      assign("results", results_tmp, envir=baseenv())
      iteration <<- iteration + 1
      if(verbose) {
        cat("L = ", res,"\n")
        cat("Iteration: ", iteration,  "\nResults:\n") 
        print(results_tmp)
      }
      
    }
    
    return(as.numeric(-1*res))
  }
  
  # Optimization:
  if(verbose) {
    cat("Lower bound:\n")
    print(bounds$lower_bound)
    cat("Upper bound:\n")
    print(bounds$upper_bound)
  }
  #tryCatch({ans <- nloptr(x0 = parameters, 
  #               eval_f = maxlik, opts = list("algorithm"=algorithm, 
  #                                            "xtol_rel"=1.0e-1, "maxeval"=maxeval),
  #               lb = bounds$lower_bound, ub = bounds$upper_bound)
  #          
  #         },  
  #         error=function(e) {if(verbose  == TRUE) {print(e)}}, 
  #         finally=NA)
  
  nloptr(x0 = parameters, 
         eval_f = maxlik, opts = list("algorithm"=algorithm, "ftol_rel"=ftol_rel, "maxeval"=maxeval),
         lb = bounds$lower_bound, ub = bounds$upper_bound)
  
  final_results <- get("results",envir=baseenv())
  
  # Check if any parameter achieved upper/lower limit and report it:
  limit <- FALSE
  for(i in 1:length(final_results)) {
    if(length(intersect(final_results[[i]],c(bounds$lower_bound[i], bounds$upper_bound[i]))) >= 1) {
      cat("Parameter", names(final_results)[i], "achieved lower/upper bound.\n")
      cat(final_results[[i]],"\n")
      limit <- TRUE
    }
  }
  final_results$limit <- limit
  #assign("results", final_results, envir=baseenv())
  class(final_results) <- "pobs.spm"
  
  invisible(final_results)
}


