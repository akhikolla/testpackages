
setMethod("mle", 
          signature = signature(object = "covAll"),
          definition =
              function(object,
                       y, X, F = NULL, beta = NULL,
                       parCovIni = coef(object),
                       parCovLower = coefLower(object),
                       parCovUpper = coefUpper(object),
                       noise = TRUE,
                       varNoiseIni = var(y) / 10,
                       varNoiseLower = 0,
                       varNoiseUpper = Inf,
                       compGrad = hasGrad(object),
                       doOptim = TRUE,
                       optimFun = c("nloptr::nloptr", "stats::optim"),
                       optimMethod = ifelse(compGrad, "NLOPT_LD_LBFGS",
                           "NLOPT_LN_COBYLA"),
                       optimCode = NULL,
                       multistart = 1,
                       ## method,
                       ## control = list(fnscale = -1, trace = 3, REPORT = 1),
                       parTrack = FALSE,
                       trace = 0,
                       checkNames = TRUE,
                       ...) {
                  
                  if (compGrad && !hasGrad(object)) {
                      stop("when 'compGrad' is given and is TRUE, 'cov' object ",
                           "must compute the gradient")
                  }
                
                  optimFun <- match.arg(optimFun)
                  if ((optimFun == "stats::optim") && missing(optimMethod)) {
                      optimMethod <- "L-BFGS-B"
                  }
                  
                  Ldots <- list(...)

                  ## cat("XXX Ongoing development for optim\n")
                  ## print(names(Ldots))

                  if ("method" %in% names(Ldots)) {
                      warning("The formal argument 'method' should no longer",
                              " be used. Use `optimFun' and 'optimMethod' to",
                              " specify the optimisation method.")
                  }
            
                  if (checkNames) {
                      X <- checkX(object, X = X)
                  }
                  
                  co <- coef(object)
                  lpar <- length(co) 
                  parNames <- names(co)
                  
                  # if (trace) {
                  #     cat("Initial values, lower and upper bounds\n")
                  #     mat <- rbind(parIni, parLower, parUpper)
                  #     rownames(mat) <- c("Initial", "lower", "upper")
                  #     print(mat)
                  # }
                  
                  if (!is.null(F)) {
                      pF <- NCOL(F)
                  } else pF <- 0L
                  
                  if (!is.null(beta)) {
                      if (pF != length(beta)) stop("'beta' and 'F' mismatch")
                      betaGiven <- TRUE
                      thisy <- y - F %*% beta
                      thisF <- NULL
                  } else {
                      thisy <- y
                      thisF <- F
                      betaGiven <- FALSE
                  }
                  
                  ##===========================================================
                  ## the loglik function to be maximised. Args, 'y',
                  ## 'X', 'F', ## 'gradEnv' comes by lexical scoping
                  ## andwill remain attached to the funs on output.
                  ## ========================================================== 
                  par.all <- rep(NA, lpar)
                  names(par.all) <- parNames
                  
                  ## Where to store the parameters and gradient
                  gradEnv <- new.env()
                  
                  if (parTrack) {
                      
                      gradEnv$parTracked <- numeric(0)
                      
                      logLikFun <- function(par) {
                          gradEnv$parTracked  <- rbind(gradEnv$parTracked, par)
                          .logLikFun0(par, object, y = thisy, X, F = thisF,
                                      compGrad = compGrad, noise = noise,
                                      gradEnv = gradEnv, trace = trace) 
                      }
                      
                  } else {
                      
                      logLikFun <- function(par) {
                          ll <- .logLikFun0(par, object, y = thisy, X, F = thisF,
                                      compGrad = compGrad, noise = noise,
                                      gradEnv = gradEnv, trace = trace) 
                      }
                      
                  }
                  
                  if (compGrad) {
                      
                      thisLogLikGrad  <- function(par) {
                          stored.par <- get("par", envir = gradEnv)
                          if ( !identical(par, stored.par) ) {
                              stop("'par' arg and 'par' stored in 'gradEnv' ",
                                   "are not identical")
                          }
                          stored.logLik.derivative <-
                              get("LLgrad", envir = gradEnv)
                          logLik.derivative <- stored.logLik.derivative 
                      }
                      
                      negLogLikFun <- function(par) {
                          
                          ll <- logLikFun(par)
                          
                          if (!is.na(ll)) {
                              return(list("objective" = -ll,
                                          "gradient" = -attr(ll, "gradient")))
                          } else {
                              return(list("objective" = NaN,
                                          "gradient" = NaN))
                          }
                          
                      }

                  } else {
                      
                      thisLogLikGrad <- NULL
                      
                      negLogLikFun <- function(par) {
  
                          ll <- logLikFun(par)
                          
                          if (!is.na(ll)) return(-ll)
                          else return(NaN)
                          
                      }
                      
                  }
                  
                  ## ===========================================================
                  ## Unless 'doOptim' is turned to FALSE, compute the
                  ## kernel part of the coefficients
                  ## ===========================================================
                  
                  if (doOptim) {
                      
                      parCovVec <- list(parCovLower, parCovUpper)
                      parCovVecNames <- c("parCovLower", "parCovUpper")
                      
                      for (i in 1L:2L) {
                          if (length(parCovVec[[i]]) != object@parN)
                              stop(parCovVecNames[i], " should be of length ",
                                   object@parN)
                      }
                      
                      needSimul <- FALSE
                      if (!missing(parCovIni)) {
                          if (is.vector(parCovIni)) {
                              parCovIni <- t(as.matrix(parCovIni))
                          }
                          if (ncol(parCovIni) != length(coef(object))) {
                              stop("Invalid number of initial covariance",
                                   " parameters in 'covParIni'")
                          }
                          if ((nrow(parCovIni) != multistart)) {
                              warning("'parCovIni' is ignored, since",
                                      " nrow(parCovIni) is not equal to",
                                      " 'multistart'. Using 'simulPar'")
                              needSimul <- TRUE
                          } else {
                              myFun <- function(x) {
                                  all((x <= parCovUpper) & (x >= parCovLower))
                              }
                              if (!all(apply(parCovIni, 1, myFun))) 
                                  stop("'parCovIni' should be in [parCovLower, ",
                                       "parCovUpper]")
                          }
                      } else {
                          if (multistart > 1) needSimul <- TRUE
                          parCovIni <- t(as.matrix(parCovIni))
                      }
                      
                      if (needSimul) {
                          parCovIni <- simulPar(object, nsim = multistart)
                      }
                      
                      parIni <- parCovIni
                      parLower <- parCovLower
                      parUpper <- parCovUpper
                      lparNN <- lpar


                      ## When 'varNoiseIni' is not given by the user but set
                      ## to its default value, we can make it fall between the
                      ## bounds.  
                      if (noise) {

                          if (varNoiseLower > varNoiseUpper) {
                              stop("'varNoiseLower' must be <= 'varNoiseUpper'") 
                          }
                          
                          if (missing(varNoiseIni)) {
                              if (varNoiseIni < varNoiseLower)
                                  varNoiseIni <- varNoiseLower
                              if (varNoiseIni > varNoiseUpper)
                                  varNoiseIni <- varNoiseUpper
                              
                          } else {
                              if ((varNoiseIni < varNoiseLower) ||
                                  (varNoiseIni > varNoiseUpper)) {
                                  stop("'varNoiseIni' must be >=",
                                       " 'varNoiseLower' and <=",
                                       " 'varNoiseUpper'.")
                              }
                          }
                          parNames <- c(parNames, "varNoise")
                          parIni <- cbind(parIni, 
                                    varNoise = rep(varNoiseIni,
                                        length.out = multistart))
                          parLower <- c(parLower, varNoiseLower)
                          parUpper <- c(parUpper, varNoiseUpper)
                          lpar <- lpar + 1
                      } else {
                          if (!missing(varNoiseIni) || !missing(varNoiseLower)
                              || !missing(varNoiseUpper)) {
                              warning("Since 'noise' is FALSE, the values of ",
                                      "the 'varNoise*' arguments are ignored")
                          }
                      }

                      if (is.null(optimCode)) {

                          if (trace) {
                              cat(sprintf("optimFun :    %s\n", optimFun))
                              cat(sprintf("optimMethod : %s\n\n", optimMethod))
                          }
                          
                          if (optimFun == "stats::optim") {
                              
                              if ("method" %in% names(Ldots)) {
                                  warning("Optional argument 'method' given in ",
                                          "'...' is ignored: ", 
                                          "'optimMethod' is used instead.")
                              }  
                              Ldots$control$fnscale <- -1
                              fitList <- list()
                              for (i in 1:multistart){
                                  args <- list(par = parIni[i, ], fn = logLikFun,
                                               method = optimMethod)
                                  if (optimMethod %in% c("L-BFGS-B", "Brent")) {
                                      args <- c(args, list(lower = parLower, upper = parUpper))
                                  }
                                  ## if (compGrad & is.element(optimMethod, c("BFGS", "L-BFGS-B"))){
                                  if (compGrad) {
                                      args <- c(args, list(gr = thisLogLikGrad))
                                  }
                                  args <- c(args, Ldots)
                                  fitList[[i]] <- try(do.call(stats::optim, args = args))
                              }
                              
                              ## if (requireNamespace("foreach", quietly = TRUE)){
                              ##   fitList <- foreach::"%dopar%"(foreach::foreach(i=1:multistart,
                              ##                                      .errorhandling='remove'), {
                              ##   Ldots$control$fnscale <- -1
                              ##   args <- list(par = parIni[i, ], fn = logLikFun, method = optimMethod)
                              ##   if (optimMethod == "L-BFGS-B") {
                              ##     args <- c(args, list(lower = parLower, upper = parUpper))
                              ##   }
                              ##   if (compGrad & is.element(optimMethod, c("BFGS", "L-BFGS-B"))){
                              ##     args <- c(args, list(gr = thisLogLikGrad))
                              ##   }
                              ##   args <- c(args, Ldots, list(envir = gradEnv))
                              ##   try(do.call(stats::optim, args = args))
                              ##   })
                              ## } # end foreach
                              
                              nlist <- length(fitList)
                              ## get the best result
                              if (nlist == 0) stop("All fit trials have failed")
                              ## print(fitList)
                              optValueVec <- sapply(fitList, function(x) x$value)
                              bestIndex <- which.max(optValueVec)
                              report <-
                                  list(parIni = parIni,
                                       par = t(sapply(fitList, function(x) x$par)),         
                                       logLik = - sapply(fitList, function(x) x$value),
                                       counts  = sapply(fitList, function(x) x$counts[1]),
                                       convergence = sapply(fitList,
                                           function(x) x$convergence == 0))
                              colnames(report$par) <- colnames(report$parIni)
                              if (trace > 0) {
                                  cat("\nOptimisation report:\n\n")
                                  print(as.data.frame(report))
                              }
                              
                              opt <- fitList[[bestIndex]]
                              
                          } else if (optimFun == "nloptr::nloptr") {
                              
                              ## if (is.element(optimMethod, c("NLOPT_LD_LBFGS", "NLOPT_LN_COBYLA"))) {
                              
                              if (!compGrad & (optimMethod == "NLOPT_LD_LBFGS")) {
                                  stop("when 'optimFun' is \"nloptr::nloptr\"",
                                       " and 'optimMethod' is \"NLOPT_LD_LBFGS\"",
                                       " 'compGrad' must be TRUE and the",
                                       " gradients must be provided")
                              }
                              
                              ## Note: A warning is done by nloptr for
                              ## NLOPT_LN_COBYLA if the gradient is supplied
                              
                              if (("opts" %in% names(Ldots)) &&
                                  ("algorithm" %in% names(Ldots$opts))) {
                                  warning("Optional argument 'opts$algorithm' given in '...' ",
                                          "is ignored: 'optimMethod' is used instead.")
                              }  
                              Ldots$opts$algorithm <- optimMethod
                              
                              if (requireNamespace("foreach", quietly = TRUE)){
                                  fitList <-
                                      foreach::"%dopar%"(foreach::foreach(
                                          i = 1:multistart, 
                                          .errorhandling='remove'), {
                                              args <- list(x0 = parIni[i, ], eval_f = negLogLikFun, 
                                                           lb = parLower, ub = parUpper)
                                              args <- c(args, Ldots)
                                              do.call(nloptr::nloptr, args = args)
                                              ## try removed, due to conflicts with
                                              ## '.errorhandling' (thanks Mickael Binois!)
                                          })
                              }
                              
                              nlist <- length(fitList)
                              ## get the best result
                              if (nlist==0) stop("All fit trials have failed")
                              optValueVec <- sapply(fitList, function(x) x$objective)
                              bestIndex <- which.min(optValueVec)
                              report <-
                                  list(parIni = parIni,
                                       par = t(sapply(fitList, function(x) x$solution)),         
                                       logLik = - sapply(fitList, function(x) x$objective),
                                       nIter  = sapply(fitList, function(x) x$iterations),
                                       convergence = sapply(fitList, function(x) x$status %in% c(1, 3, 4)))
                              colnames(report$par) <- colnames(report$parIni)
                              if (trace > 0) {
                                  cat("\nOptimisation report:\n\n")
                                  print(as.data.frame(report))
                              }
                              opt <- fitList[[bestIndex]]  # end multistart
                              
                              ## } else {
                              ## 
                              ##     warning("When 'optimFun' is \"nloptr::nloptr\" and",
                              ##             " and 'optimMethod' is not",
                              ##             " \"NLOPT_LD_LBFGS\", you must provide an",
                              ##             " 'opts' argument with an \"algorithm\"",
                              ##             " element. See the documentation of the",
                              ##             " 'nloptr' package.")
                              ##     
                              ##     opt <- try(nloptr::nloptr(x0 = parIni,
                              ##                               eval_f = negLogLikFun, ...))
                              ## }
                              
                              if (!inherits(opt, "try-error")) {
                                  opt$par <- opt$solution
                                  opt$value <- -opt$objective
                                  ## see 'opt$status' to get more info
                                  opt$convergence <- opt$status %in% c(1, 3, 4)
                              }
                          }
                          
                      } else {  ## case !is.null(optimCode)

                          if (trace) {
                              cat("Using the provided value of 'optimCode'\n")
                              cat("=======================================\n")
                          }
                          
                          report <- NULL
                          eval(parse(text = optimCode))
                          
                      } 
                      
                      if ( !inherits(opt, "try-error") ) {  
                          coef(object) <- opt$par[1:lparNN]
                          if (noise) {
                              varNoise <- opt$par[lpar]
                          } else varNoise <- NULL
                          if (opt$convergence) {
                              warning("optimisation did not converge\n")
                          }
                          
                      } else {
                          stop("error in 'optim'\n")
                          ## return(NULL)
                      }
                      
                      ## coef(object) <- coef.kernel[1L:lparNN]
                      
                  } else {
                      opt <- NULL
                      coef.kernel <- NULL
                  } ## end if (doOptim) ... else
                  
                  ## compute 'beta' if needed
                  
                  if (!betaGiven) {
                      if (pF) {
                          trendRes <- gls(object = object, y = y, X = X, F = F,
                                          varNoise = varNoise, 
                                          checkNames = FALSE) ## already checked
                          ##   coef.trend <- glsRes$betaHat
                      } else {
                          trendRes <- list(betaHat = NULL)
                      }
                  } else {
                      ##     trendRes <- list(betaHat = beta) 
                      trendRes <- gls(object = object, y = y, X = X, F = F,
                                      varNoise = varNoise,
                                      beta = beta,
                                      checkNames = FALSE) ## already checked
                  }
                  
                  ## ==========================================================
                  ## remove the 'parTrack' augmentation in the
                  ## function for output XXX keep the 'compGrad' ansd
                  ## 'gradEnv' values????
                  ## ==========================================================
                  
                  if (parTrack) {
                      logLikFun <- function(par) {
                          .logLikFun0(par, object, y = thisy, X, F = thisF,
                                      compGrad = TRUE,
                                      noise = noise,
                                      gradEnv = gradEnv,
                                      trace = FALSE) 
                      }
                      parTracked <- gradEnv$parTracked
                      rownames(parTracked) <- paste("it. ", 0:(NROW(parTracked)-1L))
                      
                      colnames(parTracked) <- parNames
                      gradEnv$parTracked <- NULL
                  } else {
                      parTracked <- NULL
                  }
                  
                  list(logLik = opt$value,
                       parIni = parIni,
                       opt = opt,
                       cov = object,
                       noise = noise,
                       varNoise = varNoise,
                       trendKnown = betaGiven,
                       trendRes = trendRes,
                       parTracked = parTracked,
                       logLikFun = logLikFun, 
                       report = report)
                  
              }
          )
