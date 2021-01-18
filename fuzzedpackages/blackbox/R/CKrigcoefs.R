# This function estimates any of the the NA parameters, or returns c and d if all parameters are given
CKrigcoefs <- function(xy,
                     nuniquerows, ## as given by selectFn
                     covfnparamA=NA, # excluding smoohtness, but output covfnparam includes it.
                     lambdaA=NA,
                     minSmoothness=blackbox.getOption("minSmoothness"),
                     maxSmoothness=blackbox.getOption("maxSmoothness"),
                     miscOptions=blackbox.getOption("miscOptions"),
                     initCovFnParam=NULL,
                     verbosity=blackbox.getOption("verbosity"),
                     optimizers=blackbox.getOption("optimizers")
){
  nrowxy <- nrow(xy)
  ncolxy <- ncol(xy)
  # Test on length(covfnparamA) assumes that when scale params are known, smoothness is too 
  # (=> testing smoothness range is irrelevant). 
  # In standard usage either covfnparamA is zero-length (and lambda is unknown)
  #  => the estimation of cov params is performed
  # or all parameters (including lambda) are known and no optimization is performed.
  # It is also possible that the user sets a scale parameter, and and smoothness through min/maxSmoothness,
  # but lambda still needs to be estimated: this is a subcase of the optimize case.
  # (we might only test lambda to distinguish the cases)
  # The user is currently not allowed to set all correlation params through the CovFnParam vector: 
  # input CovFnParam should contain scale parameters only.
  # The user is also not allowed to provide lambda (through blackbox.getOption("lambdaEst")).
  varCorrpars <- (length(covfnparamA)< ncolxy ## shorter than scale params+smoothness
                  || any(is.na(covfnparamA))) 
  optimise <- (varCorrpars || is.na(lambdaA))
  method <- intersect(optimizers,c("bobyqa","L-BFGS-B","lbfgsb3c")) ## non-default methods
  if(length(method)==0L) method <- "NLOPT_LN_BOBYQA" ## default for this part of code.
  if(length(method)!=1L) stop("length(method)!=1L in CKrigcoefs") ## incompatible with switch below
  verbosityFromObjective <- as.integer((optimise && "L-BFGS-B" %in% method))*verbosity ## FR->FR should control in the C++ code using the eval counter
  verbosityFromoptimizer <- switch(method,
                                   "bobyqa"=round(2*as.integer(optimise)*verbosity), ## default=2
                                   "L-BFGS-B"=0, ## verbosity from optim() is awkward
                                   as.integer(optimise)*verbosity
  )
  if ("Migraine" %in% blackbox.getOption("usedBy")) {
    init_factor <- 2/3
  } else init_factor <- 2/30 ## better when rugged landscapes are possible  
  success <- newCSmooth(xy= t(xy), ## newCSmooth is an Rcpp export
             nrowxy=nrowxy,
             ncolxy=ncolxy,
             nuniquerows=nuniquerows,
             GCV=0, # 0 for GCV, 1 for match_fs2hat_pure_error
             optimiseBool=optimise,
             verbosity=verbosityFromObjective)
  if ( ! success ) {
    resu <- list()
    ## go directly  to deleteCSmooth
  } else if ( ! optimise ) {
    resu <- Krig_coef_Wrapper(covfnparamA,lambdaA)[c("c","d","CKrigidx")] # throwing away u and D
    # do not deleteCSmooth(), CSmooth pointer stored in CKrigptrTable
  } else { # optimize
    if (varCorrpars) { ## need to estim scale params and/or smoohtness (not sure that it can estimate only some of them)
      if ( is.null(initCovFnParam) ) initCovFnParam <- rep(0, ncolxy)
      xx <- xy[,-ncolxy,drop=FALSE]
      KgLow <- apply(xx,2,min)
      KgUp <- apply(xx,2,max)
      maxrange <- KgUp-KgLow
      GCVlowerFactor <- 20 # or 20 ? cf comments on Migraine default in Migraine code
      GCVupperFactor <- 5
      lower <- maxrange/(GCVlowerFactor*nuniquerows)
      upper <- maxrange*GCVupperFactor
      init_factor <- rep(init_factor,length(lower))
      minSmoothness <- min(maxSmoothness,max(minSmoothness,1.001)); # >1, cf def Matern
      varSmoothness <- (maxSmoothness-minSmoothness>1e-6)
      if (varSmoothness) { # Variable Smoothness case
        maxrange <- c(maxrange,maxSmoothness-minSmoothness)
        lower <- c(lower,minSmoothness)
        upper <- c(upper,maxSmoothness)
        init_factor <- c(init_factor,2/3)
        fixedSmoothness <- numeric(0) ## not c() which is NULL...
      } else {
        initCovFnParam <- initCovFnParam[seq_len(length(lower))]
        fixedSmoothness <- maxSmoothness
      }
      for (it in seq_len(length(lower))) {
        if (initCovFnParam[it]<lower[it] ## includes default case
            || initCovFnParam[it]>upper[it]) {
          initCovFnParam[it]=lower[it]+maxrange[it]*init_factor[it]
        }
      }
      if ("L-BFGS-B" %in% method) {
        control <- list(parscale=(upper-lower), trace=verbosityFromoptimizer)
        optr <- optim(par=initCovFnParam,fn=GCV_lamVar_covFix_Wrapper,method="L-BFGS-B",
                      lower=lower,upper=upper,control=control,
                      fixedSmoothness=fixedSmoothness,returnFnvalue=TRUE)
        solution <- optr$par
      } else if ("lbfgsb3c" %in% method){ ## very slow
        control <- list(trace=verbosityFromoptimizer)
        if ( ! requireNamespace("lbfgsb3c",quietly=TRUE) ) {
          stop("Package lbfgsb3c not installed.")
        }
        optr <- lbfgsb3c::lbfgsb3c(par=initCovFnParam,fn=GCV_lamVar_covFix_Wrapper,lower=lower,upper=upper,
                        control=control,
                        fixedSmoothness=fixedSmoothness,returnFnvalue=TRUE)
        solution <- optr$par # change for lbfgsb3c
      } else if ("bobyqa" %in% method){ ## marginally better than optim ?
        if ( ! requireNamespace("minqa",quietly=TRUE) ) {
          stop("Package minqa not installed.")
        }
        control <- list(rhobeg=min(abs(upper-lower))/20,iprint=verbosityFromoptimizer)
        control$rhoend <- max(1,control$rhobeg)/1e6
        optr <- minqa::bobyqa(par=initCovFnParam,fn=GCV_lamVar_covFix_Wrapper,lower=lower,upper=upper,
                       control=control,
                       fixedSmoothness=fixedSmoothness,returnFnvalue=TRUE)
        solution <- optr$par
      } else if ("NLOPT_LN_BOBYQA" %in% method){ ## may be the fastest
        optr <- nloptr(x0=unlist(initCovFnParam),eval_f=GCV_lamVar_covFix_Wrapper,lb=lower,ub=upper,
                       opts=list(algorithm="NLOPT_LN_BOBYQA",xtol_rel=1.0e-4,maxeval=-1,print_level=verbosityFromoptimizer),
                       fixedSmoothness=fixedSmoothness,returnFnvalue=TRUE)
        solution <- optr$solution
      } else stop("Unknown 'method' in CKrigcoefs()")
      solution <- c(solution,fixedSmoothness) # required when minSmoothness=maxSmoothness
      resu <- list(covfnparam=solution, # full including smoothness whether fixed or estimated
                   fnEvalCount=getFnEvalCount(),method=method)
    } else { ## no estim of scale nor smoohtness but we will estimate lambda
      resu <- list(covfnparam=c(covfnparamA,maxSmoothness)) # also including smoothness
    }
    if (is.na(lambdaA)) resu$lambda <- GCV_lamVar_covFix_Wrapper(resu$covfnparam,fixedSmoothness=numeric(0),returnFnvalue=FALSE)
    deleteCSmooth() ## this code conditional on optimisation
  }
  return(resu) # is a list with elements depending on call:
  # if optimise covfnparam -> returns covfnparam, lambda, fnEvalCount, method
  # else if optimise lambda -> returns (input covfnparamA), lambda
  # else returns c("c","d","CKrigidx")
} ## end def CKrigcoefs

# R-style implementation of default values
R_GCV_lamVar_covFix <- function(a,fixedSmoothness=numeric(0L),returnFnvalue=TRUE) {
  ## call to R-to-C wrapper with explicit values
  GCV_lamVar_covFix_Wrapper(a=a,fixedSmoothness=fixedSmoothness,returnFnvalue=returnFnvalue)
}
