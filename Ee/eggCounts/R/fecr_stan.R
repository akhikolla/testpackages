###########################################################################
# Modelling the reduction in faecal egg count data (two-sample case) using Stan
###########################################################################

# check if the arguments in fecr_stan is sensible -------------------------

checkpars <- function(nburnin, nsamples, thinning, nchain, ncore, rawCounts, adaptDelta, verbose)
{
  if ((nburnin < 1)|(ceiling(nburnin)!=floor(nburnin)))
    stop("'nburnin' should be a positive integer", call.=FALSE)
  if ((nsamples < 1)|(ceiling(nsamples)!=floor(nsamples)))
    stop("'nsamples' should be a positive integer", call.=FALSE)
  if ((nsamples <= nburnin))
    stop("The total number of samples must be greater than the number of burn-in samples", call.=FALSE)
  if ((thinning < 1)|(ceiling(thinning)!=floor(thinning)))
    stop("'thinning' should be a positive integer", call.=FALSE)
  if ((nchain < 1)|(ceiling(nchain)!=floor(nchain)))
    stop("The number of chains must be a positive integer", call.=FALSE)
  if ((ncore < 1)|(ceiling(ncore)!=floor(ncore)))
    stop("The number of cores must be a positive integer", call.=FALSE)
  if (!is.logical(rawCounts))
    stop("The rawCounts argument must be a logical", call.=FALSE)
  if (!is.logical(verbose))
    stop("The verbose argument must be a logical", call.=FALSE)
  if ((adaptDelta <= 0)||(adaptDelta >= 1))
    stop("adaptDelta must be between 0 and 1", call.=FALSE)
  if (nburnin < 500)      warning("'nburnin' seems low.\n")
  if (nsamples < 1000)     warning("'nsamples' seems low.\n")
  invisible()
}

checkData <- function(preFEC, postFEC, rawCounts, preCF, postCF){
  # number of faecal samples
  preN <- length(preFEC)
  postN <- length(postFEC)
  
  # check correction factors
#  checkCF <- function(CF){(CF < 0)|(ceiling(CF)!=floor(CF))}
#  if (any(sapply(preCF,checkCF)))    stop("correction factor(s) should be a positive integer", call.=FALSE)
#  if (any(sapply(postCF,checkCF)))     stop("correction factor(s) should be a positive integer", call.=FALSE)
  if(length(preCF)>1 && length(preCF)!=preN) stop("Lengths of the vectors preCF and preFEC do not match\n")
  if(length(postCF)>1 && length(postCF)!=postN) stop("Lengths of the vectors postCF and postFEC do not match\n")
  
  # raw counts or EpGs?
  preDilution <- preCF; postDilution <- postCF
  if(rawCounts){
    preDilution <- postDilution <- 1
  }
  
  # divide data by correction factor
  if(any(sapply(preFEC, function(x) min(abs(c(x %% preDilution, x %% preDilution - preDilution)))) > 0.001)) stop(paste(c("Correction factor preCF does not match the given pre-treatment faecal egg counts. Non-integer detected.")))
  preFEC <- preFEC/preDilution
  if(any(sapply(postFEC, function(x) min(abs(c(x %% postDilution, x %% postDilution - postDilution)))) > 0.001))  stop(paste(c("Correction factor postCF does not match the given post-treatment faecal egg counts. Non-integer detected.")))
  postFEC <- postFEC/postDilution
  
  if (!rawCounts){
    if ( mean( preFEC) < mean( postFEC))
      warning("mean of pre-treatment is smaller of post-treatment. Results may be unreliable if default priors are used.\n")
    if ( median( preFEC) < median( postFEC))
      warning("median of pre-treatment is smaller of post-treatment. Results may be unreliable if default priors are used.\n")
  } else {
    if ( mean( preFEC*preCF) < mean( postFEC*postCF))
      warning("mean of pre-treatment is smaller of post-treatment. Results may be unreliable if default priors are used.\n")
    if ( median( preFEC*preCF) < median( postFEC*postCF))
      warning("median of pre-treatment is smaller of post-treatment. Results may be unreliable if default priors are used.\n")
  }
  # check counts
  if (preN==postN){if(sum( (preFEC-postFEC)^2) <= .Machine$double.eps)
    warning("the pre-treatment and post-treatment counts are identical, check the 'preFEC' and 'postFEC' arguments.\n")}
}

checkConvergence <- function(samples){
  rhats <- data.frame(summary(samples)$summary)
  rhats <- rhats[-nrow(rhats),] # delete the loglikelihood row
  rhatPar <- rownames(rhats)[which(rhats$Rhat > 1.1)]
  
  if (length(rhatPar) > 0){
    warning("there are evidence of non-convergence for parameter: ", rhatPar, " since their potential scale reduction factors (Brooks and Gelman, 1998) are greater than 1.1.\n")
  } else {
    cat("\nNOTE: there is no evidence of non-convergence since all parameters have potential scale reduction factors (Brooks and Gelman, 1998) less than 1.1.\n")
  }
}

checkDivergence <- function(samples, adaptDelta){
  sp <- get_sampler_params(samples, inc_warmup = FALSE)
  n_d <- sum(sapply(sp, FUN = function(x) {
    if ("divergent__" %in% colnames(x)) return(sum(x[, "divergent__"])) else return(0)
  }))
  if (n_d > 0) {
    warning("there were ", n_d, " divergent transitions after warmup, the joint posterior distribution is not sufficiently explored.\n", 
            " Re-run the model with adaptDelta > ", adaptDelta, " may help, or results can be unreliable.\n", call. = FALSE)}
}


# main function -----------------------------------------------------------



fecr_stan<-function(preFEC, postFEC, rawCounts = FALSE, preCF = 50, postCF = preCF,
                    paired = TRUE, indEfficacy = TRUE, zeroInflation = FALSE, 
                    muPrior, kappaPrior, deltaPrior, phiPrior, deltakappaPrior,
                    nsamples = 2000, nburnin=1000, thinning=1, nchain=2, ncore=1, adaptDelta=0.95,
                    saveAll = FALSE, verbose=FALSE){
# checks from FECR_PoGa.R -------------------------------------------------
#   if (sys.parent() == 0) env <- asNamespace("eggCounts") else env <- parent.frame()
#   assign(".verboselevel", verbose*.verboselevel, envir = env)
  
  checkData(preFEC, postFEC, rawCounts, preCF, postCF)
  preN <- length(preFEC)
  postN <- length(postFEC)
  if (length(preCF)==1) preCF<-rep(preCF,preN)
  if (length(postCF)==1) postCF<-rep(postCF,postN)
  
  preDilution <- preCF; postDilution <- postCF
  if(rawCounts){
    preDilution <- postDilution <- 1
  }
  # divide data by correction factor
  preFEC <- preFEC/preDilution
  postFEC <- postFEC/postDilution
  
  # check function arguments
  checkpars(nburnin, nsamples, thinning, nchain, ncore, rawCounts, adaptDelta, verbose)
  if (!is.logical(paired))
    stop("the paired argument must be a logical", call.=FALSE)
  if (!is.logical(zeroInflation))
    stop("the zeroInflation argument must be a logical", call.=FALSE)
  
  if (indEfficacy == TRUE){
    if (paired == FALSE | zeroInflation == TRUE){
      stop(paste(c("the individual efficacy model only works in combination with 'paired = TRUE, zeroInflation = FALSE'.\n")))
    }
  }
  
  # set default values
  priors <- fecr_setPrior(muPrior = muPrior, kappaPrior = kappaPrior, deltaPrior = deltaPrior,
                          phiPrior = phiPrior, deltakappaPrior = deltakappaPrior)
  
  # set update functions
  if(paired & zeroInflation){code<-ZI_paired_stan(priors);model<-"Zero-inflated Bayesian model for paired design"}
  if(paired & !zeroInflation){
    if(indEfficacy){
      code<-indeff_stan(priors);model<-"Bayesian model without zero-inflation for paired design allowing individual efficacy"
    } else {
      code<-paired_stan(priors);model<-"Bayesian model without zero-inflation for paired design"
      }
    }
  if(!paired & zeroInflation){code<-ZI_unpaired_stan(priors);model<-"Zero-inflated Bayesian model for unpaired design"}
  if(!paired & !zeroInflation){code<-unpaired_stan(priors);model<-"Bayesian model without zero-inflation for unpaired design"}

  preN <- length(preFEC)
  postN <- length(postFEC)
  if(paired && preN != postN){
    stop("post sample size different to pre sample size\n")
  }
  
  if (length(preCF)==1) preCF<-rep(preCF,preN)
  if (length(postCF)==1) postCF<-rep(postCF,postN)
  
  # create data list for stan use
  epg_data <-if(paired){
    list(J=preN, ystarbraw = as.integer(preFEC), ystararaw = as.integer(postFEC), fpre = preCF, fpost = postCF)
  } else {
      list(Jb=preN, Ja=postN, ystarbraw = as.integer(preFEC), ystararaw = as.integer(postFEC), fpre = preCF, fpost = postCF)}
  
  if (saveAll){
    savePars <- NA
  } else {
    if (indEfficacy) {
      savePars <- c("mu","delta_mu","delta_shape","kappa")
    } else {
      if (zeroInflation){
        savePars <- c("mu","delta","kappa","phi")
      } else {
        savePars <- c("mu","delta","kappa")
      }
    }
  }
  
  # whether or not to suppress progress information and errors
  if (length(setdiff(priors,fecr_setPrior()))==0){
    if(paired & zeroInflation){stanModel<-stanmodels$zipaired}
    if(paired & !zeroInflation){
      if (indEfficacy){
        stanModel<-stanmodels$indefficacy
      } else {
        stanModel<-stanmodels$paired
      }
    }
    if(!paired & zeroInflation){stanModel<-stanmodels$ziunpaired}
    if(!paired & !zeroInflation){stanModel<-stanmodels$unpaired}
  } else {
  stanModel <- stan_model(model_name=paste(model),model_code=code)
  }
  
  if (verbose){
    samples <- sampling(stanModel, data = epg_data, pars = savePars, iter = nsamples, warmup=nburnin, chains=nchain,
                        thin=thinning,control = list(adapt_delta = adaptDelta),cores=ncore)
    
   } else {
  samples <- suppressMessages(
    suppressWarnings(
      sampling(stanModel, data = epg_data, pars = savePars, iter = nsamples,
               warmup = nburnin, chains = nchain, thin = thinning,
               control = list(adapt_delta = adaptDelta),cores = ncore, refresh = nsamples/4)))}
  
 checkDivergence(samples, adaptDelta)
  
  # tree depth effect the efficiency of sampler, not validity
  # n_m <- sum(sapply(sp, FUN = function(x) {
  #   if ("treedepth__" %in% colnames(x)) return(sum(x[, "treedepth__"] >= 
  #                                                    maxTreedepth)) else return(0)
  # }))
  # if (n_m > 0) {
  #   warning("There were ", n_m, " transitions after warmup that exceeded the maximum treedepth.", 
  #           " Increasing maxTreedepth above ", maxTreedepth, " may help, or results may be unreliable.", call. = FALSE)}
  
  # generate samples according to different models
  if(paired & !zeroInflation){
    if (indEfficacy){
      meanEPG.untreated <- extract(samples,"mu")[[1]]
      deltaMeansSample <- extract(samples,"delta_mu")[[1]]
      deltaShapeSample <- extract(samples,"delta_shape")[[1]]
      deltaSample <- qgamma(0.5, shape = deltaShapeSample, scale = deltaMeansSample/deltaShapeSample)
      FECR <- 1 - deltaSample
      meanEPG.treated <-extract(samples,"mu")[[1]] * deltaSample
      result<-cbind(FECR,meanEPG.untreated,meanEPG.treated)
    } else {
      meanEPG.untreated<-extract(samples,"mu")[[1]]
      meanEPG.treated<-extract(samples,"mu")[[1]]*extract(samples,"delta")$delta
      FECR<-1-extract(samples,"delta")[[1]]
      result<-cbind(FECR,meanEPG.untreated,meanEPG.treated)
    }
  }
  if(!paired & !zeroInflation){
      meanEPG.untreated<-extract(samples,"mu")[[1]]
      meanEPG.treated<-extract(samples,"mu")[[1]]*extract(samples,"delta")$delta
      FECR<-1-extract(samples,"delta")[[1]]
      result<-cbind(FECR,meanEPG.untreated,meanEPG.treated)
  }
         
  if(paired & zeroInflation){ # this prints out the group true epg including the zero-inflated components
           meanEPG.untreated<-extract(samples,"mu")[[1]]*(1-extract(samples,"phi")$phi)
           meanEPG.treated<-extract(samples,"mu")[[1]]*extract(samples,"delta")$delta*(1-extract(samples,"phi")$phi)
           FECR<-1-extract(samples,"delta")$delta
           result<-cbind(FECR,meanEPG.untreated,meanEPG.treated)
         }
   if(!paired & zeroInflation){ # this prints out the group true epg including the zero-inflated components
           meanEPG.untreated<-extract(samples,"mu")[[1]]*(1-extract(samples,"phi")$phi)
           meanEPG.treated<-extract(samples,"mu")[[1]]*extract(samples,"delta")$delta*(1-extract(samples,"phi")$phi)
           FECR<-1-extract(samples,"delta")$delta
           result<-cbind(FECR,meanEPG.untreated,meanEPG.treated)
   }


  cat("Model: ", model,"\n","Number of Samples: ",nsamples, "\n","Warm-up samples: ",nburnin,"\n","Thinning: ",thinning,"\n","Number of Chains",nchain,"\n")
  summarys<-as.data.frame(printSummary(result))

  checkConvergence(samples)
  
  if (preN + postN < 20) warning(cat("your sample size is less than 10, consider using getPrior_mu() and getPrior_delta() to find a more informative prior for the true mean epg and reduction parameter and/or using fecr_stanSimple() for a simplified Bayesian model for paired design."))
  
  return(invisible(list(stan.samples = samples, posterior.summary = summarys)))
}