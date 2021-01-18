###########################################################################
# Modelling the reduction in faecal egg count data (two-sample case) using Stan
###########################################################################

# main function -----------------------------------------------------------

fecr_stanSimple <-function(preFEC, postFEC, rawCounts = FALSE, preCF = 50, postCF = preCF, 
                    muPrior, deltaPrior, nsamples = 2000, nburnin=1000, thinning=1, nchain=2, ncore=1, adaptDelta=0.95,
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
  
  # set default values
  priors <- fecr_setPrior(muPrior = muPrior, deltaPrior = deltaPrior)
  
  # set update functions
  code<-simple_paired_stan(priors);model<-"Simple Bayesian model without zero-inflation for paired design"
  
  preN <- length(preFEC)
  postN <- length(postFEC)
  if(preN != postN){
    stop("post sample size different to pre sample size\n")
  }
  
  if (length(preCF)==1) preCF<-rep(preCF,preN)
  if (length(postCF)==1) postCF<-rep(postCF,postN)
  
  # create data list for stan use
  epg_data <- list(J=preN, ystarbraw = preFEC, ystararaw = postFEC, fpre = preCF, fpost = postCF)

  if (saveAll){ savePars <- NA} else { savePars <- c("mu","delta")}
  
  # whether or not to suppress progress information and errors
  if (length(setdiff(priors,fecr_setPrior()))==0){
    stanModel<-stanmodels$simple
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
  
  # generate samples according to different models

      meanEPG.untreated<-extract(samples,"mu")[[1]]
      meanEPG.treated<-extract(samples,"mu")[[1]]*extract(samples,"delta")$delta
      FECR<-1-extract(samples,"delta")[[1]]
      result<-cbind(FECR,meanEPG.untreated,meanEPG.treated)
 
  
  cat("Model: ", model,"\n","Number of Samples: ",nsamples, "\n","Warm-up samples: ",nburnin,"\n","Thinning: ",thinning,"\n","Number of Chains",nchain,"\n")
  summarys<-as.data.frame(printSummary(result))
  
  checkConvergence(samples)
  
  if (preN + postN < 20) warning(cat("your sample size is less than 10, consider using getPrior_mu() and getPrior_delta() to find a more informative prior for the true mean epg and reduction parameter."))
  
  return(invisible(list(stan.samples = samples, posterior.summary = summarys)))
}