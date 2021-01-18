###########################################################################
# Modelling faecal egg count data (one-sample case) using Stan
###########################################################################

fec_stan<-function(fec, rawCounts=FALSE, CF=50,
                    zeroInflation=TRUE, muPrior, kappaPrior, phiPrior,
                    nsamples=2000, nburnin=1000, thinning=1, nchain=2, ncore=1, adaptDelta=0.95, 
                    saveAll = FALSE, verbose=FALSE){
  # checks from FECR_PoGa.R -------------------------------------------------
  #   if (sys.parent() == 0) env <- asNamespace("eggCounts") else env <- parent.frame()
  #   assign(".verboselevel", verbose*.verboselevel, envir = env)
  
  # number of faecal samples
  n <- length(fec)
  # check correction factors
  # checkCF <- function(CF){(CF < 0)|(ceiling(CF)!=floor(CF))}
  # if (any(sapply(CF,checkCF)))     stop("correction factor(s) should be a positive integer", call.=FALSE)
  if(length(CF)>1 && length(CF)!=n) stop("lengths of the vectors for FEC and correction factors do not match\n")
  
  # raw counts or EpGs?
  dilution <- CF;
  if(rawCounts){
    dilution <- 1
  }
  
  # check iteration parameters
  checkpars(nburnin, nsamples, thinning, nchain, ncore, rawCounts, adaptDelta, verbose)
  if (!is.logical(zeroInflation))
    stop("the zeroInflation argument must be a logical", call.=FALSE)
  # divide data by correction factor
  if(any(sapply(fec, function(x) min(abs(c(x %% dilution, x %% dilution - dilution)))) > 0.001)) stop(paste(c("correction factor preCF does not match the given pre-treatment faecal egg counts. Non-integer detected.")))
  fec <- fec/dilution 
  
  # check model and set default values
  priors <- fec_setPrior(muPrior=muPrior, kappaPrior=kappaPrior, phiPrior=phiPrior)
  
  # set stan code and model name
  if(zeroInflation){code<-zinb_stan(priors);model<-"Zero-inflated Bayesian model"}
  if(!zeroInflation){code<-nb_stan(priors);model<-"Bayesian model without zero-inflation"}
                 
  # create data for stan use
  if (length(CF)==1) CF<-rep(CF,n)
  epg_data <-list(J=n, ystarraw = as.integer(fec), CF=CF)
  
  if (saveAll){
    savePars <- NA
  } else {
    if (zeroInflation){
      savePars <- c("mu","kappa","phi")
    } else {
      savePars <- c("mu","kappa")
    }
  }
  
  if (length(setdiff(priors,fec_setPrior()))==0){
    if(zeroInflation){stanModel<-stanmodels$zinb}
    if(!zeroInflation){stanModel<-stanmodels$nb}
  } else {
  stanModel <- stan_model(model_name=paste(model),model_code=code)}
  # whether or not to suppress progress information and errors
  if (verbose){
    samples <- sampling(stanModel, data = epg_data, pars = savePars, iter = nsamples, warmup = nburnin, chains = nchain,
                        thin = thinning, control = list(adapt_delta = adaptDelta), cores = ncore)
  } else {
  samples <- suppressMessages(suppressWarnings(sampling(stanModel, data = epg_data, pars = savePars, iter = nsamples, warmup = nburnin, chains = nchain, thin = thinning, control = list(adapt_delta = adaptDelta), cores = ncore, refresh = nsamples/4)))
  }
  
  checkDivergence(samples, adaptDelta)
  
  # calculate samples
  if(zeroInflation){
  meanEPG<-extract(samples,"mu")$mu*(1-extract(samples,"phi")$phi)
  } else {meanEPG<-extract(samples,"mu")$mu}
  
  cat("Model: ", model,"\n","Number of Samples: ",nsamples, "\n","Warm-up Samples: ",nburnin,"\n","Thinning: ",thinning,"\n","Number of Chains",nchain,"\n")
  summarys<-as.data.frame(printSummary(cbind(meanEPG)))
  
  checkConvergence(samples)
  
  return(invisible(list(stan.samples = samples,posterior.summary = summarys)))
}