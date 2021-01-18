fecr_stanExtra <- function(preFEC, postFEC, rawCounts = FALSE, preCF = 50, postCF = preCF, 
                           modelName = NULL, modelCode = NULL, modelFile = NULL, modelData = NULL, 
                           nsamples = 2000, nburnin=1000, thinning=1, nchain=2, ncore=1, adaptDelta=0.95, 
                           verbose=FALSE){
  # check if library is available
  if (!requireNamespace("eggCountsExtra")){
    stop(cat("please install the add-on package eggCountsExtra hosted on github using the command: \n devtools::install_github(\"CraigWangatCH/eggCountsExtra\")"))} else {
      requireNamespace("eggCountsExtra")
    }
  
  # check function arguments
  checkpars(nburnin, nsamples, thinning, nchain, ncore, rawCounts, adaptDelta, verbose)
  

  if (!is.null(modelName)){
    switch (modelName,
            "Po" = {paired <- TRUE; zeroInflation <- FALSE},
            "ZIPo" = {paired <- TRUE; zeroInflation <- TRUE},
            "UPo" = {paired <- FALSE; zeroInflation <- FALSE},
            "ZIUPo" = {paired <- FALSE; zeroInflation <- TRUE},
            {stop("modelName must be one of four available models, namely ('Po', 'UPo', 'ZIPo', 'ZIUPo').")}
    )
    preN <- length(preFEC)
    postN <- length(postFEC)
    checkData(preFEC, postFEC, rawCounts, preCF, postCF)
    if (length(preCF)==1) preCF<-rep(preCF,preN)
    if (length(postCF)==1) postCF<-rep(postCF,postN)
    
    preDilution <- preCF; postDilution <- postCF
    if(rawCounts){
      preDilution <- postDilution <- 1
    }
    # divide data by correction factor
    preFEC <- preFEC/preDilution
    postFEC <- postFEC/postDilution
    
    if(paired && preN != postN){
      stop("post sample size different to pre sample size\n")
    }
    
    if(paired){
      w <- eggCountsExtra::w_ratio(preFEC, postFEC)
    } else {
      w <- eggCountsExtra::w_quant(postFEC)
    }
    
    if(is.na(w$wmo)){stop("No outliers were detected, please analyze the data again with fecr_stan() function.")}
    
    cat("NOTE: it may take up to 20s for compiling the model.\n")
   
    # create data list for stan use
    epg_data <-if(paired){
      list(J=preN, ystarbraw = preFEC, ystararaw = postFEC, fpre = preCF, fpost = postCF, w = w$weight, wmo = w$wmo, postmean = w$postmean)
    } else {
      list(Jb=preN, Ja=postN, ystarbraw = preFEC, ystararaw = postFEC, fpre = preCF, fpost = postCF, w = w$weight, wmo = w$wmo, postmean = w$postmean)}
    stanModel <- stan_model(model_name = paste(modelName), model_code = eggCountsExtra::getCode(modelName))
  } else {
    if (is.null(modelData)){stop(cat("modelData argument missing, data lists required for custom models."))}
    epg_data <- modelData
    cat("NOTE: it may take up to 20s for compiling the model.\n")
    if (!is.null(modelFile)){
      stanModel <- stan_model(file = modelFile, model_name = "Custom Model")
    }
    if (!is.null(modelCode)){
      stanModel <- stan_model(model_name = "Custom Model", model_code = modelCode)
    }
  }

  if (verbose){
    samples <- sampling(stanModel, data = epg_data, iter = nsamples, warmup=nburnin, chains=nchain,
                        thin=thinning, control = list(adapt_delta = adaptDelta), cores=ncore)
  } else {
    samples <- suppressMessages(
      suppressWarnings(
        sampling(stanModel, data = epg_data, iter = nsamples,
                 warmup = nburnin, chains = nchain, thin = thinning,
                 control = list(adapt_delta = adaptDelta),cores = ncore, refresh = nsamples/4)))}

  checkDivergence(samples, adaptDelta)
  
  if (!is.null(modelName)){
  if(paired & !zeroInflation){
      meanEPG.untreated<-extract(samples,"mu")[[1]]
      meanEPG.treated<-extract(samples,"mu")[[1]]*extract(samples,"delta")$delta
      FECR<-1-extract(samples,"delta")[[1]]
      result<-cbind(FECR,meanEPG.untreated,meanEPG.treated)
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
    summarys<-as.data.frame(printSummary(result))
    }

  cat("Number of Samples: ",nsamples, "\n","Warm-up samples: ",nburnin,"\n","Thinning: ",thinning,"\n","Number of Chains",nchain,"\n")
  checkConvergence(samples)
  
  if (!is.null(modelName)){
  return(invisible(list(stan.model = stanModel, stan.samples = samples, posterior.summary = summarys)))
  } else {
    return(invisible(list(stan.model = stanModel, stan.samples = samples)))
  }
}