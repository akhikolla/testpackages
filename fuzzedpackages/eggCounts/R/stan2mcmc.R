stan2mcmc <- function(stanFit){
  modelName <- stanFit@model_name
  switch(modelName,
         "Bayesian model without zero-inflation for paired design allowing individual efficacy"=,
         "indefficacy"={
           meanEPG.untreated <- extract(stanFit,"mu")[[1]]
           deltaMeansSample <- extract(stanFit,"delta_mu")[[1]]
           deltaShapeSample <- extract(stanFit,"delta_shape")[[1]]
           deltaSample <- qgamma(0.5, shape = deltaShapeSample, scale = deltaMeansSample/deltaShapeSample)
           FECR <- 1 - deltaSample
           meanEPG.treated <-extract(stanFit,"mu")[[1]] * deltaSample
           result <- cbind(FECR,meanEPG.untreated,meanEPG.treated)
           output <- cbind(result,as.data.frame(extract(stanFit,c("kappa","delta_mu","delta_shape"))))
         },
         "Zero-inflated Bayesian model for paired design"=,
         "zipaired"=,
         "ZIPo"={
           meanEPG.untreated<-extract(stanFit,"mu")[[1]]*(1-extract(stanFit,"phi")$phi)
           meanEPG.treated<-extract(stanFit,"mu")[[1]]*extract(stanFit,"delta")$delta*(1-extract(stanFit,"phi")$phi)
           FECR<-1-extract(stanFit,"delta")[[1]]
           result<-cbind(FECR,meanEPG.untreated,meanEPG.treated)
           output<-cbind(result,as.data.frame(extract(stanFit,c("kappa","phi","delta"))))
         },
         "Zero-inflated Bayesian model for unpaired design"=,
         "ziunpaired"=,
         "ZIUPo"={
           meanEPG.untreated<-extract(stanFit,"mu")[[1]]*(1-extract(stanFit,"phi")$phi)
           meanEPG.treated<-extract(stanFit,"mu")[[1]]*extract(stanFit,"delta")$delta*(1-extract(stanFit,"phi")$phi)
           FECR<-1-extract(stanFit,"delta")[[1]]
           result<-cbind(FECR,meanEPG.untreated,meanEPG.treated)
           output<-cbind(result,as.data.frame(extract(stanFit,c("kappa","phi","delta"))))
         },
         "Bayesian model without zero-inflation for paired design"=,
         "paired"=,
         "Po"={
           meanEPG.untreated<-extract(stanFit,"mu")[[1]]
           meanEPG.treated<-extract(stanFit,"mu")[[1]]*extract(stanFit,"delta")$delta
           FECR<-1-extract(stanFit,"delta")[[1]]
           result<-cbind(FECR,meanEPG.untreated,meanEPG.treated)
           output<-cbind(result,as.data.frame(extract(stanFit,c("kappa","delta"))))
         },
         "Simple Bayesian model without zero-inflation for paired design"=,
         "simple"={
           meanEPG.untreated<-extract(stanFit,"mu")[[1]]
           meanEPG.treated<-extract(stanFit,"mu")[[1]]*extract(stanFit,"delta")$delta
           FECR<-1-extract(stanFit,"delta")[[1]]
           result<-cbind(FECR,meanEPG.untreated,meanEPG.treated)
           output<-cbind(result,as.data.frame(extract(stanFit,c("delta"))))
         },
         "Bayesian model without zero-inflation for unpaired design"=,
         "unpaired"=,
         "UPo"={
           meanEPG.untreated<-extract(stanFit,"mu")[[1]]
           meanEPG.treated<-extract(stanFit,"mu")[[1]]*extract(stanFit,"delta")$delta
           FECR<-1-extract(stanFit,"delta")[[1]]
           result<-cbind(FECR,meanEPG.untreated,meanEPG.treated)
           output<-cbind(result,as.data.frame(extract(stanFit,c("kappa","delta"))))
         },
         "Bayesian model without zero-inflation"=,
         "nb"={
           meanEPG<-extract(stanFit,"mu")[[1]]
           kappa<-extract(stanFit,"kappa")$kappa
           output<-cbind(meanEPG=meanEPG,kappa=kappa)
         },
         "Zero-inflated Bayesian model"=,
         "zinb"={
           meanEPG<-extract(stanFit,"mu")[[1]]*(1-extract(stanFit,"phi")[[1]])
           phi<-extract(stanFit,"phi")$phi
           kappa<-extract(stanFit,"kappa")$kappa
           output<-cbind(meanEPG=meanEPG,kappa=kappa,phi=phi)
         }
        )
  return(invisible(mcmc(output,thin=stanFit@sim$thin)))
}