findglobalMLE <- function(initptinfK) {
  blob <- optimWrapper( ## purefn,
    initval=initptinfK, gr=NULL,
    chullformats=blackbox.getOption("hulls")$Kgtotal, ##
    control=list( ##  parscale is provided within optimWrapper
      fnscale=-1/blackbox.getOption("scalefactor"), trace=FALSE, maxit=10000))
  canonized <- canonizeFromKrig(blob$par)
  DemographicModel <- blackbox.getOption("DemographicModel")
  if ("IBD" %in% DemographicModel) {
    blob <- c(blob, list(latt2Ns2=canonized$latt2Ns2))
  } else if ("OnePopVarSize" %in% DemographicModel) {
    blob <- c(blob, list(Nratio=canonized$Nratio))
  } else if("OnePopFounderFlush" %in% DemographicModel) {
    blob <- c(blob, list(Nratio=canonized$Nratio), list(NactNfounderratio=canonized$NactNfounderratio), list(NfounderNancratio=canonized$NfounderNancratio))
  } else if("IM" %in% DemographicModel) {
    blob <- c(blob, list(Nratio=canonized$Nratio))
  }
  
  plotOptions <- blackbox.getOption("plotOptions")
  oneDimCIvars <- blackbox.getOption("oneDimCIvars")
  if ( length(intersect(DemographicModel, c("OnePopVarSize", "OnePopFounderFlush", "IM")))>0) {
    if ( !("IM" %in% DemographicModel) && ( ("DgmuProf" %innc% plotOptions) || ("Dgmu" %innc% oneDimCIvars) ) )  blob <- c(blob, list(Dgmu=canonized$Dgmu))
    if ( ("TgmuProf" %innc% plotOptions) || ("Tgmu" %innc% oneDimCIvars) )  blob <- c(blob, list(Tgmu=canonized$Tgmu))
  }
  if ( length(intersect(DemographicModel, c("Npop", "IM")))>0) {
    blob <- c(blob, list(NMratio=canonized$NMratio), list(mratio=canonized$mratio)) ## RL, not sur... here only to keep the goop parameter order e.g. for IM Theta, Q1, M1, M2, T,Nratio,Tmu,Mratio,miovermu,m2overmu
    if ( ("movermuProf" %innc% plotOptions) || ("m1overmu" %innc% oneDimCIvars) )  blob <- c(blob, list(m1overmu=canonized$m1overmu))
    if ( ("movermuProf" %innc% plotOptions) || ("m2overmu" %innc% oneDimCIvars) )  blob <- c(blob, list(m2overmu=canonized$m2overmu))
  }
  blob <-  c(blob, canonVP=list(canonized$canonVP))
  return(blob)
}
