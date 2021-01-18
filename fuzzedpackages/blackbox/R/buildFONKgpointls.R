pointsToFONK <- function(pointls,set_options=FALSE) {
  DemographicModel <- blackbox.getOption("DemographicModel")
  FONKgNames <- blackbox.getOption("FONKgNames")
  ParameterNames <- blackbox.getOption("ParameterNames")
  FONKgpointls <- pointls ## copy for modifications. Will write resulting FONKgpointls in global vars
  ## will be further modified in calcPredictorOK
  if ("Migraine" %in% blackbox.getOption("usedBy")) {
    if ("IBD" %in% blackbox.getOption("DemographicModel")) {
      latt2Ns2pt <- numeric(nrow(FONKgpointls))
      latt2Ns2pt <- apply(FONKgpointls[, 1:3], 1, tolatt2Ns2)[2, ]
      if("latt2Ns2" %in% FONKgNames) FONKgpointls[, 2] <- latt2Ns2pt ## but otherwise wemay still need it "globally"
      if (set_options) blackbox.options(latt2Ns2pt=latt2Ns2pt)
      if("condS2" %in% FONKgNames) {
        D2bool <- ("2D" %in% blackbox.getOption("DemographicModel"))
        FONKgpointls[, 3] <- apply(FONKgpointls[, 3, drop=FALSE], 1, condaxialS2fromg,D2bool=D2bool)
        ## pas de logscale ici !!
      }
    } else if ( ("Npop" %in% DemographicModel) || ("IM" %in% DemographicModel) )  {
      NMratiopt <- numeric(nrow(FONKgpointls))
      NMratiopt <- apply(FONKgpointls[, ParameterNames], 1, toNMratioFromCanonical)
      if (set_options) blackbox.options(NMratiopt=NMratiopt)
      mratiopt <- numeric(nrow(FONKgpointls))
      mratiopt <- apply(FONKgpointls[, ParameterNames], 1, tomratioFromCanonical)
      if (set_options) blackbox.options(mratiopt=mratiopt)
      m1overmupt <- numeric(nrow(FONKgpointls))
      m1overmupt <- apply(FONKgpointls[, ParameterNames], 1, tom1overmuFromCanonical)
      if (set_options) blackbox.options(m1overmupt=m1overmupt)
      m2overmupt <- numeric(nrow(FONKgpointls))
      m2overmupt <- apply(FONKgpointls[, ParameterNames], 1, tom2overmuFromCanonical)
      if (set_options) blackbox.options(m2overmupt=m2overmupt)
    } else if ( ("OnePopVarSize" %in% DemographicModel) || ("OnePopFounderFlush" %in% DemographicModel)  || ("IM" %in% DemographicModel) )  {
      Nratiopt <- numeric(nrow(FONKgpointls))
      Nratiopt <- apply(FONKgpointls[, ParameterNames], 1, toNratioFromCanonical)
      if (set_options) blackbox.options(Nratiopt=Nratiopt)
      if ( !("IM" %in% DemographicModel) && ("D" %in% ParameterNames) ) {
        Dgmupt <- numeric(nrow(FONKgpointls))
        Dgmupt <- apply(FONKgpointls[, ParameterNames], 1, toDgmuFromCanonical)
        if (set_options) blackbox.options(Dgmupt=Dgmupt)
      }
      Tgmupt <- numeric(nrow(FONKgpointls))
      Tgmupt <- apply(FONKgpointls[, ParameterNames], 1, toTgmuFromCanonical)
      if (set_options) blackbox.options(Tgmupt=Tgmupt)
      if ("OnePopFounderFlush" %in% DemographicModel) {
        NactNfounderratiopt <- numeric(nrow(FONKgpointls))
        NactNfounderratiopt <- apply(FONKgpointls[, ParameterNames], 1, toNactNfounderratioFromCanonical)
        if (set_options) blackbox.options(NactNfounderratiopt=NactNfounderratiopt)
        NfounderNancratiopt <- numeric(nrow(FONKgpointls))
        NfounderNancratiopt <- apply(FONKgpointls[, ParameterNames], 1, toNfounderNancratioFromCanonical)
        if (set_options) blackbox.options(NfounderNancratiopt=NfounderNancratiopt)
      }
    } 
  }
  colnames(FONKgpointls)[seq_len(length(FONKgNames))] <- FONKgNames
  if (set_options) {
    FONKgLow <- apply(FONKgpointls[, FONKgNames, drop=FALSE], 2, min)
    FONKgUp <- apply(FONKgpointls[, FONKgNames, drop=FALSE], 2, max)
    blackbox.options(FONKgLow = FONKgLow)
    blackbox.options(FONKgUp = FONKgUp)
    ##Note that FONKgLow/Up will be recomputed one the points have been selected for Kriging
    fittedNames <- FONKgNames[which((FONKgUp-FONKgLow)>0.00000001)]
    blackbox.options(fittedNames=fittedNames)
    blackbox.options(constantNames=FONKgNames %w/o% fittedNames)
    blackbox.options(fittedparamnbr=length(fittedNames)) ## variables retained in 'ptls <- FONKgpointls[first:last, fittedNames]' in generatePredictor()
  }
  for (st in FONKgNames) { ## FONKgScale in iterator bc FONKgScale is used by 'canonizeFromKrig'
    ## otherwise exp(.) may be returned instead of (.) in (highest lik, new points...) for constant parameters specified as logscale
    if(islogscale(st)) {FONKgpointls[, st] <- log(FONKgpointls[, st])}
  }
  infini <- apply(FONKgpointls, 1, function(v) {any(is.infinite(v))})
  if(any(infini)) {
    message.redef("infinite values in transformed points. Check input")
    message.redef("They will be deleted in further computations.")
    print(FONKgpointls[infini, ])
    FONKgpointls <- FONKgpointls[!infini, ]
  }
  if (set_options) {
    #blackbox.options(maxobsFONKy = min(FONKgpointls[, blackbox.getOption("ycolname")]))## will be used in pointfromR ... ugly coding
    blackbox.options(FONKgpointls = FONKgpointls) ##
  }
  invisible(FONKgpointls)
}

fromCanonToFONK <- function(pointsNonL) {
  DemographicModel <- blackbox.getOption("DemographicModel")
  FONKgNames <- blackbox.getOption("FONKgNames")
  ParameterNames <- blackbox.getOption("ParameterNames")
  FONKgpointls <- pointsNonL ## copy for modifications. Will write resulting FONKgpointls in global vars
  ## will be further modified in calcPredictorOK
  if ("Migraine" %in% blackbox.getOption("usedBy")) {
    if ("IBD" %in% blackbox.getOption("DemographicModel")) {
      latt2Ns2pt <- numeric(nrow(FONKgpointls))
      latt2Ns2pt <- apply(FONKgpointls[, 1:3], 1, tolatt2Ns2)[2, ]
      if("latt2Ns2" %in% FONKgNames) FONKgpointls[, 2] <- latt2Ns2pt ## but otherwise wemay still need it "globally"
      if("condS2" %in% FONKgNames) {
        D2bool <- ("2D" %in% blackbox.getOption("DemographicModel"))
        FONKgpointls[, 3] <- apply(FONKgpointls[, 3, drop=FALSE], 1, condaxialS2fromg,D2bool=D2bool)
        ## pas de logscale ici !!
      }
    } else if ( ("Npop" %in% DemographicModel) || ("IM" %in% DemographicModel) ) {
      NMratiopt <- numeric(nrow(FONKgpointls))
      NMratiopt <- apply(FONKgpointls[, ParameterNames], 1, toNMratioFromCanonical)
      mratiopt <- numeric(nrow(FONKgpointls))
      mratiopt <- apply(FONKgpointls[, ParameterNames], 1, tomratioFromCanonical)
      m1overmupt <- numeric(nrow(FONKgpointls))
      m1overmupt <- apply(FONKgpointls[, ParameterNames], 1, tom1overmuFromCanonical)
      m2overmupt <- numeric(nrow(FONKgpointls))
      m2overmupt <- apply(FONKgpointls[, ParameterNames], 1, tom2overmuFromCanonical)
    } else if ( ("OnePopVarSize" %in% DemographicModel) || ("OnePopFounderFlush" %in% DemographicModel) || ("IM" %in% DemographicModel) ) {
      Nratiopt <- numeric(nrow(FONKgpointls))
      Nratiopt <- apply(FONKgpointls[, ParameterNames], 1, toNratioFromCanonical)
      if ( !("IM" %in% DemographicModel) && ("D" %in% ParameterNames) ) {
        Dgmupt <- numeric(nrow(FONKgpointls))
        Dgmupt <- apply(FONKgpointls[, ParameterNames], 1, toDgmuFromCanonical)
      }
      Tgmupt <- numeric(nrow(FONKgpointls))
      Tgmupt <- apply(FONKgpointls[, ParameterNames], 1, toTgmuFromCanonical)
      if ("OnePopFounderFlush" %in% DemographicModel) {
        NactNfounderratiopt <- numeric(nrow(FONKgpointls))
        NactNfounderratiopt <- apply(FONKgpointls[, ParameterNames], 1, toNactNfounderratioFromCanonical)
        NfounderNancratiopt <- numeric(nrow(FONKgpointls))
        NfounderNancratiopt <- apply(FONKgpointls[, ParameterNames], 1, toNfounderNancratioFromCanonical)
      }
    }
  }
  colnames(FONKgpointls)[seq_len(length(FONKgNames))] <- FONKgNames
  for (st in FONKgNames) {if(islogscale(st)) {FONKgpointls[, st] <- log(FONKgpointls[, st])}    }
  infini <- apply(FONKgpointls, 1, function(v) {any(is.infinite(v))})
  if(any(infini)) {
    message.redef("infinite values in transformed points. Check input")
    message.redef("They will be deleted in further computations.")
    print(FONKgpointls[infini, ])
    FONKgpointls <- FONKgpointls[!infini, ]
  }
  invisible(FONKgpointls)
} ## nonL en fait

# pointls would need attributes usedBy, DemographicModel, FONKgNames, ParameterNames, ycolnames
buildFONKgpointls <- function(pointls) {pointsToFONK(pointls, set_options=TRUE)}
fromCanonToFONK <- function(pointls) {pointsToFONK(pointls, set_options=FALSE)}
