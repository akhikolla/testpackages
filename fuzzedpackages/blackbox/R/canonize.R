## wrapper around canonizeFromKrig, accepting matrix input and extracting canonVP 
## input has dim. Since 2016/01/05, output always has dim
toCanonical <- function(candidates, ## in (possibly incomplete) Kriging space 
                        FONKgLow,
                        otherlist=NULL ## for completion from CI info 
                        ) {
  INFO <- blackbox.options()[c("ParameterNames","fittedNames")]
  if  (is.matrix(candidates)) {
    candidates <- apply(candidates, 1, tofullKrigingspace,fixedlist=otherlist)
    if (length(INFO$fittedNames)>1L) {
      candidates <- t(candidates)
    } else candidates <- matrix(candidates, ncol=1)
    ## apply loses $canonVP names and instead copies the candidates' names...
    colnames(candidates) <- INFO$fittedNames
    candidates <- apply(candidates, 1, function(v) {canonizeFromKrig(v)$canonVP}) ## transposed // expected; except if fittedparamnbr==1...
    if (length(INFO$ParameterNames)>1L) {
      candidates <- t(candidates)
    } else candidates <- matrix(candidates, ncol=1)
    ## apply loses $canonVP names and instead copies the candidates' names...
    colnames(candidates) <- INFO$ParameterNames
  } else {
    candidates <- tofullKrigingspace(candidates,fixedlist=otherlist)
    candidates <- t(canonizeFromKrig(candidates)$canonVP) ## t() converts to matrix with the column names
  }
  return(candidates)
}



canonizeFromKrig <- function(input) { ## from vector in complete Kriging space AND in Kriging scale, also returns composite var
  INFO <- blackbox.options()[c("FONKgLow","ParameterNames","DemographicModel","plotOptions", "oneDimCIvars", "FONKgNames","FONKgScale")]
  FONKinput <- INFO$FONKgLow ## initial value
  if (length(setdiff(names(input),names(FONKinput)))>0L) {
    stop("input argument of canonizeFromKrig() is invalid (should be within kriging space).")
  }
  FONKinput[names(input)] <- input
  FONKinput <- unlist(FONKinput) ## previous line sometimes creates a list although no argument is ??
  ## unlogs
  logs <- tolower(INFO$FONKgScale)=="logscale" ## vector of T/F
  FONKinput[logs] <- exp(FONKinput[logs]) ## argument of from2Ns2Tocanon should not be logscale....
  canon <- FONKinput
  names(canon) <-INFO$ParameterNames ## anticipate future names
  DemographicModel <- INFO$DemographicModel
  FONKgNames <- INFO$FONKgNames
  if ("IBD" %in% DemographicModel) {
    if("condS2" %in% FONKgNames) {
      canon["g"] <- groot(FONKinput["condS2"], D2bool= ("2D" %in% DemographicModel) )
    }
    if("latt2Ns2" %in% FONKgNames) {
      latt2Ns2 <- FONKinput["latt2Ns2"] ## saved in return value of canonizeFromKrig
      canon["twoNm"] <- from2Ns2Tocanon(FONKinput)["twoNm"]
    } else {## constructs 2Ds2 [lattice units]
      latt2Ns2 <- (tolatt2Ns2(canon))["latt2Ns2"] ## requires that canon is indeed already canonical
    }
    return(list(canonVP=canon, latt2Ns2=latt2Ns2))
    } else if ( length(intersect(DemographicModel, c("OnePopVarSize", "OnePopFounderFlush", "IM")))>0) {
    if("Nratio" %in% FONKgNames) {
      Nratio <- FONKinput["Nratio"] ## saved in return value of canonizeFromKrig
      canon["twoNmu"] <- FONKinput[["Nratio"]]*FONKinput[["twoNancmu"]]
    } else {## constructs Nratio
      Nratio <- toNratioFromCanonical(canon) ## requires that canon is indeed already canonical
    }
    paramValList=list(canonVP=canon, Nratio=Nratio)
    if ("OnePopFounderFlush" %in% DemographicModel) {
      if("NactNfounderratio" %in% FONKgNames) {
        NactNfounderratio <- FONKinput["NactNfounderratio"] ## saved in return value of canonizeFromKrig
        canon["twoNmu"] <- FONKinput[["NactNfounderratio"]]*FONKinput[["twoNfoundermu"]]
      } else {## constructs NactNfounderratio
        NactNfounderratio <- toNactNfounderratioFromCanonical(canon) ## requires that canon is indeed already canonical
      }
      if("NfounderNancratio" %in% FONKgNames) {
        NfounderNancratio <- FONKinput["NfounderNactratio"] ## saved in return value of canonizeFromKrig
        canon["twoNancmu"] <- FONKinput[["twoNfoundermu"]]/FONKinput[["NfounderNancratio"]]
      } else {## constructs NfounderNancratio
        NfounderNancratio <- toNfounderNancratioFromCanonical(canon) ## requires that canon is indeed already canonical
      }
      paramValList=c(paramValList,list(NactNfounderratio=NactNfounderratio, NfounderNancratio=NfounderNancratio))
    }
    if ( ( !("IM" %in% DemographicModel) && ("D" %in% INFO$ParameterNames) ) && ( ("DgmuProf" %innc% INFO$plotOptions) || ("Dgmu" %innc% INFO$oneDimCIvars) ) )  {
      if("Dgmu" %in% FONKgNames) {
        Dgmu <- FONKinput["Dgmu"] ## saved in return value of canonizeFromKrig
        canon["D"] <- FONKinput[["Dgmu"]]/FONKinput[["twoNmu"]]
      } else {## constructs Dgmu
        Dgmu <- toDgmuFromCanonical(canon) ## requires that canon is indeed already canonical
      }
      paramValList=c(paramValList,list(Dgmu=Dgmu))
    }
    if ( ("TgmuProf" %innc% INFO$plotOptions) || ("Tgmu" %innc% INFO$oneDimCIvars) )  {
      if("Tgmu" %in% FONKgNames) {
        Tgmu <- FONKinput["Tgmu"] ## saved in return value of canonizeFromKrig
        canon["T"] <- FONKinput[["Tgmu"]]/FONKinput[["twoNmu"]]
      } else {## constructs Tgmu
        Tgmu <- toTgmuFromCanonical(canon) ## requires that canon is indeed already canonical
      }
      paramValList=c(paramValList,list(Tgmu=Tgmu))
    }
    return(paramValList)
    } else if ( length(intersect(DemographicModel, c("Npop", "IM")))>0) {
      if("NMratio" %in% FONKgNames) {
        NMratio <- FONKinput["NMratio"] ## saved in return value of canonizeFromKrig
        canon["M1"] <- FONKinput[["NMratio"]]*FONKinput[["M2"]]
      } else {## constructs NMratio
        NMratio <- toNMratioFromCanonical(canon) ## requires that canon is indeed already canonical
      }
      if("mratio" %in% FONKgNames) {
        mratio <- FONKinput["mratio"] ## saved in return value of canonizeFromKrig
        canon["M1"] <- FONKinput[["mratio"]]*FONKinput[["M2"]]*FONKinput[["Q1"]]/(1.0-FONKinput[["Q1"]])
      } else {## constructs mratio
        mratio <- tomratioFromCanonical(canon) ## requires that canon is indeed already canonical
      }
      paramValList=list(canonVP=canon, NMratio=NMratio, mratio=mratio)
      if ( ("movermuProf" %innc% INFO$plotOptions) || ("m1overmu" %innc% INFO$oneDimCIvars) )  {
        if("m1overmu" %in% FONKgNames) {
          m1overmu <- FONKinput["m1overmu"] ## saved in return value of canonizeFromKrig
          canon["M1"] <- FONKinput[["m1overmu"]]*FONKinput[["twoNmu"]]*FONKinput[["Q1"]]
        } else {## constructs m1overmu
          m1overmu <- tom1overmuFromCanonical(canon) ## requires that canon is indeed already canonical
        }
        paramValList=c(paramValList,list(m1overmu=m1overmu))
      }
      if ( ("movermuProf" %innc% INFO$plotOptions) || ("m2overmu" %innc% INFO$oneDimCIvars) )  {
        if("m2overmu" %in% FONKgNames) {
          m2overmu <- FONKinput["m2overmu"] ## saved in return value of canonizeFromKrig
          canon["M2"] <- FONKinput[["m2overmu"]]*FONKinput[["twoNmu"]]*(1.0-FONKinput[["Q1"]])
        } else {## constructs m2overmu
          m2overmu <- tom2overmuFromCanonical(canon) ## requires that canon is indeed already canonical
        }
        paramValList=c(paramValList,list(m2overmu=m2overmu))
      }
      return(paramValList)
    } else return(list(canonVP=canon))
} ## end def canonizeFromKrig
