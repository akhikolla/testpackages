tofullKrigingspace <- function(fittedlist, fixedlist=NULL) { ## output in kriging space (reorder columns in standard order if required)
  fittedparamnbr <- blackbox.getOption("fittedparamnbr")
  fittedNames <- blackbox.getOption("fittedNames")
  ParameterNames <- blackbox.getOption("ParameterNames")
  getvalue <- function(st) {
    if (st %in% names(fittedlist)) {
      return(as.numeric(fittedlist[st]))
    } else if (st %in% names(fixedlist)) {
      return(as.numeric(fixedlist[st]))
    } else if (st %in% blackbox.getOption("constantNames")) { ## we may need values of fixed parameters
      return(as.numeric(blackbox.getOption("FONKgLow")[st]))
    } else return(NA)
  }
  ## fittedlist and fixedlist must be in FONKgScale or else in extraScale
  ## (01/2010: 'full' Krig sp of *fittedparamnbr* (the argument of purefn)
  ## Generates a 'full' vector from fittedlist and fixedlist. fixedlist can be NULL although few subcases are implemented
  ## check fittedlist
  if(is.null(names(fittedlist))) {
    stop.redef("(!) names(fittedlist) is NULL in tofullKrigingspace. ")
  }
  ## FR 18/11/10 replaced fittedNames by ParameterNames in next line, for case where latt2Ns2 in fittedNames and twoNm in fixedlist (for LRT on twoNm...)
  ## Not sure whether this will be OK with fittedparamnbr< param nbr
  checkfitted <- names(fittedlist) %w/o% c(ParameterNames, "latt2Ns2", "NMratio", "mratio", "m1overmu", "m2overmu", "Nratio", "NactNfounderratio", "NfounderNancratio", "condS2", "Dgmu", "Tgmu")
  if (length(checkfitted)>0) {
    message.redef("(!)From tofullKrigingspace(): names(fittedlist) contains unhandled variable")
    message.redef("   or combination of variables") ##latt2Ns2 plus another variable
    message.redef(checkfitted) ## hind: dont forget to use <name>=<vector>[[<name>]], ie [[ ]] not [ ]
    stop.redef()
  }
  ## check fixedlist
  if(!is.null(fixedlist)) {
    checkfixed <- names(fixedlist) %w/o% c(ParameterNames, "latt2Ns2", "NMratio", "mratio", "m1overmu", "m2overmu", "Nratio", "NactNfounderratio", "NfounderNancratio", "condS2", "Dgmu", "Tgmu")
    if (length(checkfixed)>0) {
      message.redef("(!)From tofullKrigingspace(): names(fixedlist) contains unhandled variable")
      message.redef("   or combination of variables") ##latt2Ns2 plus another variable
      message.redef(paste(names(fixedlist)))
      stop.redef()
    }
  }
  KrigVec <- rep(NA, fittedparamnbr)
  names(KrigVec) <- fittedNames
  ## if everything is in canonical space then the next two lines fill the KrigVec vector
  finf <- intersect(names(fixedlist), fittedNames)
  KrigVec[finf] <- as.numeric(fixedlist[finf])
  finf <- intersect(names(fittedlist), fittedNames)
  KrigVec[finf] <- as.numeric(fittedlist[finf])
  if (!any(is.na(KrigVec))) { ## then we are done
    return(KrigVec)
  } ##ELSE
  ## we have to handle composite variables. We first unlog everything that is logscale
  ## we must be careful that even if say extrascale=Nb=logscale, a non-log latt2Ns2 can be given in the arguments to tofullKrigingspace...
  for(st in names(fittedlist)) if (islogscale(st)) {fittedlist[[st]] <- exp(fittedlist[[st]])}
  for(st in names(fixedlist)) if (islogscale(st)) {fixedlist[[st]] <- exp(fixedlist[[st]])}
  for(st in names(KrigVec)) if (islogscale(st)) {KrigVec[[st]] <- exp(KrigVec[[st]])} ## exp(NA)->NA, not a problem
  D2bool <- ("2D" %in% blackbox.getOption("DemographicModel"))
  ## then we operate in canonical scale
  if("twoNm" %in% fittedNames && is.na(KrigVec["twoNm"])) { ##
    latt2Ns2 <- getvalue("latt2Ns2")
    if(is.na(latt2Ns2)) stop.redef("(!) From tofullKrigingspace(): neither twoNm nor latt2Ns2 given")
    g <- getvalue("g")
    if (is.na(g)) {
      S2 <- getvalue("condS2")
      if (is.na(S2)) stop.redef("(!) From tofullKrigingspace(): neither twoNm nor g nor condS2 given")
    } else S2 <- condaxialS2fromg(g, D2bool=D2bool)
    KrigVec["twoNm"] <- latt2Ns2/S2
  }
  if("latt2Ns2" %in% fittedNames && is.na(KrigVec["latt2Ns2"])) { ## then we must have twoNm somewhere
    twoNm <- getvalue("twoNm")
    if(is.na(twoNm)) stop.redef("(!) From tofullKrigingspace(): neither twoNm nor latt2Ns2 given")
    S2 <- getvalue("condS2")
    if (is.na(S2)) {
      S2 <- condaxialS2fromg(KrigVec["g"], D2bool=D2bool) ##cond axial S2 from g
    }
    KrigVec["latt2Ns2"] <- twoNm*S2 ## 2Nm as latt2Ns2/condS2
  }
  if("condS2" %in% fittedNames && is.na(KrigVec["condS2"])) {
    g <- getvalue("g")
    if(is.na(g)) { ## need to reconstruct it from other information
      twoNm <- getvalue("twoNm")
      if(is.na(twoNm)) stop.redef("(!) From tofullKrigingspace(): neither twoNm nor condS2 given")
      latt2Ns2 <- getvalue("latt2Ns2")
      if(is.na(latt2Ns2))  stop.redef("(!) From tofullKrigingspace(): neither twoNm nor latt2Ns2 given")
      KrigVec["condS2"] <- latt2Ns2/twoNm
    } else KrigVec["condS2"] <- condaxialS2fromg(g, D2bool=D2bool) ##cond axial S2 from g
  }
  if("g" %in% fittedNames && is.na(KrigVec["g"])) { ## then we must have twoNm and latt2Ns2 somewhere
    twoNm <- getvalue("twoNm")
    S2 <- getvalue("condS2")
    if (is.na(S2)) {
      latt2Ns2 <- getvalue("latt2Ns2")
      if (is.na(twoNm) || is.na(latt2Ns2) ) stop.redef("(!) From tofullKrigingspace(): g and [either twoNm or latt2Ns2] are not given")
      S2 <- latt2Ns2/twoNm ##cond axial S2 from g
    }
    KrigVec["g"] <- groot(S2, D2bool=D2bool )
  }
  if("D" %in% fittedNames && is.na(KrigVec["D"])) { ## then we must have Dgmu somewhere
    Dgmu <- getvalue("Dgmu")
    if(is.na(Dgmu)) {
      stop.redef("(!) From tofullKrigingspace(): neither D nor Dgmu given")
    } else {
      if(is.na(KrigVec["twoNmu"])) { ## RL 022017 means that twoNmu was not fitted, yet impossible...
        KrigVec["D"] <- Dgmu/blackbox.getOption("FONKgLow")["twoNmu"]
      } else {KrigVec["D"] <- Dgmu/KrigVec["twoNmu"]} ## twoNmu is already unlog'ed above and D will be relog'ed below
    }
  }
  if("T" %in% fittedNames && is.na(KrigVec["T"])) { ## then we must have Tgmu somewhere
    Tgmu <- getvalue("Tgmu")
    if(is.na(Tgmu)) {
      stop.redef("(!) From tofullKrigingspace(): neither T nor Tgmu given")
    } else {
      if(is.na(KrigVec["twoNmu"])) { ## RL 022017 means that twoNmu was not fitted, yet impossible...
        KrigVec["T"] <- Tgmu/blackbox.getOption("FONKgLow")["twoNmu"]
      } else {KrigVec["T"] <- Tgmu/KrigVec["twoNmu"]} ## twoNmu is already unlog'ed above and T will be relog'ed below
    }
  }
  if("twoNmu" %in% fittedNames && is.na(KrigVec["twoNmu"])) { ## then we must have Nratio somewhere
    Nratio <- getvalue("Nratio")
    if(is.na(Nratio)) {
      NactNfounderratio <- getvalue("NactNfounderratio")
      if(is.na(NactNfounderratio)) stop.redef("(!) From tofullKrigingspace(): neither twoNmu nor Nratio nor NactNfounderratio given")
      if(is.na(KrigVec["twoNfoundermu"])) { ## RL 052013 means that twoNfoundermu was not fitted
        KrigVec["twoNmu"] <- NactNfounderratio*blackbox.getOption("FONKgLow")["twoNfoundermu"]
      } else {KrigVec["twoNmu"] <- NactNfounderratio*KrigVec["twoNfoundermu"]} ## twoNfoundermu is already unlog'ed above and twoNmu will be relog'ed below
    } else {
      if(is.na(KrigVec["twoNancmu"])) { ## RL 052013 means that twoNancmu was not fitted
        KrigVec["twoNmu"] <- Nratio*blackbox.getOption("FONKgLow")["twoNancmu"]
      } else {KrigVec["twoNmu"] <- Nratio*KrigVec["twoNancmu"]} ## twoNancmu is already unlog'ed above and twoNmu will be relog'ed below
    }
  }
  if("twoNancmu" %in% fittedNames && is.na(KrigVec["twoNancmu"])) { ## then we must have Nratio somewhere
    Nratio <- getvalue("Nratio")
    if(is.na(Nratio)) {
      NfounderNancratio <- getvalue("NfounderNancratio")
      if(is.na(NfounderNancratio)) stop.redef("(!) From tofullKrigingspace(): neither twoNancmu nor Nratio nor NfounderNancratio given")
      if(is.na(KrigVec["twoNfoundermu"])) { ## RL 022016 means that twoNfoundermu was not fitted
        KrigVec["twoNancmu"] <- blackbox.getOption("FONKgLow")["twoNfoundermu"]/NfounderNancratio
      } else {KrigVec["twoNancmu"] <- KrigVec["twoNfoundermu"]/NfounderNancratio} ## twoNfoundermu is already unlog'ed above and twoNmu will be relog'ed below
    } else {
      if(is.na(KrigVec["twoNmu"])) { ## RL 022016 means that twoNmu was not fitted
        KrigVec["twoNancmu"] <- blackbox.getOption("FONKgLow")["twoNmu"]/Nratio
      } else {KrigVec["twoNancmu"] <- KrigVec["twoNmu"]/Nratio} ## twoNmu is already unlog'ed above and twoNancmu will be relog'ed below
    }
  }
  if("M1" %in% fittedNames && is.na(KrigVec["M1"])) { ## then we must have m1overmu somewhere
    m1overmu <- getvalue("m1overmu")
    if(is.na(m1overmu)) {
      NMratio <- getvalue("NMratio")
      if(is.na(NMratio)) {
        mratio <- getvalue("mratio")
        if(is.na(mratio)) {
          stop.redef("(!) From tofullKrigingspace(): neither M1 nor m1overmu nor NMratio not mratio given")
          } else {
            if(is.na(KrigVec["M2"]) && is.na(KrigVec["Q1"])) { ## RL 022017 means that neither twoNmu nor Q1 was not fitted, yet impossible...
              KrigVec["M1"] <- mratio*blackbox.getOption("FONKgLow")["M2"]*blackbox.getOption("FONKgLow")["Q1"]/(1.0-blackbox.getOption("FONKgLow")["Q1"])
            } else if(is.na(KrigVec["M2"])) { ## RL 022017 means that twoNmu was not fitted, yet impossible...
              KrigVec["M1"] <- mratio*blackbox.getOption("FONKgLow")["M2"]*KrigVec["Q1"]/(1.0-KrigVec["Q1"])
            } else if(is.na(KrigVec["Q1"])) { ## RL 022017 means that twoNmu was not fitted, yet impossible...
              KrigVec["M1"] <- mratio*KrigVec["M2"]*blackbox.getOption("FONKgLow")["Q1"]/(1.0-blackbox.getOption("FONKgLow")["Q1"])
            } else {KrigVec["M1"] <- mratio*KrigVec["M2"]*KrigVec["Q1"]/(1.0-KrigVec["Q1"])} ## M2 and Q1 are already unlog'ed above and M1 will be relog'ed below
          }
      } else {
        if(is.na(KrigVec["M2"])) { ## RL 022017 means that twoNmu was not fitted, yet impossible...
          KrigVec["M1"] <- NMratio*blackbox.getOption("FONKgLow")["M2"]
        } else {KrigVec["M1"] <- NMratio*KrigVec["M2"]} ## M2 is already unlog'ed above and M1 will be relog'ed below
      }
    } else {
      if(is.na(KrigVec["twoNmu"]) && is.na(KrigVec["Q1"])) { ## RL 022017 means that twoNmu nor Q1 was not fitted, yet impossible...
        KrigVec["M1"] <- m1overmu*blackbox.getOption("FONKgLow")["twoNmu"]*blackbox.getOption("FONKgLow")["Q1"]
      } else if(is.na(KrigVec["twoNmu"])) { ## RL 022017 means that twoNmu was not fitted, yet impossible...
        KrigVec["M1"] <- m1overmu*blackbox.getOption("FONKgLow")["twoNmu"]*KrigVec["Q1"]
      } else if(is.na(KrigVec["Q1"])) { ## RL 022017 means that twoNmu was not fitted, yet impossible...
        KrigVec["M1"] <- m1overmu*KrigVec["twoNmu"]*blackbox.getOption("FONKgLow")["Q1"]
      } else {KrigVec["M1"] <- m1overmu*KrigVec["twoNmu"]*KrigVec["Q1"]} ## twoNmu is already unlog'ed above and M1 will be relog'ed below
    }
  }
  if("M2" %in% fittedNames && is.na(KrigVec["M2"])) { ## then we must have m2overmu somewhere
    m2overmu <- getvalue("m2overmu")
    if(is.na(m2overmu)) {
      stop.redef("(!) From tofullKrigingspace(): neither M2 nor m2overmu given")
    } else {
      if(is.na(KrigVec["twoNmu"]) && is.na(KrigVec["Q1"])) { ## RL 022017 means that twoNmu nor Q1 was not fitted, yet impossible...
        KrigVec["M2"] <- m2overmu*blackbox.getOption("FONKgLow")["twoNmu"]*(1.0-blackbox.getOption("FONKgLow")["Q1"])
      } else if(is.na(KrigVec["twoNmu"])) { ## RL 022017 means that twoNmu was not fitted, yet impossible...
        KrigVec["M2"] <- m2overmu*blackbox.getOption("FONKgLow")["twoNmu"]*(1.0-KrigVec["Q1"])
      } else if(is.na(KrigVec["Q1"])) { ## RL 022017 means that twoNmu was not fitted, yet impossible...
        KrigVec["M2"] <- m2overmu*KrigVec["twoNmu"]*(1.0-blackbox.getOption("FONKgLow")["Q1"])
      } else {KrigVec["M2"] <- m2overmu*KrigVec["twoNmu"]*(1.0-KrigVec["Q1"])} ## twoNmu is already unlog'ed above and M1 will be relog'ed below
    }
  }
  if(!all(is.numeric(KrigVec))) {
    message.redef("(!) From tofullKrigingspacefn(): !all(is.numeric(KrigVec))")
    llocalst <- paste("fittedlist was ", fittedlist)
    message.redef(llocalst)
    llocalst <- paste("fixedlist was ", fixedlist)
    message.redef(llocalst)
    llocalst <- paste("KrigVec was ", KrigVec)
    message.redef(llocalst)
  }
  ## Finally we relog everything that is logscale
  for(st in names(KrigVec)) if ((st %in% fittedNames) && islogscale(st)) {KrigVec[[st]] <- log(KrigVec[[st]])}
  return(KrigVec)
} ## end tofullKrigingspace()
