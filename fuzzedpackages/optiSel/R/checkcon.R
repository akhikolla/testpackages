
"checkcon" <- function(con, phen, TraitNames, KinshipNames, considerSexes, quiet){
  BreedNames <- names(considerSexes)
  
  validConstraints <- c(outer(c("lb.","ub."),   TraitNames, paste0), paste0("ub.",  KinshipNames), "ub", "lb", "uniform")
  if(!all(names(con) %in% validConstraints)){
    report <- paste(setdiff(names(con), validConstraints), collapse=", ")
    stop(paste0("The following constraints do not have valid names: ", report, "\n."))
  }
  
  available4OCS <- phen$isCandidate
  
  ### Check constraint 'uniform' ###
  unif <- NULL
  if("uniform" %in% names(con)){
    unif <- con[["uniform"]]
    if(!is.vector(unif)){stop("Constraint 'uniform' must be a character vector.\n")}
    unif <- as.character(unif)
    if(("male"%in% unif) && any(!considerSexes)){
      stop("String 'male' in constraint 'uniform' is only permitted if there are no NA-values in column Sex of cand$phen.\n")
    }
    if(("female"%in% unif) && any(!considerSexes)){
      stop("String 'female' in constraint 'uniform' is only permitted if there are no NA-values in column Sex of cand$phen.\n")
    }  
    if(any( ((paste0(BreedNames, ".female") %in% unif) & !considerSexes))){
      stop("String 'BREEDNAME.female' in constraint 'uniform' is only permitted for breeds that have no NA-values in column Sex of cand$phen.\n")
    }
    if(any( ((paste0(BreedNames, ".male") %in% unif) & !considerSexes))){
      stop("String 'BREEDNAME.male' in constraint 'uniform' is only permitted for breeds that have no NA-values in column Sex of cand$phen.\n")
    }
    validUnif <- c(BreedNames, "male", "female")
    if(any(considerSexes)){
      validUnif <- c(validUnif, outer(BreedNames[considerSexes], c(".female",".male"), paste0))
    }
    if(!all(unif %in% validUnif)){
      report <- paste(setdiff(unif, validUnif), collapse="', '")
      stop(paste0("The following strings in constraint 'uniform' are not meaningful: '", report, "'\n"))
    }
    if("female" %in% unif){
      unif <- c(setdiff(unif, "female"), paste0(BreedNames,".female"))
    }
    if("male" %in% unif){
      unif <- c(setdiff(unif, "male"), paste0(BreedNames,".male"))
    }
    for(b in BreedNames){
       if((b %in% unif) && (considerSexes[b])){
        unif <- c(paste0(b,".male"), paste0(b,".female"), setdiff(unif, b))
       }
      if((b %in% unif) && (!considerSexes[b])){
        unif <- c(paste0(b,".NA"), setdiff(unif, b))
      }
    }
    unif <- as.data.frame(tstrsplit(unique(unif), "\\."), 
                          col.names        = c("Breed","Sex"), 
                          stringsAsFactors = FALSE)
    unif$Sex[unif$Sex=="NA"] <- NA
    
    for(i in seq_along(unif$Breed)){
      ### Check if all individuals with uniform    ###
      ### contributions are selection candidates   ###
      use <- phen$Breed==unif$Breed[i]
      if(!is.na(unif$Sex[i])){use <- use & phen$Sex == unif$Sex[i]}
      available4OCS <- available4OCS & !use
      #if(any(!phen$isCandidate[use])){
      #  if(!quiet){cat(paste0("Column 'isCandidate' of phen is ignored for ", unif$Sex[i], " individuals from ", unif$Breed[i], ".\n"))}
      #}
    }
    con[["uniform"]] <- unif
  }
  
  
  Indiv4OCS <- phen$Indiv[available4OCS]
  
  ### Check constraint 'ub' ###
  ub <- NULL
  if("ub" %in% names(con)){
    ub <- con[["ub"]]
    if(!is.vector(ub) || !is.numeric(ub)){stop("Constraint ub must be a named numeric vector.\n")}
    if(is.null(names(ub))){stop("Constraint ub must be a named numeric vector.\n")}
    ub <- ub[!is.na(ub)]
    if(!all(names(ub) %in% phen$Indiv[phen$isCandidate])){
      reportIndiv <- paste(setdiff(names(ub), phen$Indiv[phen$isCandidate]), collapse=", ")
      stop(paste0("The following individuals from constraint 'ub' do not appear as candidates in data frame cand$phen:", reportIndiv,".\n"))
    }
    if(!all(names(ub) %in% Indiv4OCS)){
      reportIndiv <- paste(setdiff(names(ub), Indiv4OCS), collapse=", ")
      stop(paste0("For the following individuals from constraint 'ub' an upper bound cannot be defined due to constraint 'uniform':", reportIndiv,".\n"))
    }
    if(any(ub<0)){stop("Some upper bounds are negative.\n")}
    if(length(ub)==0){ub <- NULL}
    con[["ub"]] <- ub
  }

  ### Check constraint 'lb' ###
  lb <- NULL
  if("lb" %in% names(con)){
    lb <- con[["lb"]]
    if(!is.vector(lb) || !is.numeric(lb)){stop("Constraint lb must be a named numeric vector.\n")}
    if(is.null(names(lb))){stop("Constraint lb must be a named numeric vector.\n")}
    lb <- lb[!is.na(lb)]
    if(!all(names(lb) %in% phen$Indiv[phen$isCandidate])){
      reportIndiv <- paste(setdiff(names(lb), phen$Indiv[phen$isCandidate]), collapse=", ")
      stop(paste0("The following individuals from constraint 'lb' do not appear as candidates in data frame cand$phen:", reportIndiv,".\n"))
    }
    if(!all(names(lb) %in% Indiv4OCS)){
      reportIndiv <- paste(setdiff(names(lb), Indiv4OCS), collapse=", ")
      stop(paste0("For the following individuals from constraint 'lb' an lower bound cannot be defined due to constraint 'uniform':", reportIndiv,".\n"))
    }
    if(any(lb<0)){stop("Some lower bounds are negative.\n")}
    if(length(lb)==0){lb <- NULL}
    con[["lb"]] <- lb
  }
  
  if(("lb" %in% names(con)) && ("ub" %in% names(con))){
    both <- intersect(names(lb), names(ub))
    if(length(both)>0 && any(lb[both]>ub[both])){
      reportIndiv <- paste(both[lb[both]>ub[both]], collapse=", ")
      stop(paste0("For the following individuals the lower bound is larger than the upper bound:", reportIndiv,".\n"))
    }
  }
    
  con <- list(
    uniform = unif,
    lb      = lb,
    ub      = ub,
    df      = con2df(con, TraitNames)
  )
  
  return(con)
}