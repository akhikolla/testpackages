
"bounds4cand" <- function(phen, con, considerSexes, bc, ageClass, quiet=FALSE){
  lb   <- con[["lb"]]
  ub   <- con[["ub"]]
  unif <- con[["uniform"]]
  BreedNames <- names(considerSexes)
  pOff <- ageClass[phen$Class,"propOff"]/ageClass[phen$Class,"n"]
  pOff[is.na(pOff)] <- 0
  
  ### If males, females or individuals from one breed should have unif (equal)
  ### numbers of offspring, then the contributions they have as parents
  ### are added to the lower bounds. The upper bounds are set equal to the lower bounds.
  lbval <- setNames(phen$c1,            phen$Indiv)
  ubval <- setNames(rep(1, nrow(phen)), phen$Indiv)
  
  for(b in BreedNames){
    isBreedb <- phen$Breed==b
    ra       <- 1 - sum(phen$c1[isBreedb])
    if(!is.null(unif)){
      if(considerSexes[b]  && any(unif$Breed==b & unif$Sex %in% "female")){
        use        <- isBreedb & phen$Sex=="female" 
        lbval[use] <- lbval[use] + ra*0.5*pOff[use]/sum(pOff[use])
        ubval[use] <- lbval[use]
      }
      if(considerSexes[b]  && any(unif$Breed==b & unif$Sex %in% "male")){
        use        <- isBreedb & phen$Sex=="male"
        lbval[use] <- lbval[use] + ra*0.5*pOff[use]/sum(pOff[use])
        ubval[use] <- lbval[use]
      }
      if((!considerSexes[b])  && (b %in% unif$Breed)){
        use        <- isBreedb
        lbval[use] <- lbval[use] + ra*pOff[use]/sum(pOff[use])
        ubval[use] <- lbval[use]
      }
    }
    if(!is.null(ub)){
      ubb <- ub[names(ub) %in% phen$Indiv[isBreedb]]
      ubval[names(ubb)] <- phen[names(ubb),"c1"] + ra*ubb
    }
    if(!is.null(lb)){
      lbb <- lb[names(lb) %in% phen$Indiv[isBreedb]]
      lbval[names(lbb)] <- phen[names(lbb),"c1"] + ra*lbb
    }
  }

  ubval[!phen$isCandidate] <- lbval[!phen$isCandidate]
  
  return(list(upper=ubval*bc[phen$Breed], lower=lbval*bc[phen$Breed]))
}