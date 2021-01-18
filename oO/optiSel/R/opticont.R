

"opticont"<-function(method, cand, con, bc=NULL,  solver="default", quiet=FALSE, make.definite=FALSE, ...){
  ################################
  ###   Check the parameters   ###
  ################################
  
  if(!("candes" %in% class(cand))){stop("Parameter 'cand' must be created with function candes.")}
  
  BreedNames   <- names(cand$breed)
  KinshipNames <- names(cand$kinship)
  TraitNames   <- gettraits(cand$phen)
  
  validMethods     <- c(outer(c("min.","max."), TraitNames, paste0), paste0("min.", KinshipNames))
  if(!(method %in% validMethods)){
    stop(paste0("The method must be one of ", paste(validMethods, collapse=", "), "\n."))
  }
  
  ####################################################
  ### Determine for which breeds the contributions ###
  ###    of males and females should be equal.     ###
  ####################################################

  considerSexes <- sapply(map(cand$breed,"Groups"), length)>1
  
  if((!all(considerSexes)) && (!quiet)){
    cat("Assuming: Sexes do not have equal contributions for \n          ")
    cat(paste0(BreedNames[!considerSexes], collapse=", "))
    cat(", because some sexes are NA.\n")
  }
  
  #####################################
  ### Get breed composition for     ###  
  ### multi-breed evaluations       ###  
  #####################################
  
  if(is.null(bc)){bc <- cand$bc}
  bc <- checkbc(BreedNames, bc)
  
  ##########################################
  ### Check constraints and              ###
  ### convert constraints to data frames ###
  ##########################################
  
  con   <- checkcon(con, cand$phen, TraitNames, KinshipNames, considerSexes, quiet=quiet)
  
  #######################################
  ### Define the optimization problem ###
  ### and solve it                    ###
  #######################################
  
  mycop     <- qcrc4cand(cand$kinship, con$df, bc)
  mycop$f   <- f4cand(cand, method, TraitNames, bc)
  mycop$max <- str_sub(method, 1, 3)=="max"
  bound     <- bounds4cand(cand$phen, con, considerSexes, bc, cand$classes, quiet=quiet)
  mycop$lb  <- lbcon(val=bound$lower, id=cand$phen$Indiv)
  mycop$ub  <- ubcon(val=bound$upper, id=cand$phen$Indiv)
  mycop$lc  <- lc4cand(cand$phen, bc, con$df, TraitNames, considerSexes)
  mycop     <- do.call(cop, mycop)
  res       <- solvecop(mycop, solver=solver, make.definite=make.definite, quiet=quiet, ...)
  sy        <- validate(mycop, res, quiet=TRUE)
  
  ######################################
  ### Transform the results back     ###
  ######################################
  
  res$info    <- sy$info
  res$summary <- summary4opticont(sy$summary, res$x, cand, con$df)
  res$solver  <- NULL
  res$status  <- NULL
  res$mean    <- setNames(res$summary$Val[-1], res$summary$Var[-1])
  res$mean    <- res$mean[names(res$mean) %in% c(KinshipNames, TraitNames)]
  res$mean    <- as.data.frame(as.list(res$mean))
  res$bc      <- as.data.frame(as.list(bc))
  res$obj.fun <- setNames(res$summary$Val[1], res$summary$Var[1])
  
  if(!quiet){
    Sy     <- res$summary
    if(identical(BreedNames, "missing")){Sy$Breed[Sy$Breed=="missing"] <- ""}
    Sy$Var    <- paste0(str_pad(Sy$Name, max(nchar(Sy$Name))+1, "right"), Sy$Breed)
    Sy        <- list(summary=Sy, info=res$info)
    class(Sy) <- "copValidation"
    print.copValidation(Sy)
  }
  
  res$parent  <- data.frame(cand$phen, lb=bound$lower, oc=res$x, ub=bound$upper)
  for(b in BreedNames){
    use <- res$parent$Breed==b
    res$parent$lb[use] <- res$parent$lb[use]/bc[b]-res$parent$c1[use]
    res$parent$oc[use] <- res$parent$oc[use]/bc[b]-res$parent$c1[use]
    res$parent$ub[use] <- res$parent$ub[use]/bc[b]-res$parent$c1[use]
    ra <- sum(res$parent$oc[use])
    if(ra>0.00001){
      res$parent$lb[use] <- res$parent$lb[use]/ra
      res$parent$oc[use] <- res$parent$oc[use]/ra
      res$parent$ub[use] <- res$parent$ub[use]/ra
    }
  }
  res$parent$ub[!is.na(res$parent$ub) & res$parent$ub>=1] <- NA
  
  res$x  <- NULL
  return(res)
  
}