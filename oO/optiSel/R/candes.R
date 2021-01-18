### Individuals not contributing to the population estimates can be selection candidates

"candes"<-function(phen, cont=NULL, N=1000, quiet=FALSE, t=NA, bc=NULL, reduce.data=TRUE, ...){

  ### Check parameter cont (contributions of age classes) #######
  if(is.null(cont)){
    cont <- data.frame(age=1, male=0.5, female=0.5, cohort=1.0)
  }else{
    cont <- checkcont(cont)
  }
  
  ### Check parameter N (population size) #######################
  N <- as.numeric(N)
  if(!is.vector(N) || length(N)>1 || any(is.na(N)) || N<=0){
    stop("N must be positive.\n")
  }
  
  ########## Define some parameters ###########################
  ra          <- cont$cohort[1]
  N0          <- ra*N
  overlap     <- ra < 0.9999
  
  
  ### Check data frame phen ###################################
  ### append columns Age, Class, c0, c1, and 
  ### an individual trait for each breed
  #############################################################
  phenAsDataTable <- "data.table" %in% class(phen)
  if(!("Born"  %in% colnames(phen)) && (!overlap)){phen$Born  <- 1}
  phen <- checkphen(phen, columns=c("Indiv", "Breed", "Born", "Sex", if(overlap){c("Sire","Dam")}), quiet=quiet, na.Sex=TRUE)
  if(is.na(t)){
    t <- max(floor(phen$Born), na.rm=TRUE)
    if(!quiet){cat(paste0("The population is evaluated at time ",t,"\n"))}
  }
  if((!overlap) && (t>=1900) && (reduce.data)){
    stop("The generations are non-overlapping, \n  so column 'Born' must contain the generation number, \n  not the year-of-birth.\n")
  }
  if(!("isCandidate" %in% colnames(phen))){
    phen$isCandidate <- TRUE
  }
  
  breed       <- characterizeBreeds(phen) 
  BreedNames  <- names(breed)
  singleBreed <- length(BreedNames)==1
  phen        <- addClasses(phen, breed, t=t)
  if(reduce.data){phen <- phen[(!is.na(phen$Age) & phen$Age>=1 & phen$Age<=nrow(cont)) | phen$isCandidate,]}
  ageClass    <- characterizeClasses(cont, breed, phen)
  phen        <- addContributions(phen, ageClass)

  Traits <- gettraits(phen, quiet=quiet)
  for(trait in Traits){
    if(all(!is.na(phen[[trait]]))&&(!singleBreed)){
      phenT <- phen[[trait]]
      phen[[trait]] <- NULL
      phen[[trait]] <- phenT
      for(b in BreedNames){
        phen[[paste(trait,b,sep=".")]] <- ifelse(phen$Breed==b, phenT, NA)
      }
    }
  }

  #### Check '...' parameters and create list with kinships ###

  kinship <- extractKinships(list(...), phen, BreedNames)

  ### Get mean kinship in each class ###
  for(i in seq_along(kinship)){
    b <- kinship[[i]]$breed
    if(b != "across breeds"){
      Classes <- ageClass$Class[ageClass$Breed==b]
      if("quadFun" %in% class(kinship[[i]])){
        kinship[[i]]$mkin  <- meankin(phen, kinship[[i]]$Q,  Classes)
      }else{
        kinship[[i]]$mkin1 <- meankin(phen, kinship[[i]]$Q1, Classes)
        kinship[[i]]$mkin2 <- meankin(phen, kinship[[i]]$Q2, Classes)
      }
    }
  }

  ### Define functions for prediciting Kin and KinatN at time t+1 ###########

  for(i in seq_along(kinship)){
    if("quadFun" %in% class(kinship[[i]])){
      kinship[[i]] <- fun4Kin(kinship[[i]], ageClass, phen, N)
    }else{
      kinship[[i]] <- fun4KinatN(kinship[[i]], ageClass, phen, N)
    }
  }
  
  ### get breed composition for multi-breed evaluations ###
  bc <- bc4candes(kinship, phen, BreedNames, bc, quiet)

  ### get current mean values of kinships and traits ###
  
  cur <- param4candes(phen, kinship, ageClass, bc, N, quiet)
  
  if(!quiet){printParam4candes(cur, BreedNames)}
  
  obj <- list(
    kinship = kinship,
    phen    = phen,
    current = cur,
    mean    = as.data.frame(as.list(setNames(cur$Val, cur$Var)),stringsAsFactors=FALSE),
    bc      = bc,
    classes = ageClass,
    breed   = breed 
  )
  
  class(obj)  <- "candes"
  
  attr(obj,"singleBreed") <- singleBreed
  
  return(obj)
}