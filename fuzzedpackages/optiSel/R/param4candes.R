
"param4candes" <- function(phen, kinship, ageClass, bc, N, quiet){
  BreedNames <- names(bc)
  Trait <- gettraits(phen, quiet=TRUE)

  ### Data frame with current trait values ###
  curT <- NULL
  if(length(Trait)>0){
    if(length(BreedNames)==1){
      BreedforTrait <- rep(BreedNames, length(Trait))
    }else{
      BreedforTrait <- rep("across breeds", length(Trait))
      for(t in seq_along(Trait)){
        if(any(is.na(phen[[Trait[t]]]))){
          BreedforTrait[t] <- unique(phen$Breed[!is.na(phen[[Trait[t]]])])
        }
      }
    }   
    
    curT <- data.frame(
      Type     = "trait",
      Var       = Trait, 
      Breed     = BreedforTrait, 
      Val       = NA,
      obj.fun.and.constraints = '',
      stringsAsFactors=FALSE)

    for(i in seq_along(curT$Var)){
      v <- curT$Var[i]
      b <- curT$Breed[i]
      if(b == "across breeds"){
        curT[i,"Val"] <- sum(bc[phen$Breed]*phen$c0*phen[[v]])
      }else{
        curT[i,"Val"] <- sum((phen$c0*phen[[v]])[phen$Breed==b])
      }
      curT[i,"obj.fun.and.constraints"] <- paste0(c("min.", "max.", "lb.","eq.", "ub."), curT$Var[i], collapse=", ")
    }
  }
  
  
  ### Data frame with current kinships and native kinships ###
  
  curK <- data.frame(
    Type  = unlist(lapply(kinship,function(x){class(x)})),
    Var   = names(kinship),
    Breed = unlist(map(kinship,"breed")),
    Val   = NA,
    obj.fun.and.constraints = '',
    stringsAsFactors=FALSE)
  
  curK$Type <- mapvalues(curK$Type, c("quadFun", "ratioFun"), c("kinship", "nat. kin."), warn_missing = FALSE)

  for(i in seq_along(curK$Var)){
    v <- curK$Var[i]
    b <- curK$Breed[i]
    
    if(b == "across breeds"){
      u <- bc[phen$Breed]*phen$c0
      curK[i,"Val"] <- c(t(u)%*%(kinship[[v]]$Q)%*%u)
    }else{
      if(curK$Type[i]=="kinship"){
        ra  <- sum(ageClass$rcont0[ageClass$Breed==b & ageClass$age==1])
        r0  <- ageClass$rcont0[ageClass$Breed==b]
        nr  <- ageClass$Class[ageClass$Breed==b]
        NPO <- pPOpop(ageClass, Breed=b, N0=N*ra)[nr, nr, drop=FALSE]
        diag(NPO)    <- 1/pmax(r0*N,1)
        mkin         <- kinship[[v]]$mkin
        curK[i,"Val"] <- c(r0%*%(mkin$f + NPO*mkin$diffF)%*%r0)
      }
      
      if(curK$Type[i]=="nat. kin."){
        ra  <- sum(ageClass$rcont0[ageClass$Breed==b & ageClass$age==1])
        r0  <- ageClass$rcont0[ageClass$Breed==b]
        nr  <- ageClass$Class[ageClass$Breed==b]
        NPO <- pPOpop(ageClass, Breed=b, N0=N*ra)[nr, nr, drop=FALSE]
        diag(NPO) <- 1/pmax(r0*N,1)
        mkin <- kinship[[v]]$mkin1
        x1   <- c(r0%*%(mkin$f + NPO*mkin$diffF)%*%r0)
        mkin <- kinship[[v]]$mkin2
        x2   <- c(r0%*%(mkin$f + NPO*mkin$diffF)%*%r0)
        curK[i,"Val"] <- x1/x2
      }
    }
    curK[i,"obj.fun.and.constraints"] <- paste0(c("min.", "ub."), curK$Var[i], collapse=", ")
  }
  
  cur <- rbind(curT, curK)
  
  ### Append data frame with current breed contributions ###
  
  if(length(BreedNames)>1){
    curBC <- data.frame(
      Type  = rep("breed cont",length(bc)),
      Var   = paste0("bc.",names(bc)),
      Breed = names(bc),
      Val   = bc,
      obj.fun.and.constraints = '',
      stringsAsFactors=FALSE)
    
    cur <- rbind(cur, curBC)
  }
  
  cur$Name <- str_replace(cur$Var, paste0(paste0(".", BreedNames),collapse="|"),"")
  cur$Name <- str_replace(cur$Name, "bc", "Breed Contribution")
  rownames(cur) <- 1:nrow(cur)
  return(cur)
}