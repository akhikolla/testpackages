### This function estimates the proportion of parent-offspring pairs
### within each pair of Sex X Age classes

"pPOpop" <- function(ageClass, Breed=NA, N0){
  b         <- Breed
  ageClass  <- ageClass[ageClass$Breed==b,]
  withSexes <- all(ageClass$Group %in% c("male", "female"))
  
  if(withSexes){
    cont <- data.frame(age    = 1:max(ageClass$age), 
                       female = ageClass$cont0[ageClass$Group=="male"], 
                       male   = ageClass$cont0[ageClass$Group=="female"])
  
    pOff <- data.frame(age    = 1:max(ageClass$age), 
                       female = ageClass$propOff[ageClass$Group=="male"], 
                       male   = ageClass$propOff[ageClass$Group=="female"])
  
    fClasses <- paste0(Breed,".Age", cont$age, ".female")
    mClasses <- paste0(Breed,".Age", cont$age, ".male")
    Classes  <- c(fClasses, mClasses)

    if(nrow(cont)==1){
      NPO <- matrix(0, nrow=2, ncol=2, dimnames=list(Classes, Classes))
      return(NPO)
    }
    
    nAges <- nrow(cont)
    Xm <- matrix(0, nrow=nAges, ncol=nAges)
    Xf <- matrix(0, nrow=nAges, ncol=nAges)
    for(i in 1:(nAges-1)){
      Xm[i, (i+1):nAges] <- pOff$male[1:(nAges-i)]
      Xf[i, (i+1):nAges] <- pOff$female[1:(nAges-i)]
    }
    
    ra        <- cont$female[1] + cont$male[1]  
    N         <- round(N0/ra)
    Nm        <- round(cont$male[1]*N)
    Nf        <- round(cont$female[1]*N)
    Ncohort   <- (cont$male+cont$female)*N
    Nclass    <- c(cont$female*N, cont$male*N)
    psurvived <- c(cont$female*N/Nf, cont$male*N/Nm)
    pmale     <- cont$male/(cont$male+cont$female)
    
    NPO  <- ((Ncohort) %*% t(psurvived)) * cbind(Xf, Xm)
    NPO  <- rbind(diag(1-pmale)%*%NPO, diag(pmale)%*%NPO)
    NPO  <- ((1/Nclass) %*% t(1/Nclass)) * NPO 
    
  }else{

    cohort  <- ageClass$cont0
    pOff    <- ageClass$propOff
    Classes <- paste0(Breed,".Age", 1:max(ageClass$age), ".cohort")
    if(length(cohort)==1){
      NPO <- matrix(0, nrow=1, ncol=1, dimnames=list(Classes, Classes))
      return(NPO)
    }
    nAges <- length(cohort)
    X     <- matrix(0, nrow=nAges, ncol=nAges)
    for(i in 1:(nAges-1)){
      X[i, (i+1):nAges] <- pOff[1:(nAges-i)]
    }
    N         <- round(N0/cohort[1])
    Ncohort   <- cohort*N
    psurvived <- cohort/cohort[1]
    NPO      <- ((Ncohort) %*% t(psurvived)) * X
    NPO      <- ((1/Ncohort) %*% t(1/Ncohort)) * NPO 
  }
  
  NPO[is.na(NPO)] <- 0
  NPO <- NPO + t(NPO)
  rownames(NPO) <- Classes
  colnames(NPO) <- Classes
  return(NPO)
}