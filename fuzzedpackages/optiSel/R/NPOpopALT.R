### This function estimates the number of parent-offspring pairs
### within each pair of the Sex X Age classes defined in argument cont

"NPOpop" <- function(ageClass, phen=NULL, pOff=NULL, N0=NA, Breed=NA, only.alive=TRUE, by.sex=TRUE, proportion=TRUE, symmetric=FALSE, quiet=TRUE){
  b <- Breed
  AgeSym <- paste0(Breed,".Age")

  cont <- dcast(ageClass[ageClass$Breed==b,], age ~ Group, value.var="cont0")
  setDF(cont)
  
  cont     <- checkcont(cont)
  cont     <- rbind(cont, 0)
  cont$age <- 1:nrow(cont)
  
  ra       <- cont$cohort[1]
  overlapping  <- ra < 0.9999
  withoutSexes <- !is.na(Breed) && !is.null(phen) && any(is.na(phen$Sex[phen$Breed==Breed]))
  Ages     <- paste0(AgeSym, cont$age)
  mClasses <- paste0(Ages, ".male")
  fClasses <- paste0(Ages, ".female")
  
  if(withoutSexes){
    Classes  <- paste0(Ages, ".cohort") #Not yet implemented, function returns Zero-Matrix
  }else{
    Classes  <- c(mClasses, fClasses)
  }
  
  
  if(withoutSexes || !overlapping){
    if(by.sex || withoutSexes){
      NPO <- matrix(0, nrow=length(Classes), ncol=length(Classes))
      rownames(NPO) <- Classes
    }else{
      NPO <- matrix(0, nrow=length(Ages), ncol=length(Classes))
      rownames(NPO) <- Ages
    }
    colnames(NPO) <- Classes
    return(NPO)
  }
  
  if(!is.null(phen)){
      phen <- phen[phen$Breed %in% Breed,]
  }
  
  if(is.na(N0)){
    if(is.null(phen)){
      stop("Argument N0 is NA.\n")
    }else{
      N0 <- max(table(phen$Born))
      if(!quiet){cat(paste0("Each age cohort is assumed to consist of N0=",N0," individuals.\n"))}
    }
  }
  
  
  N    <- round(N0/ra)
  Nm   <- round(cont$male[1]*N)
  Nf   <- round(cont$female[1]*N)
  
  if(is.null(pOff)){
    pOff <- dcast(ageClass[ageClass$Breed==b,], age ~ Group, value.var="propOff")
    setDF(pOff)
  }else{
    if(nrow(pOff)<length(Ages)-1){
      stop("Not enough rows in argument pOff.\n")
    }
  }
  
  
  #### Define contL with columns #####################################################
  ### Class: Name of the AgeXSex Class                                               #
  ### r:     Contribution of each AgeXSex class to the population at time t          #
  ### N:     Expected number of individuals in the population for each AgeXSex class #
  ### n:     Number of individuals in the sample from each AgeXSex class             #
  ####################################################################################
  contL           <- melt(cont[,c("age", "male", "female")], id="age")
  colnames(contL) <- c("Age", "Sex", "r")
  contL$Class     <- paste0(AgeSym, contL$Age, ".", contL$Sex)
  contL$N         <- N*contL$r
  if(is.null(phen)){
    nTab <- data.frame(Class=contL$Class, n=0,stringsAsFactors=FALSE)
  }else{
    nTab            <- melt(table(phen$Class))
    colnames(nTab)  <- c("Class", "n")
  }
  contL           <- merge(contL, nTab,by=c("Class"), all.x=TRUE)
  contL$n[is.na(contL$n)]<-0
  rownames(contL) <- contL$Class
  
  #### Define matrix nPO  #################################
  ### containing the number of parent-offspring pairs   ### 
  ### in the data set for each pair of classes          ###
  #########################################################
  nPO <- matrix(0,nrow=length(Ages), ncol=length(Classes))
  rownames(nPO) <- Ages
  colnames(nPO) <- Classes
  if(!is.null(phen)){
    phen <- phen[phen$Class %in% Classes, ]
    phen$Sire.Class <- phen[phen$Sire,"Class"]
    phen$Dam.Class  <- phen[phen$Dam, "Class"]
    TabS <- table(paste0(AgeSym, phen$Age), phen$Sire.Class)
    TabD <- table(paste0(AgeSym, phen$Age), phen$Dam.Class)
    nPO[rownames(TabS), colnames(TabS)] <- TabS
    nPO[rownames(TabD), colnames(TabD)] <- TabD
  }

  
  ###############################################################
  # alphaM[k1,k2] = Conditional probability that                #
  #         an indivdual was born at time t-k1+1                #
  #         if its sire was born at time  t-k2+1                #
  ###############################################################
  
  alphaM <- matrix(0, nrow=length(Ages), ncol=length(mClasses))
  alphaF <- matrix(0, nrow=length(Ages), ncol=length(fClasses))
  rownames(alphaM) <- Ages
  rownames(alphaF) <- Ages
  colnames(alphaM) <- mClasses
  colnames(alphaF) <- fClasses 
  
  
  for(k1 in 1:(nrow(alphaM)-1)){
    for(k2 in (k1+1):ncol(alphaM)){
      alphaM[k1, k2] <- pOff[k2-k1, "male"]/sum(pOff[1:(k2-1), "male"])
      alphaF[k1, k2] <- pOff[k2-k1, "female"]/sum(pOff[1:(k2-1), "female"])
    }
  }
  
  ##################################################################
  # Xm[k1,k2] = Probability that the sire was born at time t-k2+1  #
  #             if the indivdiual is in age cohort k1              #
  # Pm[k2] =  Probability that the sire of an individual randomly  #
  #           chosen from the population was born at time t-k2+1   #
  ##################################################################
  d  <- nrow(pOff)
  r  <- setNames(contL[mClasses,"r"]+contL[fClasses,"r"], Ages)
  
  ### Get approximation of Xm and Xf ####
  Xm  <- matrix(0, nrow=length(Ages), ncol=length(mClasses))
  Xf  <- matrix(0, nrow=length(Ages), ncol=length(fClasses))
  
  XmSum <- sum(alphaM[1:d,d+1]/r[1:d])
  XfSum <- sum(alphaF[1:d,d+1]/r[1:d])
  for(k1 in 1:d){
    k2 <- (k1+1):(d+1)
    Xm[k1,k2]<- alphaM[d+1-k2+k1,d+1]/r[d+1-k2+k1]/XmSum
    Xf[k1,k2]<- alphaF[d+1-k2+k1,d+1]/r[d+1-k2+k1]/XfSum
  }
  
  ### Get approximation of vectors Pm and Pf
  ### containing the probaobilies that the sire (or dam)
  ### of a randomly chosen individual is in the specified age class
  
  Lm  <- c(pOff$male%*%(1:d))
  if(floor(Lm)==ceiling(Lm)){
    Pm <- c(0,rep(0.01,Lm-1), r)[1:(d+1)]
  }else{
    rm1 <- c(0,rep(0.01,floor(Lm)-1),  r)[1:(d+1)]
    rm2 <- c(0,rep(0.01,ceiling(Lm)-1),r)[1:(d+1)]
    Lm1 <- c(rm1%*%(1:length(rm1))) - c(r%*%(1:length(r)))
    Lm2 <- c(rm2%*%(1:length(rm2))) - c(r%*%(1:length(r)))
    Pm  <- ((Lm-Lm2)/(Lm1-Lm2))*rm1+(1-(Lm-Lm2)/(Lm1-Lm2))*rm2
  }
  
  Lf  <- c(pOff$female%*%(1:d))
  if(floor(Lm)==ceiling(Lm)){
    Pf <- c(0,rep(0.01,Lf-1), r)[1:(d+1)]
  }else{
    rf1 <- c(0,rep(0.01,floor(Lf)-1),  r)[1:(d+1)]
    rf2 <- c(0,rep(0.01,ceiling(Lf)-1),r)[1:(d+1)]
    Lf1 <- c(rf1%*%(1:length(rf1))) - c(r%*%(1:length(r)))
    Lf2 <- c(rf2%*%(1:length(rf2))) - c(r%*%(1:length(r)))
    Pf  <- ((Lf-Lf2)/(Lf1-Lf2))*rf1+(1-(Lf-Lf2)/(Lf1-Lf2))*rf2
  }
  
  ### Improve approximations for pM, Pf, Xm and Xf iteratively
  am   <- apply(Xm,1,sum)
  af   <- apply(Xf,1,sum)
  
  GammaM   <- matrix(0, nrow=length(Ages), ncol=length(mClasses))
  GammaF   <- matrix(0, nrow=length(Ages), ncol=length(fClasses))
  rownames(GammaM) <- Ages
  rownames(GammaF) <- Ages
  colnames(GammaM) <- mClasses
  colnames(GammaF) <- fClasses
  for(k1 in 1:d){
    for(k2 in (k1+1):(d+1)){
      GammaM[k1,k2]<- alphaM[k1,k2]/r[k1]
      GammaF[k1,k2]<- alphaF[k1,k2]/r[k1]
    }
  }  
  
  Xm1 <- Xm
  for(j in 1:5000){
    Xm <- GammaM %*% diag(Pm)
    Dm <- am/apply(Xm,1,sum)
    Xm <- diag(Dm) %*% Xm 
    Xm[is.na(Xm)] <- 0
    #print(paste0(j, "  ", max(abs(Xm-Xm1))))
    if(max(abs(Xm-Xm1))<0.0001)break
    Xm1 <- Xm
    Pm  <- c(r%*%Xm)
  }
  rownames(Xm) <- Ages
  colnames(Xm) <- mClasses
  names(Pm)    <- mClasses
  
  Xf1 <- Xf
  for(j in 1:5000){
    Xf <- GammaF %*% diag(Pf)
    Df <- af/apply(Xf,1,sum)
    Xf <- diag(Df) %*% Xf 
    Xf[is.na(Xf)] <- 0
    if(max(abs(Xf-Xf1))<0.0001)break
    Xf1 <- Xf
    Pf  <- c(r%*%Xf)
  }
  rownames(Xf) <- Ages
  colnames(Xf) <- fClasses
  names(Pf)    <- fClasses
  
  ###### Create Matrix NPO ################################
  ### Containing the estimated numbers of parent-offspring#
  ###  pairs in the idealized population                  #
  #########################################################
  nvec <- tapply(contL$n,contL$Age,sum)
  Nvec <- tapply(contL$N,contL$Age,sum)
  names(nvec) <- paste0(AgeSym, names(nvec))
  names(Nvec) <- paste0(AgeSym, names(Nvec))
  
  pr  <- ifelse(is.na(N0/nvec), 1, N0/nvec)
  NPO <- diag(pmin(pr,1))%*%nPO + diag(pmax(N0-nvec,0))%*%cbind(Xm,Xf)
  
  if(only.alive){
    rs  <- contL[Classes,"N"]/ifelse(contL[Classes,"Sex"]=="male", Nm, Nf)
    NPO <- diag(Nvec/N0) %*% NPO %*% diag(rs)
    rownames(NPO) <- Ages
    colnames(NPO) <- Classes
  }
 
  if(by.sex){
    prop        <- cont
    prop$male   <- cont$male/(cont$male+cont$female)
    prop$female <- cont$female/(cont$male+cont$female)
    prop        <- melt(prop, id="age")
    prop        <- setNames(prop[[3]], paste0(AgeSym, prop[[1]], ".", prop[[2]]))
    prop[is.na(prop)] <- 0.5
    NPO <- rbind(diag(prop[mClasses])%*%NPO, diag(prop[fClasses])%*%NPO)
    rownames(NPO) <- Classes
    colnames(NPO) <- Classes
  }

  if(proportion){
    if(only.alive & by.sex){
      NvecL <- setNames(contL[Classes,"N"], Classes)
      NvecR <- setNames(contL[Classes,"N"], Classes)
    }
    if(only.alive & !by.sex){
      NvecL <- Nvec
      NvecR <- setNames(contL[Classes,"N"], Classes)
    }
    if(!only.alive & by.sex){
      NvecL <- setNames(ifelse(contL[Classes,"Sex"]=="male",Nm, Nf), Classes)
      NvecR <- setNames(ifelse(contL[Classes,"Sex"]=="male",Nm, Nf), Classes)
    }
    if(!only.alive & !by.sex){
      NvecL <- setNames(rep(N0, length(Ages)), Ages)
      NvecR <- setNames(ifelse(contL[Classes,"Sex"]=="male",Nm, Nf), Classes)
    }
    
    NPO <- diag(1/NvecL) %*% NPO %*% diag(1/NvecR)
    NPO[is.na(NPO)] <- 0
    rownames(NPO) <- names(NvecL)
    colnames(NPO) <- Classes
  }
  
  if(symmetric && by.sex){
    NPO[is.na(NPO)]<-0
    NPO <- NPO + t(NPO)
    diag(NPO) <- 0
  }
    

  
  NPO
}