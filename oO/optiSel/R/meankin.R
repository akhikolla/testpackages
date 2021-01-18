
"meankin"<-function(phen, Q, Classes){
  
  if(!is.matrix(Q) || !is.numeric(Q)){stop("Q must be a numeric matrix.\n")}
  if(nrow(Q)!=ncol(Q)){stop("Q must be a quadratic matrix.\n")}

  use <- phen$Class %in% Classes
  if(!all(use)){
    Q    <- Q[use, use]
    phen <- phen[use, ]
  }

  by      <- phen$Class
  selfkin <- setNames(diag(Q), colnames(Q))
  
  diag(Q) <- NA
  
  POkin <- NULL
  if(max(phen$Age)>1){
    Parents <- intersect(c("Sire", "Dam"), colnames(phen))
    for(Par in Parents){
      use       <- phen[[Par]] %in% colnames(Q)
      POpair    <- as.matrix(phen[use, c("Indiv", Par)])
      X         <- data.frame(Class1 = phen[POpair[,"Indiv"],"Class"],
                              Class2 = phen[POpair[,Par],    "Class"],
                              kin    = Q[POpair], 
                              stringsAsFactors = FALSE)
      POkin     <- rbind(POkin, X)
      Q[POpair] <- NA
      Q[POpair[,c(2,1)]] <- NA
    }
  }
  
  Q <- aggregate(Q, by=list(by), FUN=mean, na.rm=TRUE)
  rownames(Q) <- Q$Group.1
  Q$Group.1 <- NULL
  Q <- t(Q)
  
  cnumber <- match(phen$Class, colnames(Q))
  kinwac  <- Q[cbind(1:nrow(Q), cnumber)]
  kinwac  <- setNames(kinwac, phen$Indiv)

  Q <- aggregate(Q, by=list(by), FUN=mean, na.rm=TRUE)
  rownames(Q) <- Q$Group.1
  Q$Group.1 <- NULL
  Q <- as.matrix(Q)

  f <- matrix(NA, nrow=length(Classes), ncol=length(Classes))
  rownames(f) <- Classes
  colnames(f) <- Classes
  f[rownames(Q), colnames(Q)] <- Q
  f[is.na(f)] <- mean(f, na.rm=TRUE)

  if(!is.null(POkin) && nrow(POkin)>0){
    POkin <- aggregate(POkin$kin,list(POkin$Class1, POkin$Class2), mean, na.rm=TRUE)
    diffF                 <- 0*f+mean(POkin$x)
    diffF[as.matrix(POkin[,c(1,2)])] <- POkin$x
    diffF[as.matrix(POkin[,c(2,1)])] <- POkin$x
    diffF                 <- diffF - f
  }else{
    diffF <- mean(selfkin)/2 - f
  }
  
  ddiffF      <- tapply(selfkin, by, mean  )[Classes] - diag(f)
  ddiffF[is.na(ddiffF)] <- mean(ddiffF, na.rm=TRUE)
  diag(diffF) <- ddiffF
 
  return(list(f=f, diffF=diffF, kinwac=kinwac, selfkin=selfkin))
}