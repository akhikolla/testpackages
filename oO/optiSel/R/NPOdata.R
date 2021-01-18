
"NPOdata" <- function(phen, Classes, proportion=TRUE, symmetric=FALSE, quiet=FALSE){
  phen <- phen[phen$Class %in% Classes, ]
  
  withoutSexes <-  any(is.na(phen$Sex))
  if(withoutSexes){ #Not yet implemented. Function returns Zero-Matrix
    nPO     <- matrix(0, nrow=length(Classes), ncol=length(Classes))
    rownames(nPO) <- Classes
    colnames(nPO) <- Classes
    return(nPO)
  }
  
  phen$Sire.Class <- phen[phen$Sire,"Class"]
  phen$Dam.Class  <- phen[phen$Dam, "Class"]
  
  TabS <- table(phen$Class, phen$Sire.Class)
  TabD <- table(phen$Class, phen$Dam.Class)
  TabS <- TabS[rownames(TabS) %in% Classes, colnames(TabS) %in% Classes]
  TabD <- TabD[rownames(TabD) %in% Classes, colnames(TabD) %in% Classes]

  nPO <- matrix(0,nrow=length(Classes), ncol=length(Classes))
  rownames(nPO) <- Classes
  colnames(nPO) <- Classes
  nPO[rownames(TabS), colnames(TabS)] <- TabS
  nPO[rownames(TabD), colnames(TabD)] <- TabD

  if(proportion){
    nvec <- table(phen$Class)[Classes]
    nvec[is.na(nvec)] <- 1
    nPO <- diag(1/nvec) %*% nPO %*% diag(1/nvec)
    rownames(nPO) <- Classes
    colnames(nPO) <- Classes
  }
  
  if(symmetric){
    nPO[is.na(nPO)]<-0
    nPO <- nPO + t(nPO)
    diag(nPO) <- 0
  }
  
  nPO
}