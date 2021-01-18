
"fun4KinatN"<-function(obj, ageClass, phen, N){
  if(obj$breed=="across breeds"){
    obj$d1 <- 0
    obj$d2 <- 0
    obj$a1 <- 0*obj$a1
    obj$a2 <- 0*obj$a2
    return(obj)
    }
  if(is.na(N)){stop("Parameter N is NA.\n")}

  ageClass  <- ageClass[ageClass$Breed==obj$breed,]
  phen      <- phen[obj$id,]
  phen$NC   <- obj$NC
  
  classes   <- ageClass$Class
  r0        <- ageClass$rcont0
  r1        <- ageClass$rcont1
  NIndiv0   <- pmax(N*r0,1)
  NIndiv1   <- pmax(N*r1,1)
  nIndiv    <- pmax(ageClass$n,1)
  ra        <- 1 - sum(r1)

  Delta00   <- ra^2*sum((phen$c0*(phen$NC-diag(obj$Q1)))[phen$Breed==obj$breed])/(2*ra*N)
  
  if(ra>0.9999){
    obj$d1    <- Delta00
    obj$a1    <- setNames(rep(0,nrow(phen)), phen$Indiv)
  }else{
    diffF1    <- obj$mkin1$diffF
    NPO       <- pPOpop(ageClass, Breed=obj$breed, N0=ra*N)
    nPO       <- NPOdata(phen, Classes=classes, symmetric=TRUE, quiet=TRUE)
    NPO       <- NPO[classes,classes, drop=FALSE]
    diag(NPO) <- 1/NIndiv0
    diag(nPO) <- 1/nIndiv
    
    u <- ageClass[phen$Class,"rcont1"]
    u <- u*(1/ageClass[phen$Class,"n"] - 1/pmax(N*ageClass[phen$Class,"rcont0"],1))
    u <- u*(obj$mkin1$selfkin[phen$Indiv]-obj$mkin1$kinwac[phen$Indiv])
    u[is.na(u)|phen$Breed!=obj$breed] <- 0
    
    Delta01   <- 2*sum(phen$c1*u)
    Delta11   <- c(r1%*%(diag(1/NIndiv1-1/NIndiv0, nrow=length(NIndiv1))*diffF1 - (nPO-NPO)*diffF1)%*%r1)
    obj$d1    <- Delta00 + Delta01 + Delta11
    obj$a1    <- setNames(-2*c(u), phen$Indiv)
    }
  
  Delta00   <- ra^2*sum((phen$c0*(phen$NC-diag(obj$Q2)))[phen$Breed==obj$breed])/(2*ra*N)

  if(ra>0.9999){
    obj$d2    <- Delta00
    obj$a2    <- setNames(rep(0,nrow(phen)), phen$Indiv)
  }else{
    diffF2 <- obj$mkin2$diffF
    u <- ageClass[phen$Class,"rcont1"]
    u <- u*(1/ageClass[phen$Class,"n"] - 1/pmax(N*ageClass[phen$Class,"rcont0"],1))
    u <- u*(obj$mkin2$selfkin[phen$Indiv]-obj$mkin2$kinwac[phen$Indiv])
    u[is.na(u)|phen$Breed!=obj$breed] <- 0
    
    Delta01   <- 2*sum(phen$c1*u)
    Delta11   <- c(r1%*%(diag(1/NIndiv1-1/NIndiv0, nrow=length(NIndiv1))*diffF2 - (nPO-NPO)*diffF2)%*%r1)
    obj$d2    <- Delta00 + Delta01 + Delta11
    obj$a2    <- setNames(-2*c(u), phen$Indiv)
  }

  obj
}