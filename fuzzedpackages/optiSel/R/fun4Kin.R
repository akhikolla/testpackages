
### This function defines the linear component (a) and the constant term (d).
### The quadratic term (Q) remains unchanged.

"fun4Kin" <- function(obj, ageClass, phen, N){
  if(obj$breed=="across breeds"){
    obj$d <- 0
    obj$a <- 0*obj$a
    return(obj)
    }
  if(is.na(N)){stop("Parameter N is NA.\n")}
  
  ageClass  <- ageClass[ageClass$Breed==obj$breed,]
  classes   <- ageClass$Class
  r0        <- ageClass$rcont0
  r1        <- ageClass$rcont1
  diffF     <- obj$mkin$diffF
  NIndiv0   <- pmax(N*r0,1)
  NIndiv1   <- pmax(N*r1,1)
  nIndiv    <- pmax(ageClass$n,1)
  ra        <- 1 - sum(r1)

  Delta00  <- ra^2*(1-sum((phen$c0*diag(obj$Q))[phen$Breed==obj$breed]))/(2*ra*N)
  
  if(ra>0.9999){
    obj$a    <- setNames(rep(0,nrow(phen)), phen$Indiv)
    obj$d    <- Delta00
  }else{
    NPO       <- pPOpop(ageClass, Breed=obj$breed, N0=ra*N)
    nPO       <- NPOdata(phen, Classes=classes, symmetric=TRUE, quiet=TRUE)
    NPO       <- NPO[classes,classes, drop=FALSE]
    diag(NPO) <- 1/NIndiv0
    diag(nPO) <- 1/nIndiv
    
    u <- ageClass[phen$Class,"rcont1"]
    u <- u*(1/pmax(ageClass[phen$Class,"n"],1) - 1/pmax(N*ageClass[phen$Class,"rcont0"],1))
    u <- u*(obj$mkin$selfkin[phen$Indiv]-obj$mkin$kinwac[phen$Indiv])
    u[is.na(u)|phen$Breed!=obj$breed] <- 0
    
    Delta01  <- 2*sum(phen$c1*u)
    Delta11  <- c(r1%*%(diag(1/NIndiv1-1/NIndiv0, nrow=length(NIndiv1))*diffF - (nPO-NPO)*diffF)%*%r1)
    obj$a    <- setNames(-2*c(u), phen$Indiv)
    obj$d    <- Delta00 + Delta01 + Delta11
    }
  
  #cat(paste0("Delta00 =", Delta00, "\n"))
  #cat(paste0("Delta01 =", Delta01, "\n"))
  #cat(paste0("Delta11 =", Delta11, "\n\n"))



  
  return(obj)
}