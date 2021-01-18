#' @title Tuning funtion of splmm object
#'
#' @description This package fits linear mixed-effects models for high-dimensional data (p>>n) with penalty for both the fixed effects and random effects for variable selection.
#'
#' @import penalized
#' 
#' @import methods
#' 
#' @import emulator 
#' 
#' @import miscTools
#' 
#' @import penalized
#' 
#' @param x
#'
#' @param y
#'
#' @param z
#'
#' @param grp
#'
#' @param lam1
#'
#' @param lam2
#'
#' @param nonpen.b
#'
#' @param nonpen.L
#'
#' @param penalty.b
#'
#' @param penalty.L
#'
#' @param CovOpt
#'
#' @param standardize
#'
#' @param control
#'
#' @return splmm
#'
#' @examples
#'
#' @export

splmm <- function(x,y,z,grp,lam1,lam2,nonpen.b=1,nonpen.L=1,penalty.b=c("lasso","scad"),
                  penalty.L=c("lasso","scad"),CovOpt=c("nlminb","optimize"),standardize=TRUE,control=splmmControl())
  UseMethod("splmm")

splmm.default <- function(x,y,z,grp,lam1,lam2,nonpen.b=1,nonpen.L=1,penalty.b=c("lasso","scad"),
                          penalty.L=c("lasso","scad"),CovOpt=c("nlminb","optimize"),standardize=TRUE,control=splmmControl()){
  
  # --- Introductory checks and reformulations ---
  # ----------------------------------------------
  
  CovOpt <- match.arg(CovOpt)
  penalty.b <- match.arg(penalty.b)
  penalty.L <- match.arg(penalty.L)
  
  # transform data
  if(is.data.frame(x)) x = as.matrix(x)
  if(is.data.frame(z)) z = as.matrix(z)
  
  # do some checks
  if (!is.matrix(x)) stop("x has to be a matrix or data frame")
  if (any(is.na(x))) stop("Missing values in x not allowed")
  
  if (any(is.na(y))) stop("Missing values in y not allowed")
  if (!is.numeric(y)) stop("y has to be of type 'numeric'")
  if (nrow(x)!=length(y)) stop("x and y have not correct dimensions")
  
  if (!is.matrix(z)) stop("z has to be a matrix")
  if (any(is.na(z))) stop("Missing values in z not allowed")

  if (any(x[,1]!=rep(1,dim(x)[[1]]))) stop("first column is not the intercept")
  
  if (length(levels(grp))==1) stop("Only one group. No covariance parameters!")
  
  if (!all(nonpen.b%in%1:dim(x)[[2]])) stop("Error with the argument nonpen for beta")
  if (!all(nonpen.L%in%1:dim(z)[[2]])) stop("Error with the argument nonpen for L")
  
  if (length(lam1)!=1|length(lam2)!=1) stop("lambda can only be one number")
  if (lam1<0|lam2<0) stop("regularization parameter must be positive")
  
  
  ##### Standardize covariates
  
  if (standardize)
  {
    xOr <- x
    meanx <- apply(x[,-1,drop = FALSE],2,mean)
    sdx <- apply(x[,-1,drop = FALSE],2,sd)
    x <- cbind(1,scale(x[,-1],center=meanx,scale=sdx))
    
    zOr <- z
    meanz <- apply(z[,-1,drop = FALSE],2,mean)
    sdz <- apply(z[,-1,drop = FALSE],2,sd)
    z <- cbind(1,scale(z[,-1,drop = FALSE],center=meanz,scale=sdz))
  }
  
  ##### allocate variables
  
  grp <- factor(grp)
  N <- length(levels(grp)) # N is the number of groups
  p <- dim(x)[[2]]         # p is the number of covariates
  q <- dim(z)[[2]]         # q is the number of random effects variables
  ntot <- length(y)        # ntot is the total number of observations
  Q <- q*(q+1)/2           # maximum number of variance components parameters
  
  lambda1 = lam1*N
  lambda2 = lam2*N
  ###### save the grouped information as components of a list
  yGrp <- split(y,grp)
  xGrp <- split.data.frame(x,grp)
  zGrp <- split.data.frame(z,grp)
  zIdGrp <- mapply(ZIdentity,zGrp)
  
  
  ##### Initiation
  
  ###### beta
  init <- optL1(y,x[,-1],model="linear",fold=10,trace=FALSE)
  betaStart <- c(init$fullfit@unpenalized,init$fullfit@penalized)
  
  ###### covariance
  covStart <- covStartingValues(xGrp,yGrp,zGrp,zIdGrp,betaStart,ntot,N)
  sigmaStart <- covStart$sigma
  parsStart <- vecli(covStart$tau*diag(q))
  
  LStart <- t(triang(parsStart,q))
  DStart <- crossprod(LStart)
  
  # --- Calculate objective function for the starting values ---
  # ------------------------------------------------------------
  
  VInvGrp <- mapply(VInv,Z=zGrp,ZId=zIdGrp,MoreArgs=list(D=DStart,sigma=sigmaStart))
  ll1 <- 1/2*ntot*log(2*pi)
  fctStart <- ObjFunction(xGroup=xGrp,yGroup=yGrp,LGroup=VInvGrp,b_nonpen=betaStart[-nonpen.b],L_nonpen=LStart[-nonpen.L,,drop=FALSE],
                          lambda1=lambda1,penalty_b=penalty.b,lambda2=lambda2,penalty_L=penalty.L,ll1=ll1)
  
  
  # --- Coordinate Gradient Descent-iteration ---
  # ---------------------------------------------
  
  # some necessary allocations:
  betaIter <- betaStart
  sigmaIter <- sigmaStart
  LIter <- LStart
  DIter <- DStart
  
  LvecIter <- LStart[lower.tri(LStart,diag = TRUE)]
  convPar <- crossprod(betaIter)
  convCov <- crossprod(c(sigmaStart,LvecIter))
  fctIter <- convFct <- fctStart
  hessian0 <- rep(0,p)
  mat0 <- matrix(0,ncol=p,nrow=N)
  
  
  ##### algorithm parameters
  
  stopped <- FALSE
  doAll <- FALSE
  converged <- 0
  counterIn <- 0
  counter <- 0       # counts the number of outer iterations
  
  
  while ((counter<control$maxIter)&((convPar>control$tol|convFct>control$tol|convCov>control$tol|!doAll ))) {
    
    
    counter <- counter + 1 
    
    
    
    betaIterOld <- betaIter
    LIterOld <- LIter
    fctIterOld <- fctIter
    covIterOld <- c(sigmaIter,LvecIter)
    
    # --- optimization w.r.t the fixed effects vector beta ---
    # --------------------------------------------------------
    
    
    activeSet <- which(betaIter!=0)
    if ((length(activeSet)>min(p,ntot))&(lambda1>0)&(counter>2)) {stopped <- TRUE ; break}
    
    if (counterIn==0 | counterIn>control$number)
    {
      doAll <- TRUE
      activeSet <- 1:p
      counterIn <- 1    
    } else
    {
      doAll <- FALSE
      counterIn <- counterIn+1
    }
    
    
    HessIter <- HessIterTrunc <- HessianMatrix(xGroup=xGrp,LGroup=VInvGrp,activeSet=activeSet,N=N,hessian=hessian0,mat=mat0[,activeSet,drop=FALSE])
    HessIter[activeSet] <- pmin(pmax(HessIter[activeSet],control$lower),control$upper)
    
    LxGrp <- as1(xGrp,VInvGrp,activeSet,N=N)
    
    ll2 <- nlogdet(LGroup=VInvGrp)
    
    for (j in activeSet)
    {
      cut1 <- as2(x=x,y=y,b=betaIter,j=j,activeSet=activeSet,group=grp,sGroup=LxGrp)
      JinNonpen <- j%in%nonpen.b
      
      # optimum can be calculated analytically
      if (HessIterTrunc[j]==HessIter[j])
      {
        if (JinNonpen) {betaIter[j] <- cut1/HessIter[j]} else {
          if(penalty.b=="scad"){
            scada = 3.7
            betaIter[j] <- SoftThreshold(cut1,lambda1)/(HessIter[j]*(1-1/scada))
            
          }else if(penalty.b=="lasso"){
            betaIter[j] <- SoftThreshold(cut1,lambda1)/HessIter[j]
          }
          
        }
      }else
        
        # optimimum is determined by the armijo rule
      {
        armijo <- ArmijoRule_b(xGroup=xGrp,yGroup=yGrp,LGroup=VInvGrp,b=betaIter,j=j,cut=cut1,HkOldJ=HessIterTrunc[j],
                                   HkJ=HessIter[j],JinNonpen=JinNonpen,lambda=lambda1,nonpen=nonpen.b,penalty=penalty.b,
                                   ll1=ll1,ll2=ll2,converged=converged,control=control)
        
        
        betaIter <- armijo$b
        converged <- armijo$converged
        fctIter <- armijo$fct
      }
      
    } 
    
    
    # --- optimization w.r.t the variance components parameters ---
    # -------------------------------------------------------------
    
    # calculations before the covariance optimization
    activeset <- which(betaIter!=0)
    resGrp <- ResAsSplit(x=x,y=y,b=betaIter,f=grp,activeset=activeset)
    ll4 <- lambda1*sum(abs(betaIter[-nonpen.b]))
    
    # optimization of L
    
    activeSet.L = which(rowSums(abs(LIter))!=0)
    
    # calculate the hessian matrices for k in the activeSet
    
    D.grad = D_Gradient(xGroup=xGrp,zGroup=zGrp,LGroup=VInvGrp,yGroup=yGrp,b=betaIter,N=N,q=q)
    L.grad = t(LIter%*%(D.grad+t(D.grad)))
    
    D.hessian = D_HessianMatrix(xGroup=xGrp,zGroup=zGrp,LGroup=VInvGrp,yGroup=yGrp,b=betaIter,N=N,q=q)
    D.hessian.submatrix = matsplitter(D.hessian,q)
    
    for (k in 1:q) {
      
      #L.hessian.sub = sapply(D.hessian.submatrix[k,],function(x) x[,k])
      L.hessian.sub = sapply(D.hessian.submatrix[((k-1)*q+1):(k*q)],function(x) x[,k])
      #L.hessian = diag(2*D.grad[k,k],q,q)+2*LIter%*%(matrix(unlist(D.hessian.submatrix[k,k]),ncol = q,byrow = TRUE)+L.hessian.sub)%*%t(LIter)
      L.hessian = diag(2*D.grad[k,k],q,q)+2*LIter%*%(matrix(unlist(D.hessian.submatrix[[(k-1)*q+k]]),ncol = q,byrow = TRUE)+L.hessian.sub)%*%t(LIter)
      
      for (l in intersect(k:q,activeSet.L)) {
        
        L.lk.grad <- L.grad[l,k]
        L.lk.Hess <- L.hessian[l,l]
        L.lk.Hess <- min(max(L.lk.Hess,control$lower),control$upper)
        
        linNonpen <- l%in%nonpen.L
        
        armijo <- ArmijoRule_L(xGroup=xGrp, yGroup=yGrp, zGroup=zGrp, L=LIter, l=l-1,k=k-1,grad=L.lk.grad,hessian=L.lk.Hess, 
                               b=betaIter,sigma=sigmaIter,zIdGrp=zIdGrp, linNonpen=linNonpen, nonpen=nonpen.L-1, lambda=lambda2, penalty=penalty.L, 
                               ll1=ll1, gamma=control$gamma, maxArmijo=control$maxArmijo, a_init=control$a_init, delta=control$delta, rho=control$rho, converged=converged)
        
        #armijo <- armijoRule_L(xGroup=xGrp,yGroup=yGrp,zGroup=zGrp,L=LIter,l=l,k=k,grad=L.lk.grad,hessian=L.lk.Hess,
        #                       b=betaIter,sigma=sigmaIter,zIdGrp=zIdGrp,linNonpen=linNonpen,lambda=lambda2,nonpen=nonpen.L,penalty=penalty.L,ll1=ll1,converged=converged,control=control)
        
        
        LIter <- armijo$L
        converged <- armijo$converged
        fctIter <- armijo$fct
        
      }
      
      
    }
    
    LIter[LIter < 1e-4] <- 0
    
    DIter = LIter%*%t(LIter)
    LvecIter = LIter[lower.tri(LIter,diag = TRUE)]
    
    
    # optimization of the error variance \sigma^2
    covParOpt <- MLsigma(zGroup=zGrp,zIdGroup=zIdGrp,resGroup=resGrp,q=q,ll1=ll1,ll4=ll4,true.sigma=sigmaIter,D=DIter,
                         trace=control$trace,CovOpt=CovOpt,VarInt=control$VarInt)
    sigmaIter <- covParOpt$sigma
    fctIter <- covParOpt$fct
    
    VInvGrp <- mapply(VInv,Z=zGrp,ZId=zIdGrp,MoreArgs=list(D=DIter,sigma=sigmaIter))
    covIter <- c(sigmaIter,LvecIter)
    
    # --- check convergence ---
    convPar <- sqrt(crossprod(betaIter-betaIterOld))/(1+sqrt(crossprod(betaIter)))
    convFct <- abs((fctIterOld-fctIter)/(1+abs(fctIter)))
    convCov <- sqrt(crossprod(covIter-covIterOld))/(1+sqrt(crossprod(covIter)))
    
    if ((convPar <= control$tol) & (convFct <= control$tol) & (convCov <= control$tol)) counterIn <- 0
    
  }
  
  if (standardize)
  {
    betaIter[-1] <- betaIter[-1]/sdx
    betaIter[1] <- betaIter[1] - sum(meanx*betaIter[-1])
    
    x <- xOr
    xGrp <- split.data.frame(xOr,grp)
    
    z <- zOr
    zGrp <- split.data.frame(zOr,grp)
  }
  
  # --- prediction of the random effects ---
  # ----------------------------------------
  
  biGroup <- u <- list() ; length(biGroup) <- length(u) <- N
  D <- DIter
  #corD <- cov2cor(D)
  
  cholD <- t(triang(LvecIter,q))
  
  for (i in 1:N)
  {
    u[[i]] <- sigmaIter*solve(t(zGrp[[i]]%*%cholD)%*%zGrp[[i]]%*%cholD+sigmaIter^2*diag(q),tol = 1e-40)%*%t(zGrp[[i]]%*%cholD)%*%resGrp[[i]]
    biGroup[[i]] <- 1/sigmaIter*cholD%*%u[[i]]
  }
  
  
  #  --- final calculations ---
  # ---------------------------
  
  # fitted values and residuals
  residGrp <- fittedGrp <- list() ; length(residGrp) <- length(fittedGrp) <- N
  for (i in 1:N)
  {
    fittedGrp[[i]] <- xGrp[[i]][,activeSet,drop=FALSE]%*%betaIter[activeSet,drop=FALSE] + zGrp[[i]]%*%biGroup[[i]]
    residGrp[[i]] <- yGrp[[i]]-fittedGrp[[i]]
  }
  residuals <- unlist(residGrp)
  fitted <- unlist(fittedGrp)
  
  # random effects, sorted per subject
  u <- unlist(u) # corresponds to lmer@u
  bi <- unlist(biGroup) # unsorted random effects
  
  # random effects, sorted per effect
  ranef <- bi[order(rep(1:q,N))] # corresponds to lmer@ranef
  
  # fixed effects without names
  fixef <- betaIter
  names(fixef) <- NULL
  
  # --- summary information ---
  # ---------------------------
  npar <- sum(betaIter!=0) + length(c(sigmaIter,LvecIter))
  logLik <- MLloglik(xGroup=xGrp,yGroup=yGrp,LGroup=VInvGrp,b=betaIter,ntot=ntot,N=N,activeSet=which(betaIter!=0))
  deviance <- -2*logLik
  aic <- -2* logLik + 2*npar
  bic <- -2* logLik + log(ntot)*npar
  
  p <- sum(betaIter!=0)
  q <- sum(diag(D)!=0)
  bbic <- -2*logLik + max(1,log(log(p+q)))*log(ntot)*npar
  ebic <- -2*logLik + (log(ntot)+2*log(p+q))*npar
  
  #if (any(LvecIter==0)) cat("Redundant covariance parameters.","\n")
  
  if (converged>0) cat("Algorithm does not properly converge.","\n")
  if (stopped) {cat("|activeSet|>=min(p,ntot): Increase lambda or set stopSat=FALSE.","\n")
    ; sigmaIter <- LvecIter <- nlogLik <- aic <- bic <- NA ; betaIter <- rep(NA,p) ; bi <- fitted <- residuals <- NULL}
  D[D < 1e-4] <- 0
  out <- list(data=list(x=x,y=y,z=z,grp=grp),coefInit=list(betaStart=betaStart,parsStart=parsStart,sigmaStart=sigmaStart),penalty.b=penalty.b,penalty.L=penalty.L,
              nonpen.b=nonpen.b,nonpen.L=nonpen.L,lambda1=lambda1,lambda2=lambda2,sigma=abs(sigmaIter),Lvec=LvecIter,coefficients=betaIter,D=D,random=bi,u=u,
              ranef=ranef,fixef=fixef,fitted.values=fitted,
              residuals=residuals,converged=converged,logLik=logLik,npar=npar,deviance=deviance,
              aic=aic,bic=bic,bbic=bbic,ebic=ebic,counter=counter,CovOpt=CovOpt,control=control,call=match.call(),
              stopped=stopped,objective=fctIter)
  
  out
  structure(out,class="splmm")
  
  
}