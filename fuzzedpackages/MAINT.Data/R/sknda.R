Ilocsnda <- function(Conf,x,prior,W,B,mu,gamma1,alpha,ksi,egvtol,limlnk2)
{
  p <- 2*x@NIVar
  nk <- as.numeric(table(x@grouping))
  N <- sum(nk)
  k <- length(nk) 
  if (k==1) {
    stop("The data belongs to one single group. A partition into at least two different groups is required\n")
  }

  if (prior[1]=="proportions") { prior <- nk/N }
  names(prior) <- levels(x@grouping)
  Wi <- Safepdsolve(W,maxlnk2=limlnk2,scale=TRUE)
  if (is.null(Wi)) {
    warning("Ilocsnda function received a singular matrix in the  W argument\n")
    return(NULL)
  }
  WiBdecp <- eigen(Wi%*%B)
  if (Conf==1) {
    r <- min(p,k-1)
  }  else  {
    WiBegval <- Re(WiBdecp$values)	
    posWiBegval <- WiBegval[WiBegval>egvtol]
    r <- length(posWiBegval)
  }
  eigvct <- Re(WiBdecp$vectors[,1:r])
  if (r==1) { dim(eigvct) <- c(p,1) }
  sclvar <- apply(eigvct,2,function(v) v%*%W%*%v)	
  if (r>1) {	
    scaling <- scale(eigvct,center=FALSE,scale=sqrt(sclvar))
    dimnames(scaling) <- list(rownames(W),paste("LD",1:r,sep=""))	
    attr(scaling,"scaled:scale") <- NULL
  } else {	
    scaling <- matrix(eigvct/sqrt(sclvar),ncol=1,dimnames=list(rownames(W),"LD1"))
  }
  eta <- alpha / sqrt(diag(W))
  names(eta) <- rownames(W)

  new("IdtSNlocda",prior=prior,ksi=ksi,eta=eta,scaling=scaling,mu=mu,gamma1=gamma1,N=N,CovCase=Conf) 
}

setMethod("snda",
  signature(x = "IdtLocSNMANOVA"),
  function(x, prior="proportions", selmodel=BestModel(H1res(x)), egvtol=1.0e-10, silent=FALSE, k2max=1e6, ...)
  {
    limlnk2 <- log(k2max)
    if (is.character(selmodel)) { 
      selmodel <- sapply(selmodel,function(mod) which(mod==x@H0res@ModelNames))
    }
    if (!is.finite(x@H0res@logLiks[selmodel]) || !is.finite(x@H1res@logLiks[selmodel]))
    {
      if (silent) { 
        return(NULL)
      } else {
        stop("Trying to compute discriminant functions from a model with non-finite log-likelihood\n")
      }
    }
    W <- x@H1res@CovConfCases[[selmodel]]$OmegaE	

    Ilocsnda(Conf=selmodel,x=x,prior=prior,W=W,B=x@H0res@CovConfCases[[selmodel]]$OmegaE-W,
      mu=x@H1res@CovConfCases[[selmodel]]$muE,gamma1=x@H1res@CovConfCases[[selmodel]]$gamma1E,
      alpha=x@H1res@CovConfCases[[selmodel]]$alphaE,ksi=x@H1res@CovConfCases[[selmodel]]$ksiE,egvtol=egvtol,limlnk2=limlnk2)
  }
)

setMethod("snda",
  signature(x = "IdtLocNSNMANOVA"),
  function(x, prior="proportions", selmodel=BestModel(H1res(x)@SNMod), egvtol=1.0e-10, silent=FALSE, k2max=1e6, ...)
  {
    limlnk2 <- log(k2max)
    H0res <- H0res(x)@SNMod
    H1res <- H1res(x)@SNMod
    if (is.character(selmodel)) {
      selmodel <- sapply(selmodel,function(mod) which(mod==H1res@ModelNames))
    }
    if (!is.finite(H0res@logLiks[selmodel]) || !is.finite(H1res@logLiks[selmodel]))
    {
      if (silent) {
        return(NULL)
      }  else  {
        stop("Trying to compute discriminant functions from a model with non-finite log-likelihood\n")
      }
    }
    W <- H1res@CovConfCases[[selmodel]]$OmegaE	

    Ilocsnda(Conf=selmodel,x=x,prior=prior,W=W,B=H0res@CovConfCases[[selmodel]]$OmegaE-W,
      mu=H1res@CovConfCases[[selmodel]]$muE,gamma1=H1res@CovConfCases[[selmodel]]$gamma1E,
      alpha=H1res@CovConfCases[[selmodel]]$alphaE,ksi=H1res@CovConfCases[[selmodel]]$ksiE,egvtol=egvtol,limlnk2=limlnk2)
  }
)

setMethod("snda",
  signature(x = "IData"),
  function( x, grouping, prior="proportions", CVtol=1.0e-5, subset=1:nrow(x), CovCase=1:4, SelCrit=c("BIC","AIC"),
    Mxt=c("Loc","Gen"), k2max=1e6, ... )
  {
    SelCrit <- match.arg(SelCrit)
    Mxt <- match.arg(Mxt)

#    Config <- getConfig(...)
#    if (is.null(Config))  
#    {
#      Config <- ifelse(CovCase==1,1,CovCase+1)
#      CovCaseArg <- TRUE	
#    } else {  
#      CovCaseArg <- FALSE
#    }	
#    if (x@NIVar==1) {
#      CovCase <- q1CovCase(CovCase) 
#      Config <- q1Config(Config)
#    }  

    if (length(subset) < x@NObs)
    {
      x <- x[subset,]
      grouping <- grouping[subset]
    }
    grouping <- factor(grouping,exclude=NULL)
    grplvls <- levels(grouping)
    if (length(grplvls)==1) {
      stop("The data belongs to one single group. A partition into at least two different groups is required\n")
    }

#    snda(MANOVA(x,grouping,Mxt=Mxt,Model="SKNormal",Config=Config,SelCrit=SelCrit,CVtol=CVtol,k2max=k2max,...),
    snda(MANOVA(x,grouping,Mxt=Mxt,Model="SKNormal",CovCase=CovCase,SelCrit=SelCrit,CVtol=CVtol,k2max=k2max,...),
      prior=prior)
  }
)

setMethod("predict",
  signature(object = "IdtSNlocda"),
  function(object,newdata,prior=object@prior,...)
  {
    if (is(newdata,"IData")) { newdata <- as.matrix(cbind.data.frame(newdata@MidP,newdata@LogR)) }
    if (is(newdata,"data.frame")) { newdata <- as.matrix(newdata) }
    n <- nrow(newdata)
    k <- length(prior) 
    if (k==1) {
      stop("The data belongs to one single group. A partition into at least two different groups is required\n")
    }

    sphdata <- newdata %*% object@scaling 
    sphksi <- object@ksi %*% object@scaling 
    Mahdistovertwo <- apply(sphdata, 1, function(x) apply(sphksi, 1, function(ksi) (sum((x-ksi)^2)/2)))
    skewcorr <- apply(newdata, 1, function(x) apply(object@ksi, 1, function(ksi) pnorm(object@eta%*%(x-ksi))))
    minhlfMD2 <- apply(Mahdistovertwo,2,min)
    wghtdensities <- sweep(exp(sweep(-Mahdistovertwo,2,minhlfMD2,"+"))*skewcorr,1,prior,"*")
    ncnst <- apply(wghtdensities,2,sum)  			# normalizing constants
    posterior <- sweep(wghtdensities,2,STATS=ncnst,FUN="/")
    clres <- apply(posterior, 2, function(pst) return(dimnames(sphksi)[[1]][which.max(pst)]))

    list(class=factor(clres,levels=dimnames(object@ksi)[[1]]),posterior=t(posterior))
  }
)

setMethod("show",					
  signature(object = "IdtSNlocda"),
  function(object)
  {
     cat("Prior probabilities of groups:\n") ; print(object@prior) ; cat("\n")
     cat("Centred parameters:\n")
     cat("Group means:\n") ; print(object@mu) ; cat("\n")
     cat("Skewness coeficients:\n") ; print(object@gamma1) ; cat("\n")
     cat("\nCoefficients of linear discriminants:\n") ; print(object@scaling) ; cat("\n")
     cat("\nDirect parameters:\n")
     cat("Location parameters (ksi):\n") ; print(object@ksi) ; cat("\n")
     cat("Scaled skewness parameters (eta):\n") ; print(object@eta) ; cat("\n")
  }
)

Igensnda <- function(Conf,x,prior,Wg,mu,gamma1,alpha,ksi,limlnk2,silent=FALSE)
{
  p <- 2*x@NIVar
  nk <- as.numeric(table(x@grouping))
  N <- sum(nk)
  k <- length(nk) 
  vnames <- colnames(ksi)
  lev <- levels(x@grouping)
  scaling <- array(dim=c(p,p,k),dimnames=list(vnames,paste("LD",1:p,sep=""),lev))
  ldet <- numeric(k)
  eta <- matrix(nrow=k,ncol=p,dimnames=list(lev,vnames))
  if (prior[1]=="proportions") { prior <- nk/N }
  names(prior) <- lev
  for (g in 1:k)
  {
    if (any(!is.finite(Wg))) scalingg <- NULL
    else scalingg <- Safepdsolve(Wg[,,g],maxlnk2=limlnk2,scale=TRUE)
    if (is.null(scalingg)) {
      warning("Found a non positive definite within groups scatter matrix\n")
      return(NULL) 
    }    
    scaling[,,g] <- scalingg
    ldet[g] <- as.numeric(determinant(Wg[,,g],logarithm=TRUE)$modulus)/2
    eta[g,] <- alpha[g,] / sqrt(diag(Wg[,,g]))
  }
  
  new("IdtSNgenda",prior=prior,ksi=ksi,eta=eta,scaling=scaling,ldet=ldet,lev=lev,mu=mu,gamma1=gamma1,CovCase=Conf) 
}

setMethod("snda",
  signature(x = "IdtGenSNMANOVA"),
  function(x, prior="proportions", selmodel=BestModel(H1res(x)), silent=FALSE, k2max=1e6, ...)
  {
    limlnk2 <- log(k2max)
    if (is.character(selmodel))  { selmodel <- sapply(selmodel,function(mod) which(mod==x@H1res@ModelNames)) }
    if (!is.finite(x@H0res@logLiks[selmodel]) || !is.finite(x@H1res@logLiks[selmodel]))
    {
      if (silent) { 
         return(NULL)
      }  else  {
        stop("Trying to compute discriminant functions from a model with non-finite log-likelihood\n")
      }
    }	
	
    Igensnda(Conf=selmodel,x=x,prior=prior,Wg=x@H1res@CovConfCases[[selmodel]]$OmegaE,
      mu=x@H1res@CovConfCases[[selmodel]]$muE,gamma1=x@H1res@CovConfCases[[selmodel]]$gamma1E,
      alpha=x@H1res@CovConfCases[[selmodel]]$alphaE,ksi=x@H1res@CovConfCases[[selmodel]]$ksiE,limlnk2=limlnk2)
  }
)

setMethod("snda",
  signature(x = "IdtGenNSNMANOVA"),
  function(x, prior="proportions", selmodel=BestModel(H1res(x)@SNMod), silent=FALSE, k2max=1e6, ...)
  {
    limlnk2 <- log(k2max)
    H0res <- H0res(x)@SNMod
    H1res <- H1res(x)@SNMod
    if (is.character(selmodel))  { selmodel <- sapply(selmodel,function(mod) which(mod==H1res@ModelNames)) }
    if (!is.finite(H0res@logLiks[selmodel]) || !is.finite(H1res@logLiks[selmodel]))
    {
      if (silent) {
        return(NULL)
      }  else  {
        stop("Trying to compute discriminant functions from a model with non-finite log-likelihood\n")
      }
    }

    Igensnda(Conf=selmodel,x=x,prior=prior,Wg=H1res@CovConfCases[[selmodel]]$OmegaE,
      mu=H1res@CovConfCases[[selmodel]]$muE,gamma1=H1res@CovConfCases[[selmodel]]$gamma1E,
      alpha=H1res@CovConfCases[[selmodel]]$alphaE,ksi=H1res@CovConfCases[[selmodel]]$ksiE,limlnk2=limlnk2)
  }
)

setMethod("predict",
  signature(object = "IdtSNgenda"),
  function(object,newdata,prior=object@prior,...)
  {
    if (is(newdata,"IData")) { newdata <- as.matrix(cbind.data.frame(newdata@MidP,newdata@LogR)) }
    if (is(newdata,"data.frame")) { newdata <- as.matrix(newdata) }
    n <- nrow(newdata)
    k <- length(prior)
    grpnames <- dimnames(object@ksi)[[1]]
        
    Mahdistovertwo <- matrix(nrow=k,ncol=n,dimnames=list(grpnames,rownames(newdata)))
    skewcorr <- matrix(nrow=k,ncol=n,dimnames=list(grpnames,rownames(newdata)))
    for (g in 1:k)
    {
      sphdata <- newdata %*% object@scaling[,,g] 
      sphksig <- object@ksi[g,] %*% object@scaling[,,g]
      Mahdistovertwo[g,] <- apply( sphdata, 1, function(x) sum((x-sphksig)^2)/2 )
      skewcorr[g,] <- apply( newdata, 1, function(x) pnorm(object@eta[g,]%*%(x-object@ksi[g,])) )
    } 
    minhlfMD2 <- apply(Mahdistovertwo,2,min)
    nrmhlfMD2 <- sweep(Mahdistovertwo,2,minhlfMD2)
    wghtdensities <- sweep(exp(sweep(-nrmhlfMD2,1,object@ldet))*skewcorr,1,prior,"*")
    ncnst <- apply(wghtdensities,2,sum)  			# normalizing constants
    posterior <- sweep(wghtdensities,2,STATS=ncnst,FUN="/")
    clres <- apply(posterior, 2, function(pst) return(grpnames[which.max(pst)]))
    
    list(class=factor(clres,levels=grpnames),posterior=t(posterior))
  }
) 

setMethod("show",					
  signature(object = "IdtSNgenda"),
  function(object)
  {
    cat("Prior probabilities of groups:\n") ; print(object@prior) ; cat("\n")
    cat("Centred parameters:\n")
    cat("Group means:\n") ; print(object@mu) ; cat("\n")
    cat("Group skewness coeficients:\n") ; print(object@gamma1) ; cat("\n")
    cat("\nDirect parameters:\n")
    cat("Location parameters (ksi):\n") ; print(object@ksi) ; cat("\n")
    cat("Scaled skewness parameters (eta):\n") ; print(object@eta) ; cat("\n")
  }
)

setMethod("CovCase",					
  signature(object = "IdtSNlocda"),
  function(object)
  {
    object@CovCase
  }
)

setMethod("CovCase",					
  signature(object = "IdtSNgenda"),
  function(object)
  {
    object@CovCase
  }
)
