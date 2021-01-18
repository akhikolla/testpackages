#Ilda <- function(Conf,p,nk,prior,means,W,B,egvtol)
Ilda <- function(Conf,p,nk,prior,means,W,B,egvtol,limlnk2)
{
  N <- sum(nk)
  k <- length(nk) 
  if (k==1) {
    stop("The data belongs to one single group. A partition into at least two different groups is required\n")
  }

  if (prior[1]=="proportions") { prior <- nk/N }
  names(prior) <- rownames(means)
#  Wi <- pdwt.solve(W)					
#  Wi <- pdwt.solve(W,silent=TRUE)
  Wi <- Safepdsolve(W,maxlnk2=limlnk2,scale=TRUE)
  if (is.null(Wi)) {
    warning("Ilda function received a singular matrix in the  W argument\n")
    return(NULL)
  }
  WiBdecp <- eigen(Wi%*%B)
  if (Conf!=5) 
  {
    if (Conf==1)  {
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
  }  else  {
    scaling <- diag(1/sqrt(diag(W)))
    dimnames(scaling) <- list(rownames(W),paste("LD",1:p,sep=""))
  }
  new("Idtlda",prior=prior,means=means,scaling=scaling,N=N,CovCase=Conf) 
}

setMethod("lda",
  signature(x = "IdtMxtNDE"),
#  function(x,prior="proportions",selmodel=BestModel(x),egvtol=1.0e-10,silent=FALSE,...)
  function(x,prior="proportions",selmodel=BestModel(x),egvtol=1.0e-10,silent=FALSE,k2max=1e6,...)
  {
    limlnk2 <- log(k2max)
    if (!x@Hmcdt) 
    {
       if (silent)  {
          return(NULL)
       }  else { stop("Trying to compute a linear discriminant function from an estimate of a heteroscedastic mixture\n") }
    }
    if (is.character(selmodel))  { selmodel <- sapply(selmodel,function(mod) which(mod==x@ModelNames)) }
    if (!is.finite(x@logLiks[selmodel]))
    {
       if (silent)  {
          return(NULL)
       }  else { stop("Trying to compute a linear discriminant function from a model with non-finite log-likelihood\n") }
    }
    grouping <- factor(x@grouping,exclude=NULL)
    nk <- as.numeric(table(grouping))
    n <- sum(nk)
    p <- 2*x@NIVar
    k <- length(nk)
    grpmeans <- coef(x)$mu
    glbmeans <- colSums(nk*grpmeans)
    mugdev <- scale(grpmeans,center=glbmeans,scale=FALSE)
    vnames <- unlist(dimnames(grpmeans)[2]) 
    B <- matrix(0.,nrow=p,ncol=p,dimnames=list(vnames,vnames))
    for (g in 1:k) B <- B + (nk[g]/n)*outer(mugdev[g,],mugdev[g,]) 

#    Ilda(Conf=selmodel,p=p,nk=nk,prior=prior,means=grpmeans,W=coef(x,selmodel)$Sigma,B=B,egvtol=egvtol,...)
    Ilda(Conf=selmodel,p=p,nk=nk,prior=prior,means=grpmeans,W=coef(x,selmodel)$Sigma,B=B,egvtol=egvtol,limlnk2=limlnk2,...)
  }
)

setMethod("lda",
  signature(x = "IdtClMANOVA"),
#  function(x,prior="proportions",selmodel=BestModel(H1res(x)),egvtol=1.0e-10,silent=FALSE,...)
  function(x,prior="proportions",selmodel=BestModel(H1res(x)),egvtol=1.0e-10,silent=FALSE,k2max=1e6,...)
  {
    limlnk2 <- log(k2max)
    if (is.character(selmodel))  { selmodel <- sapply(selmodel,function(mod) which(mod==x@H0res@ModelNames)) }
    if (!is.finite(x@H0res@logLiks[selmodel]) || !is.finite(x@H1res@logLiks[selmodel]))
    {
       if (silent)  {
          return(NULL)
       }  else { stop("Trying to compute a linear discriminant function from a model with non-finite log-likelihood\n") }
    }
    W <- coef(H1res(x),selmodel)$Sigma
    Ilda(Conf=selmodel,p=2*x@NIVar,nk=as.numeric(table(x@grouping)),prior=prior,
#      means=coef(H1res(x))$mu,W=W,B=coef(H0res(x),selmodel)$Sigma-W,egvtol=egvtol)
      means=coef(H1res(x))$mu,W=W,B=coef(H0res(x),selmodel)$Sigma-W,egvtol=egvtol,limlnk2=limlnk2)
  }
)

setMethod("lda",
  signature(x = "IdtLocNSNMANOVA"),
#  function(x,prior="proportions",selmodel=BestModel(H1res(x)@NMod),egvtol=1.0e-10,silent=FALSE,...)
  function(x,prior="proportions",selmodel=BestModel(H1res(x)@NMod),egvtol=1.0e-10,silent=FALSE,k2max=1e6,...)
  {
    limlnk2 <- log(k2max)
    H0res <- H0res(x)@NMod
    H1res <- H1res(x)@NMod
    if (is.character(selmodel)) { selmodel <- sapply(selmodel,function(mod) which(mod==H0res@ModelNames)) }
    if ( !is.finite(H0res@logLiks[selmodel]) || !is.finite(H1res@logLiks[selmodel]) )
    {
      if (silent)  {
        return(NULL)
        }  else { stop("Trying to compute a linear discriminant function from a model with non-finite log-likelihood\n") }
     }
     W <- coef(H1res(x),selmodel)$Sigma
     Ilda(Conf=selmodel,p=2*x@NIVar,nk=as.numeric(table(x@grouping)),prior=prior,
#       means=coef(H1res)$mu,W=W,B=coef(H0res(x),selmodel)$Sigma-W,egvtol=egvtol)
       means=coef(H1res)$mu,W=W,B=coef(H0res(x),selmodel)$Sigma-W,egvtol=egvtol,limlnk2=limlnk2)
  }
)

setMethod("lda",
  signature(x = "IData"),
  function(x, grouping, prior="proportions", CVtol=1.0e-5, egvtol=1.0e-10, subset=1:nrow(x), CovCase=1:4,
#    SelCrit=c("BIC","AIC"), silent=FALSE, ...)
    SelCrit=c("BIC","AIC"), silent=FALSE, k2max=1e6, ...)
  {
    limlnk2 <- log(k2max)
    Config <- getConfig(...)
    if (is.null(Config))  
    {
      Config <- ifelse(CovCase==1,1,CovCase+1)
      CovCaseArg <- TRUE	
    } else {  
      CovCaseArg <- FALSE
    }	
    if (x@NIVar==1) {
      CovCase <- q1CovCase(CovCase) 
      Config <- q1Config(Config)
    }  

    SelCrit <- match.arg(SelCrit)
#    if (length(subset) < nrow(x))
    if (length(subset) < x@NObs)
    {
      x <- x[subset,]
      grouping <- grouping[subset]
    }
    grouping <- factor(grouping,exclude=NULL)
    grplvls <- levels(grouping)
    n <- x@NObs
    p <- 2*x@NIVar
    k <- length(grplvls)
    if (k==1)  { 
      errmsg <- "The data belongs to one single group. A partition into at least two different groups is required\n"
      if (silent==FALSE) {
        stop(errmsg)
      } else {
        warning(errmsg)
        return(NULL)
      }
    }
    nk <- as.numeric(table(grouping))
    if (sum(nk)!=n)  { stop("Dimensions of the x and grouping arguments do not agree with each other\n") }
 
    MxtDEst <- IdtNmle(x,grouping,Type="HomMxt",CVtol=CVtol,CovCaseArg=CovCaseArg,Config=Config,SelCrit=SelCrit,...)
#    glbmeans <- colMeans(cbind(x@MidP,x@LogR))
    glbmeans <- colMeans(cbind.data.frame(x@MidP,x@LogR))
    grpmeans <- coef(MxtDEst)$mu
    mugdev <- scale(grpmeans,center=glbmeans,scale=FALSE)
    vnames <- unlist(dimnames(grpmeans)[2]) 
    B <- matrix(0.,nrow=p,ncol=p,dimnames=list(vnames,vnames))
    for (g in 1:k) B <- B + (nk[g]/n) * outer(mugdev[g,],mugdev[g,]) 
    selmodel <- BestModel(MxtDEst)

    Ilda(Conf=selmodel,p=p,nk=nk,prior=prior,means=grpmeans,W=coef(MxtDEst,selmodel)$Sigma,B=B,egvtol=egvtol,limlnk2=limlnk2,...)
  }
)

setMethod("predict",
  signature(object = "Idtlda"),
  function(object,newdata,prior=object@prior,...)
  {
#    if (is(newdata,"IData")) newdata <- as.matrix(cbind(newdata@MidP,newdata@LogR))
    if (is(newdata,"IData")) newdata <- as.matrix(cbind.data.frame(newdata@MidP,newdata@LogR))
    if (is(newdata,"data.frame")) newdata <- as.matrix(newdata)
    n <- nrow(newdata)
    k <- length(prior) 
    if (k==1) {
      stop("The data belongs to one single group. A partition into at least two different groups is required\n")
    }

    sphdata <- newdata %*% object@scaling 
    sphmeans <- object@means %*% object@scaling 
    Mahdistovertwo <- apply(sphdata, 1, function(x) apply(sphmeans, 1, function(mu) (sum((mu-x)^2)/2)))
#    wghtdensities <- sweep(exp(-Mahdistovertwo),1,STATS=prior,FUN="*")
    minhlfMD2 <- apply(Mahdistovertwo,2,min)
    wghtdensities <- sweep(exp(sweep(-Mahdistovertwo,2,minhlfMD2,"+")),1,prior,"*")
    ncnst <- apply(wghtdensities,2,sum)  			# normalizing constants
    posterior <- sweep(wghtdensities,2,STATS=ncnst,FUN="/")
    clres <- apply(posterior, 2, function(pst) return(dimnames(sphmeans)[[1]][which.max(pst)]))
    list(class=factor(clres,levels=dimnames(object@means)[[1]]),posterior=t(posterior))
  }
)

setMethod("show",					
  signature(object = "Idtlda"),
  function(object)
  {
    cat("Prior probabilities of groups:\n") ; print(object@prior) ; cat("\n")
    cat("Group means:\n") ; print(object@means) ; cat("\n")
    cat("Coefficients of linear discriminants:\n") ; print(object@scaling) ; cat("\n")
  }
)

Iqda <- function(Conf,p,nk,lev,prior,means,Wg,limlnk2)
{
  N <- sum(nk)
  k <- length(nk) 
  vnames <- colnames(means)
  scaling <- array(dim=c(p,p,k),dimnames=list(vnames,paste("LD",1:p,sep=""),lev))
  ldet <- numeric(k)

  if (prior[1]=="proportions") prior <- nk/N
  names(prior) <- lev
  Ip <- diag(p)
  for (g in 1:k)
  {
    if (Conf != 5)
    {
#      scaling[,,g] <- backsolve(chol(Wg[,,g]),Ip)
      scalingg <- Safepdsolve(Wg[,,g],maxlnk2=limlnk2,scale=TRUE)
      if (is.null(scalingg)) {
        warning("Found a non positive definite within group matrix\n")
        return(NULL) 
      }    
      scaling[,,g] <- scalingg
      ldet[g] <- as.numeric(determinant(Wg[,,g],logarithm=TRUE)$modulus)/2
    }  else {
      Wd <- diag(Wg[,,g])
      scaling[,,g] <- diag(1/sqrt(Wd))
      ldet[g] <- sum(log(Wd))/2
    }
  }  
  new("Idtqda",prior=prior,means=means,scaling=scaling,ldet=ldet,lev=lev,CovCase=Conf) 
}

setMethod("qda",
  signature(x = "IdtMxtNDE"),
  function(x,prior="proportions",selmodel=BestModel(x),silent=FALSE,k2max=1e6,...)
  {
    limlnk2 <- log(k2max)
    if (x@Hmcdt) 
    {
       if (silent)  {
          return(NULL)
       }  else { stop("Trying to compute a quadratic discriminant function from an estimate of a homoscedastic mixture\n") }
    }
    if (is.character(selmodel))  { selmodel <- sapply(selmodel,function(mod) which(mod==x@ModelNames)) }
    if (!is.finite(x@logLiks[selmodel]))
    {
       if (silent)  {
          return(NULL)
       }  else { stop("Trying to compute a quadratic discriminant function from a model with non-finite log-likelihood\n") }
    }

    Iqda(Conf=selmodel,p=2*x@NIVar,nk=as.numeric(table(x@grouping)),lev=levels(x@grouping),
      prior=prior,means=coef(x)$mu,Wg=coef(x,selmodel)$Sigma,limlnk2=limlnk2)
  }
)

setMethod("qda",
  signature(x = "IdtHetNMANOVA"),
  function(x, prior="proportions", selmodel=BestModel(H1res(x)),silent=FALSE,k2max=1e6, ...)
  {
    limlnk2 <- log(k2max)
    if (is.character(selmodel)) { selmodel <- sapply(selmodel,function(mod) which(mod==x@H1res@ModelNames)) }
    if ( !is.finite(x@H0res@logLiks[selmodel]) || !is.finite(x@H1res@logLiks[selmodel]) )
    {
      if (silent) {
        return(NULL)
      } else { 
        stop("Trying to compute a quadratic discriminant function from a model with non-finite log-likelihood\n")
      }
    }
    Iqda(Conf=selmodel,p=2*x@NIVar,nk=as.numeric(table(x@grouping)),lev=levels(x@grouping),
      prior=prior,means=coef(H1res(x))$mu,Wg=coef(H1res(x),selmodel)$Sigma,limlnk2=limlnk2)
  }
)

setMethod("qda",
  signature(x = "IdtGenNSNMANOVA"),
  function(x, prior="proportions", selmodel=BestModel(H1res(x)@NMod),silent=FALSE,k2max=1e6, ...)
  {
    limlnk2 <- log(k2max)
    H0res <- H0res(x)@NMod
    H1res <- H1res(x)@NMod
    if (is.character(selmodel)) { selmodel <- sapply(selmodel,function(mod) which(mod==H1res@ModelNames)) }
    if (!is.finite(H0res@logLiks[selmodel]) || !is.finite(H1res@logLiks[selmodel]))
    {
      if (silent) {
        return(NULL)
      }  else {
        stop("Trying to compute a quadratic discriminant function from a model with non-finite log-likelihood\n")
      }
    }
    Iqda(Conf=selmodel,p=2*x@NIVar,nk=as.numeric(table(x@grouping)),lev=levels(x@grouping),
      prior=prior,means=coef(H1res)$mu,Wg=coef(H1res,selmodel)$Sigma,limlnk2=limlnk2)
  }
)

setMethod("qda",
  signature(x = "IData"),
  function(x, grouping, prior="proportions", CVtol=1.0e-5, subset=1:nrow(x),
    CovCase=1:4, SelCrit=c("BIC","AIC"), silent=FALSE, k2max=1e6, ...)
  {
    limlnk2 <- log(k2max)
    SelCrit <- match.arg(SelCrit)
    if (x@NIVar==2) CovCase <- q1CovCase(CovCase) 

    Config <- getConfig(...)
    if (is.null(Config))  
    {
      Config <- ifelse(CovCase==1,1,CovCase+1)
      CovCaseArg <- TRUE	
    } else {  
      CovCaseArg <- FALSE
    }	
    if (x@NIVar==1) {
      CovCase <- q1CovCase(CovCase) 
      Config <- q1Config(Config)
    }  

    if (length(subset) < x@NObs)
    {
      x <- x[subset,]
      grouping <- grouping[subset]
    }
    grouping <- factor(grouping,exclude=NULL)
    grplvls <- levels(grouping)
    if (length(grplvls)==1)
    { 
      errmsg <- "The data belongs to one single group. A partition into at least two different groups is required\n"
      if (silent==FALSE) {
        stop(errmsg)
      }  else {
        warning(errmsg)
        return(NULL)
      }
    }
    MxtDEst <- IdtHetMxtNmle(x,grouping,CVtol=CVtol,CovCaseArg=CovCaseArg,Config=Config,SelCrit=SelCrit,...)
    selmodel <- BestModel(MxtDEst)

    Iqda(Conf=selmodel,p=2*x@NIVar,nk=as.numeric(table(grouping)),lev=grplvls,
      prior=prior,means=coef(MxtDEst)$mu,Wg=coef(MxtDEst,selmodel)$Sigma,limlnk2=limlnk2)
  }
)

setMethod("predict",
  signature(object = "Idtqda"),
  function(object,newdata,prior=object@prior,...)
  {
#    if (is(newdata,"IData")) { newdata <- as.matrix(cbind(newdata@MidP,newdata@LogR)) }
    if (is(newdata,"IData")) { newdata <- as.matrix(cbind.data.frame(newdata@MidP,newdata@LogR)) }
    if (is(newdata,"data.frame")) { newdata <- as.matrix(newdata) }
    n <- nrow(newdata)
    p <- ncol(newdata)
    k <- length(prior) 
    grpnames <- dimnames(object@means)[[1]]
    Mahdistovertwo <- matrix(nrow=k,ncol=n,dimnames=list(grpnames,rownames(newdata)))
    for (g in 1:k)
    {
      sphdata <- newdata %*% object@scaling[,,g] 
      sphmeang <- object@means[g,] %*% object@scaling[,,g]
      Mahdistovertwo[g,] <- apply(sphdata, 1, function(x) sum((x-sphmeang)^2)/2)
    } 
#    wghtdensities <- sweep(exp(sweep(-Mahdistovertwo,1,STATS=object@ldet,FUN="-")),1,STATS=prior,FUN="*")
    minhlfMD2 <- apply(Mahdistovertwo,2,min)
    nrmhlfMD2 <- sweep(Mahdistovertwo,2,minhlfMD2)
#    wghtdensities <- sweep(exp(sweep(-nrmhlfMD2,1,object@ldet-minhlfMD2)),1,prior,"*")
    wghtdensities <- sweep(exp(sweep(-nrmhlfMD2,1,object@ldet)),1,prior,"*")
    ncnst <- apply(wghtdensities,2,sum)  			# normalizing constants
    posterior <- sweep(wghtdensities,2,STATS=ncnst,FUN="/")
    NAind <- which(apply(posterior,2,function(x) any(!is.finite(x))))
    if (length(NAind)>0) {
      NAMDover2 <- Mahdistovertwo[,NAind,drop=FALSE]
      minNAMDover2 <- apply(NAMDover2,2,min)
      normNAMDover2 <- scale(NAMDover2,center=minNAMDover2,scale=FALSE)
      NAwghtdensities <- sweep(exp(sweep(-normNAMDover2,1,STATS=object@ldet,FUN="-")),1,STATS=prior,FUN="*")
      NAncnst <- apply(NAwghtdensities,2,sum)  			# normalizing constants
      posterior[,NAind] <- sweep(NAwghtdensities,2,STATS=NAncnst,FUN="/")
    }
    clres <- apply(posterior, 2, function(pst) return(grpnames[which.max(pst)]))

    list(class=factor(clres,levels=grpnames),posterior=t(posterior))
  }
) 

setMethod("show",					
  signature(object = "Idtqda"),
  function(object)
  {
    cat("Prior probabilities of groups:\n") ; print(object@prior) ; cat("\n")
    cat("Group means:\n") ; print(object@means) ; cat("\n")
  }
)

setMethod("CovCase",					
  signature(object = "Idtlda"),
  function(object)
  {
    object@CovCase
  }
)

setMethod("CovCase",					
  signature(object = "Idtqda"),
  function(object)
  {
    object@CovCase
  }
)
