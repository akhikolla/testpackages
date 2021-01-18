setMethod("Idtmclust",
  signature(Idt = "IData"),
  function(Idt, G=1:9, CovCase=1:4, SelCrit=c("BIC","AIC"), Mxt=c("Hom","Het","HomandHet"), control=EMControl())   
  {
    pertubzVct <- function(z,maxprt) {
      maxzind <- which.max(z)
      change <- runif(1,0.,min(z[maxzind],maxprt))
      z[maxzind] <- z[maxzind] - change
      z[-maxzind] <- z[-maxzind] + change/(length(z)-1)
      z
    } 
    pertubzMat <- function(OrigzMat,maxprt=0.1) t(apply(OrigzMat,1,pertubzVct,maxprt=maxprt))
    
    call <- match.call()
    SelCrit <- match.arg(SelCrit)
    Mxt <- match.arg(Mxt)
    k2max <- control@k2max
    maxiter <- control@maxiter 
    set.seed(control@seed)
    pertubfct <- control@pertubfct
    protol <- control@protol
    convtol <- control@convtol
    nrep <- control@nrep
    MaxVarGRt <- control@MaxVarGRt
                                 
    MaxSctEgvlRt <- Inf                                     # Not implemented yet!!
    if (!is.finite(MaxSctEgvlRt)) SctEgvCnstr <- FALSE 
    else SctEgvCnstr <- TRUE 

    if (class(G)!="integer") G <- as.integer(G)
    if (class(CovCase)!="integer") CovCase <- as.integer(CovCase)

    X <- data.frame(c(Idt@MidP,Idt@LogR),row.names=rownames(Idt))
    n <- nrow(X)
    p <- ncol(X)
    q <- p/2
    if (q==1) CovCase <- q1CovCase(CovCase) 
    nCases <- length(CovCase)
    nG <- length(G)
    ntot <- nG*nCases
    GCnames <- paste(rep(paste("G",G,sep=""),each=nCases),paste("C",CovCase,sep=""),sep="")
    Vnames <- names(X)
    Onames <- rownames(X)

    if (Mxt=="HomandHet" || Mxt=="Hom") {
      RepresHom <- vector("list",ntot)
      HomlogLiks <- rep(-Inf,ntot)
      HomBICs <- rep(Inf,ntot) 
      HomAICs <- rep(Inf,ntot)
      names(RepresHom) <- names(HomlogLiks) <- names(HomBICs) <- names(HomAICs) <- paste("Hom",GCnames,sep="")
    } else {
      RepresHom <- HomlogLiks <- HomBICs <- HomAICs <- NULL
    }
    if (Mxt=="HomandHet" || Mxt=="Het") {
      RepresHet <- vector("list",ntot)
      HetlogLiks <- rep(-Inf,ntot) 
      HetBICs <- rep(Inf,ntot)  
      HetAICs <- rep(Inf,ntot) 
      names(RepresHet) <- names(HetlogLiks) <- names(HetBICs) <- names(HetAICs) <- paste("Het",GCnames,sep="")
    } else {
      RepresHet <- HetlogLiks <- HetBICs <- HetAICs <- NULL
    }
    ind0 <- 0
    for (k in G) {

      if (k==1) {
        for (Cfi in nCases:1) {
          ind <- ind0+Cfi
          Cf <- CovCase[Cfi]
          mleEst <- mle(Idt,CovCase=Cf)
          muSig <- coef(mleEst)
          if ( CheckSigmaSing(Cf,muSig$Sigma,limlnk2=log(k2max),scale=TRUE) == TRUE ) {
            warning(paste("No one group solution was produced for configuration case",Cf,
                    "\nbecause the covariance estimate was found to be numerically singular",
                    "\nYou can try to force a solution by increasing the maximum allowed",
                     "\ncovariance condition number (EMControl argument k2max, currently set at",k2max,")",
                     "\nbut be awere that a larger k2max value may lead to spurious 'optimal' solutions.\n"))
            Lik <- -Inf
            BIC <- AIC <- Inf
          } else {     
            Lik <- mleEst@logLiks[mleEst@BestModel]
            BIC <- mleEst@BICs[mleEst@BestModel]
            AIC <- mleEst@AICs[mleEst@BestModel]
          }
          if (Cf==1) Sigmapar <- p*(p+1)/2 
          else if (Cf==2) Sigmapar <- 3*p/2	
          else if (Cf==3) Sigmapar <- p*(p/2+1)/2   
          else if (Cf==4) Sigmapar <- p
          if (Mxt=="HomandHet" || Mxt=="Hom") {
            clusters <- rep("CP1",Idt@NObs)
            names(clusters) <- Idt@ObsNames
            RepresHom[[ind]] <- new("IdtMclustEl",NObs=Idt@NObs,NIVar=Idt@NIVar,SelCrit=SelCrit,Hmcdt=TRUE,Conf=Cf,nG=k,
              logLik=Lik,alllnLik=Lik,bic=BIC,aic=AIC,z=matrix(1.,nrow=n,ncol=1,dimnames=list(Onames,"CP1")),classification=clusters,
              parameters=list(pro=NULL,mean=matrix(muSig$mu,ncol=1,dimnames=list(Vnames,"CP1")),
                covariance=array(muSig$Sigma,dim=c(p,p,1),dimnames=list(Vnames,Vnames,NULL)))
              )
            HomlogLiks[ind] <- Lik
            HomBICs[ind] <- BIC
            HomAICs[ind] <- AIC
          }
          if (Mxt=="HomandHet" || Mxt=="Het") {
            clusters <- rep("CP1",Idt@NObs)
            names(clusters) <- Idt@ObsNames
            RepresHet[[ind]] <- new("IdtMclustEl",NObs=Idt@NObs,NIVar=Idt@NIVar,SelCrit=SelCrit,Hmcdt=FALSE,Conf=Cf,nG=k,
              logLik=Lik,alllnLik=Lik,bic=BIC,aic=AIC,z=matrix(1.,nrow=n,ncol=1,dimnames=list(Onames,"CP1")),classification=clusters,
              parameters=list(pro=NULL,mean=matrix(muSig$mu,ncol=1,dimnames=list(Vnames,"CP1")),
                covariance=array(muSig$Sigma,dim=c(p,p,1),dimnames=list(Vnames,Vnames,"CP1")))
              )
            HetlogLiks[ind] <- Lik
            HetBICs[ind] <- BIC
            HetAICs[ind] <- AIC
          }
        }  

      } else { 
        ptaufct <- 0.5   
        pR <- 0.1
        ptau <- pertubfct*rep(ptaufct/k,k)
        
        z <- matrix(nrow=n,ncol=k)
        if (Mxt=="HomandHet" || Mxt=="Hom") {
           McHomCf <- 1
           Mcsol <- Mclust(X,k,modelNames="EEE",verbose=FALSE)
           if ( is.null(Mcsol) || CheckSigmaSing(1,Mcsol$parameters$variance$Sigma,limlnk2=log(k2max),scale=TRUE) == TRUE ) {
             Mcsol <- Mclust(X,k,modelNames="EEI",verbose=FALSE)  # If unrestriced covariance is singular, try a diagonal covariance estimate 
             if ( is.null(Mcsol) || CheckSigmaSing(4,Mcsol$parameters$variance$Sigma,limlnk2=log(k2max),scale=TRUE) == TRUE )  {
               Mcsol <- Mclust(X,k,modelNames="EII",verbose=FALSE)  # If still singular, try a spheric covariance estimate 
               if ( is.null(Mcsol) || CheckSigmaSing(4,Mcsol$parameters$variance$Sigma,limlnk2=log(k2max),scale=TRUE) == TRUE )  
                 stop("Idtmclust ran into numerical difficulties to find a valid solution.",
                      "\nMaybe your data set is too small for the specified models.",
                      " \nTry running Idtmclust with a smaller number of components\n")
             }        
             McHomCf <- 4
           }

           HomISol <- list(tau=Mcsol$parameters$pro,muk=Mcsol$parameters$mean,Sigma=Mcsol$parameters$variance$Sigma,Sigmak=NULL,
                           z=Mcsol$z,clusters=Mcsol$classification,LnLik=Mcsol$loglik,npar=NULL,BIC=NULL,AIC=NULL)

        }  
        if (Mxt=="HomandHet" || Mxt=="Het") {
           TrueHetSol <- TRUE
           McHetCf <- 1
           Mcsol <- Mclust(X,k,modelNames="VVV",verbose=FALSE)
           if ( is.null(Mcsol) || CheckSigmak(1,Mcsol$parameters$variance$sigma,MaxVarGRt=MaxVarGRt,limlnk2=log(k2max),scale=TRUE)[[1]] != "Valid" ) {
             Mcsol <- Mclust(X,k,modelNames="VVI",verbose=FALSE)  # If unrestriced covariances are singular, try a diagonal covariance estimates
             if ( is.null(Mcsol) || CheckSigmak(4,Mcsol$parameters$variance$sigma,MaxVarGRt=MaxVarGRt,limlnk2=log(k2max),scale=TRUE)[[1]] != "Valid" )  {
               Mcsol <- Mclust(X,k,modelNames="VII",verbose=FALSE)  # If still singular, try a spheric covariance estimate 
               if (is.null(Mcsol) || CheckSigmak(4,Mcsol$parameters$variance$sigma,MaxVarGRt=MaxVarGRt,limlnk2=log(k2max),scale=TRUE)[[1]] != "Valid" ) { 
               # If it still finds numerical difficulties, try to initialize with homoscedastic solutions
                 TrueHetSol <- FALSE
                 Mcsol <- Mclust(X,k,modelNames="EEE",verbose=FALSE)
                 if ( is.null(Mcsol) || CheckSigmaSing(1,Mcsol$parameters$variance$Sigma,limlnk2=log(k2max),scale=TRUE) == TRUE ) {
                   Mcsol <- Mclust(X,k,modelNames="EEI",verbose=FALSE)
                   if ( is.null(Mcsol) || CheckSigmaSing(1,Mcsol$parameters$variance$Sigma,limlnk2=log(k2max),scale=TRUE) == TRUE ) {
                     Mcsol <- Mclust(X,k,modelNames="EII",verbose=FALSE)
                   }
                   McHetCf <- 4
                 }  
                 Sigmak <- array(rep(Mcsol$parameters$variance$Sigma,k),dim=c(p,p,k))
               } else {
                 McHetCf <- 4
               }
             }
           }
           if (TrueHetSol) {                   # Retrieve Sigmak array for genuine heteroscedastic solutions
             Sigmak <- Mcsol$parameters$variance$sigma
           }
           HetISol <- list(tau=Mcsol$parameters$pro,muk=Mcsol$parameters$mean,Sigma=NULL,Sigmak=Sigmak,
                           z=Mcsol$z,clusters=Mcsol$classification,LnLik=Mcsol$loglik,npar=NULL,BIC=NULL,AIC=NULL)
        }  

        Vnames <- names(X)
        Onames <- rownames(X)
        Cnames <- paste("CP",1:k,sep="")

        for (Cfi in nCases:1) {
          Cf <- CovCase[Cfi]
          ind <- ind0+Cfi

          maxzpert <- 0.1
          maxinittrials <- 100
          inittrial <- 0
          if (Mxt=="HomandHet" || Mxt=="Hom") {

            if (Cf==McHomCf) LnLik <- HomISol$LnLik
            else {
              if (Cf==1) Sigma <- HomISol$Sigma
              else Sigma <- SetCovZeros(Cf,HomISol$Sigma) 
              LnLik <- MClusLikpar(X,n,p,k,Homoc=TRUE,tau=HomISol$tau,muk=HomISol$muk,Sigma=Sigma,k2max=k2max) 
            }

            ISol <- NULL 
            startwithM <- FALSE

            stdv <- sqrt(diag(HomISol$Sigma))
            sdmuk <- pertubfct*rep(stdv,k)
            sdStdev <- pertubfct*stdv/4
            if (Cf==1) sdR <- pertubfct*rep(pR,p*(p-1)/2)
            else if (Cf==2) sdR <- pertubfct*rep(pR,p)
            else if (Cf==3) sdR <- pertubfct*rep(pR,q*(q-1))
            else if (Cf==4) sdR <- NULL

            while ( is.null(ISol)  || !is.finite(ISol$LnLik) || any(!is.finite(ISol$tau)) || any(!is.finite(ISol$muk)) || any(!is.finite(ISol$Sigma)) )  
            {
              inittrial <- inittrial + 1

              ISol <- .Call( "CEMGauss", as.matrix(X), k, Cf, TRUE, maxiter, protol, convtol, k2max, MaxSctEgvlRt,
                 HomISol$z, HomISol$tau, HomISol$muk, HomISol$Sigma, NULL, LnLik, startwithM, SctEgvCnstr, MaxVarGRt,
                 PACKAGE = "MAINT.Data" )

              if ( k!=min(G) ) {
                prvk <- max(G[G<k])
                cmpGspos <- which(G==prvk)
                cmpGInd <- (cmpGspos-1)*nCases+Cfi
              } else cmpGInd <- NULL

              if ( (Cf==2 || Cf==3) && is.element(4,CovCase) ) cmpCfInd <- ind0 + nCases
              else if (Cf==1 && nCases > 1) cmpCfInd <- ind0 + which.max(HomlogLiks[ind0+CovCase])
                   else cmpCfInd <- NULL

              CfCrctd <- CGCrctd <- FALSE
              while (  ( !CfCrctd && (!is.null(cmpCfInd) && (is.null(ISol) || !is.finite(ISol$LnLik) || ISol$LnLik < HomlogLiks[cmpCfInd])) ) ||
                       ( !CGCrctd && (!is.null(cmpGInd) && (is.null(ISol) || !is.finite(ISol$LnLik) || ISol$LnLik < HomlogLiks[cmpGInd])) ) 
                  )
              {

                if ( !CfCrctd && (CGCrctd || (is.null(cmpGInd)  || (!is.null(cmpCfInd) && HomlogLiks[cmpCfInd] > HomlogLiks[cmpGInd]))) )
                {

                  ISol <- .Call( "CEMGauss", as.matrix(X), k, Cf, TRUE, maxiter, protol, convtol, k2max, MaxSctEgvlRt,
                    RepresHom[[cmpCfInd]]@z, RepresHom[[cmpCfInd]]@parameters$pro,
                    RepresHom[[cmpCfInd]]@parameters$mean, RepresHom[[cmpCfInd]]@parameters$covariance[,,1], NULL, HomlogLiks[cmpCfInd],
                    startwithM, SctEgvCnstr, MaxVarGRt,
                    PACKAGE = "MAINT.Data" )

                  CfCrctd <- TRUE 

                }  else if (!CGCrctd) {

                  newSol <- Addgrp(X,RepresHom[[cmpGInd]]@z,Cf,TRUE,k2max,protol,MaxVarGRt,k-prvk)
                  if (!is.null(newSol) && newSol$LnLik > ISol$LnLik) { 

                    ISol <- .Call( "CEMGauss", as.matrix(X), k, Cf, TRUE, maxiter, protol, convtol, k2max, MaxSctEgvlRt,
                      newSol$z, newSol$tau, newSol$muk, newSol$Sigma, NULL, newSol$LnLik, startwithM, SctEgvCnstr, MaxVarGRt,
                      PACKAGE = "MAINT.Data" )
                  }
                  CGCrctd <- TRUE 
                }                
              }

              if (is.null(ISol) || !is.finite(ISol$LnLik) || any(!is.finite(ISol$tau)) || any(!is.finite(ISol$muk)) || any(!is.finite(ISol$Sigma)) ) {
                if (inittrial >= maxinittrials) {
                  stop(paste("Idtmclust failed to find a valid Homoscedastic ",k,"-group solution for configuration C",CovCase," after ",inittrial," trials.\n",sep=""))
                } else { 
                  HomISol$z <- Pertubz(HomISol$z,pertubz=matrix(rep(ptau*1.02^inittrial,n),nrow=n,ncol=k))
                  startwithM <- TRUE
                }  
              }  
            } 

            rownames(ISol$muk) <- rownames(ISol$Sigma) <- colnames(ISol$Sigma) <- Vnames
            rownames(ISol$z) <- Onames
            names(ISol$tau) <- colnames(ISol$muk) <- colnames(ISol$z) <- Cnames
            ISol$clusters <- Cnames[ISol$clusters] 
            names(ISol$clusters) <- Onames

            if (nrep > 0) {

              FSol <- RepEMGauss(X,n,p,k,ISol,pertub=list(z=NULL,tau=ptau,muk=sdmuk,Stdev=sdStdev,cor=sdR),
                nrep=nrep,Cf=Cf,Homoc=TRUE,maxiter=maxiter,protol=protol,k2max=k2max,
                MaxSctEgvlRt=MaxSctEgvlRt,SctEgvCnstr=SctEgvCnstr,MaxVarGRt=MaxVarGRt,convtol=convtol)

              RepresHom[[ind]] <- new("IdtMclustEl",NObs=Idt@NObs,NIVar=Idt@NIVar,SelCrit=SelCrit,Hmcdt=TRUE,Conf=Cf,nG=k,
                logLik=FSol$BestSol$LnLik,alllnLik=FSol$alllnLik,bic=FSol$BestSol$BIC,aic=FSol$BestSol$AIC,
                parameters=list(
                  pro=FSol$BestSol$tau,mean=FSol$BestSol$muk,covariance=array(FSol$BestSol$Sigma,dim=c(p,p,1),dimnames=list(Vnames,Vnames,NULL))
                ),
                z=FSol$BestSol$z,classification=FSol$BestSol$clusters
              )
              HomlogLiks[ind] <- FSol$BestSol$LnLik
              HomBICs[ind] <- FSol$BestSol$BIC
              HomAICs[ind] <- FSol$BestSol$AIC

            } else {

              RepresHom[[ind]] <- new("IdtMclustEl",NObs=Idt@NObs,NIVar=Idt@NIVar,SelCrit=SelCrit,Hmcdt=TRUE,Conf=Cf,nG=k,
                logLik=ISol$LnLik,alllnLik=ISol$LnLik,bic=ISol$BIC,aic=ISol$AIC,
                parameters=list(
                  pro=ISol$tau,mean=ISol$muk,covariance=array(ISol$Sigma,dim=c(p,p,1),dimnames=list(Vnames,Vnames,NULL))
                ),
                z=ISol$z,classification=ISol$clusters
              )
              HomlogLiks[ind] <- ISol$LnLik
              HomBICs[ind] <- ISol$BIC
              HomAICs[ind] <- ISol$AIC
            }
              
          }

          if (Mxt=="HomandHet" || Mxt=="Het") {

              if (Cf==McHetCf) LnLik <- HetISol$LnLik
              else {
                Sigmak <- HetISol$Sigmak
                if (Cf!=1) for (g in 1:k) Sigmak[,,g] <- SetCovZeros(Cf,Sigmak[,,g]) 
                LnLik <- MClusLikpar(X,n,p,k,Homoc=FALSE,tau=HetISol$tau,muk=HetISol$muk,Sigmak=Sigmak,k2max=k2max)
              } 
           
            ISol <- NULL 
            startwithM <- FALSE

            sdmuk <- numeric(p*k)
            sdStdev <- numeric(p*k)
            for (g in 1:k) {
              stdv <- sqrt(diag(HetISol$Sigmak[,,g]))
              sdmuk[(g-1)*p+1:p] <- pertubfct*stdv
              sdStdev[(g-1)*p+1:p] <- pertubfct*stdv/4
            }
            if (Cf==1) sdR <- pertubfct*rep(pR,k*p*(p-1)/2)
            else if (Cf==2) sdR <- pertubfct*rep(pR,k*p)
            else if (Cf==3) sdR <- pertubfct*rep(pR,k*q*(q-1))
            else if (Cf==4) sdR <- NULL

            while ( is.null(ISol)  || !is.finite(ISol$LnLik) || any(!is.finite(ISol$tau)) || any(!is.finite(ISol$muk)) || any(!is.finite(ISol$Sigmak)) )  
            {
              inittrial <- inittrial + 1

              ISol <- .Call( "CEMGauss", as.matrix(X), k, Cf, FALSE, maxiter, protol, convtol, k2max, MaxSctEgvlRt,
                as.matrix(HetISol$z), HetISol$tau, as.matrix(HetISol$muk), NULL, HetISol$Sigmak, LnLik, startwithM, SctEgvCnstr, MaxVarGRt,
                PACKAGE = "MAINT.Data" )

              if ( k!=min(G) ) {
                prvk <- max(G[G<k])
                cmpGspos <- which(G==prvk)
                cmpGInd <- (cmpGspos-1)*nCases+Cfi
              } else cmpGInd <- NULL
              
              if ( (Cf==2 || Cf==3) && is.element(4,CovCase) ) cmpCfInd <- ind0 + nCases
              else if (Cf==1 && nCases > 1) cmpCfInd <- ind0 + which.max(HetlogLiks[ind0+CovCase])
                   else cmpCfInd <- NULL

              HomCrctd <- CfCrctd <- CGCrctd <- FALSE
              while (  ( !HomCrctd && (!is.null(HomlogLiks) && (is.null(ISol) || ISol$LnLik < HomlogLiks[ind])) )  ||
                       ( !CfCrctd && (!is.null(cmpCfInd) && (is.null(ISol) || !is.finite(ISol$LnLik) || ISol$LnLik < HetlogLiks[cmpCfInd])) ) ||
                       ( !CGCrctd && (!is.null(cmpGInd) && (is.null(ISol) || !is.finite(ISol$LnLik) || ISol$LnLik < HetlogLiks[cmpGInd])) ) 
                  )
              {
              if ( !HomCrctd && 
                   ( CfCrctd || is.null(cmpCfInd) || (!is.null(HomlogLiks) && HomlogLiks[ind] > HetlogLiks[cmpCfInd]) ) &&
                   ( CGCrctd || is.null(cmpGInd)  || (!is.null(HomlogLiks) && HomlogLiks[ind] > HetlogLiks[cmpGInd]) ) 
                 ) {
                  
                  Sigmak <- array(rep(RepresHom[[ind]]@parameters$covariance,k),dim=c(p,p,k))

                  ISol <- .Call( "CEMGauss", as.matrix(X), k, Cf, FALSE, maxiter, protol, convtol, k2max, MaxSctEgvlRt,
                    RepresHom[[ind]]@z, RepresHom[[ind]]@parameters$pro,
                    RepresHom[[ind]]@parameters$mean, NULL, Sigmak, HomlogLiks[ind], FALSE, SctEgvCnstr, MaxVarGRt,
                    PACKAGE = "MAINT.Data" )

                  HomCrctd <- TRUE

                }  else if ( !CfCrctd && (CGCrctd || (is.null(cmpGInd)  || (!is.null(cmpCfInd) && HetlogLiks[cmpCfInd] > HetlogLiks[cmpGInd]))) ) {
                
                  ISol <- .Call( "CEMGauss", as.matrix(X), k, Cf, FALSE, maxiter, protol, convtol,
                    k2max, MaxSctEgvlRt,
                    RepresHet[[cmpCfInd]]@z, RepresHet[[cmpCfInd]]@parameters$pro,
                    RepresHet[[cmpCfInd]]@parameters$mean, NULL, RepresHet[[cmpCfInd]]@parameters$covariance, HetlogLiks[cmpCfInd], 
                    startwithM, SctEgvCnstr, MaxVarGRt, 
                    PACKAGE = "MAINT.Data" )

                  CfCrctd <- TRUE 

                }  else if (!CGCrctd) {
                  
                    newSol <- Addgrp(X,RepresHet[[cmpGInd]]@z,Cf,FALSE,k2max,protol,MaxVarGRt,k-prvk)
                    if (!is.null(newSol) && newSol$LnLik > ISol$LnLik) { 
                      ISol <- .Call( "CEMGauss", as.matrix(X), k, Cf, FALSE, maxiter, protol, convtol, k2max, MaxSctEgvlRt,
                        newSol$z, newSol$tau, newSol$muk, NULL, newSol$Sigmak, newSol$LnLik, startwithM, SctEgvCnstr, MaxVarGRt,
                        PACKAGE = "MAINT.Data" )
                    }
                   CGCrctd <- TRUE 
                }                
              }

              if (is.null(ISol) || !is.finite(ISol$LnLik) || any(!is.finite(ISol$tau)) || any(!is.finite(ISol$muk)) || any(!is.finite(ISol$Sigmak)) ) {
                if (inittrial >= maxinittrials) {
                  stop(paste("Idtmclust failed to find a valid Heteroscedastic ",k,"-group solution for configuration C",CovCase," after ",inittrial," trials.\n",sep=""))
                } else {
                  HetISol$z <- Pertubz(HetISol$z,pertubz=matrix(rep(ptau*1.02^inittrial,n),nrow=n,ncol=k))
                  startwithM <- TRUE
                }
              }  
            }

            rownames(ISol$muk) <- dimnames(ISol$Sigmak)[[1]]  <- dimnames(ISol$Sigmak)[[2]] <- Vnames
            dimnames(ISol$Sigmak)[[3]] <- Cnames
            rownames(ISol$z) <- Onames
            names(ISol$tau) <- colnames(ISol$muk) <- colnames(ISol$z) <- Cnames
            ISol$clusters <- Cnames[ISol$clusters] 
            names(ISol$clusters) <- Onames

            if (nrep>0) {

              FSol <- RepEMGauss(X,n,p,k,ISol,pertub=list(z=NULL,tau=ptau,muk=sdmuk,Stdev=sdStdev,cor=sdR),
                nrep=nrep,Cf=Cf,Homoc=FALSE,maxiter=maxiter,protol=protol,k2max=k2max,
                MaxSctEgvlRt=MaxSctEgvlRt,SctEgvCnstr=SctEgvCnstr,MaxVarGRt=MaxVarGRt,convtol=convtol)
            
              RepresHet[[ind]] <- new("IdtMclustEl",NObs=Idt@NObs,NIVar=Idt@NIVar,SelCrit=SelCrit,Hmcdt=FALSE,Conf=Cf,nG=k,
                logLik=FSol$BestSol$LnLik,alllnLik=FSol$alllnLik,bic=FSol$BestSol$BIC,aic=FSol$BestSol$AIC,
                parameters=list(pro=FSol$BestSol$tau,mean=FSol$BestSol$muk,covariance=FSol$BestSol$Sigmak),
                z=FSol$BestSol$z,classification=FSol$BestSol$clusters
              )
              HetlogLiks[ind] <- FSol$BestSol$LnLik
              HetBICs[ind] <- FSol$BestSol$BIC
              HetAICs[ind] <- FSol$BestSol$AIC
             
            } else {  

              RepresHet[[ind]] <- new("IdtMclustEl",NObs=Idt@NObs,NIVar=Idt@NIVar,SelCrit=SelCrit,Hmcdt=FALSE,Conf=Cf,nG=k,
                logLik=ISol$LnLik,alllnLik=ISol$LnLik,bic=ISol$BIC,aic=ISol$AIC,
                parameters=list(pro=ISol$tau,mean=ISol$muk,covariance=ISol$Sigmak),
                z=ISol$z,classification=ISol$clusters
              )
              HetlogLiks[ind] <- ISol$LnLik
              HetBICs[ind] <- ISol$BIC
              HetAICs[ind] <- ISol$AIC
             
            }
          }
        }
      }
      ind0 <- ind0 + nCases
    }  
    if (Mxt=="HomandHet" || Mxt=="Hom") {
      if (SelCrit=="BIC") {
        bestHomMod <- which.min(HomBICs)
        bestHomCrit <- HomBICs[bestHomMod]
      } else if (SelCrit=="AIC") {
        bestHomMod <- which.min(HomAICs)
        bestHomCrit <- HomAICs[bestHomMod]
      }
    } else {
      bestHomMod <- bestHomCrit <- NULL
    }
    if (Mxt=="HomandHet" || Mxt=="Het") {
      if (SelCrit=="BIC") {
        bestHetMod <- which.min(HetBICs)
        bestHetCrit <- HetBICs[bestHetMod]
      } else if (SelCrit=="AIC") {
        bestHetMod <- which.min(HetAICs)
        bestHetCrit <- HetAICs[bestHomMod]
      }
    } else {
      bestHetMod <- bestHetCrit <- NULL
    }
    if ( is.null(bestHetMod) || (!is.null(bestHomMod) && bestHomCrit < bestHetCrit) )
    {
      BestMxt <- "Hom"
    } else {
      BestMxt <- "Het"
    }
     if (BestMxt=="Hom") {
       return (
         new("IdtMclust",call=call,data=Idt,NObs=Idt@NObs,NIVar=Idt@NIVar,SelCrit=SelCrit,Hmcdt=TRUE,
           BestC=RepresHom[[bestHomMod]]@Conf,BestG=RepresHom[[bestHomMod]]@nG,
           logLiks=c(HomlogLiks,HetlogLiks),BICs=c(HomBICs,HetBICs),AICs=c(HomAICs,HetAICs),
           logLik=RepresHom[[bestHomMod]]@logLik,bic=RepresHom[[bestHomMod]]@bic,aic=RepresHom[[bestHomMod]]@aic,
           parameters=RepresHom[[bestHomMod]]@parameters,z=RepresHom[[bestHomMod]]@z,classification=RepresHom[[bestHomMod]]@classification,
           allres=list(RepresHom=RepresHom,RepresHet=RepresHet) )
         )
     }  else if (BestMxt=="Het") {  
       return (
         new("IdtMclust",call=call,data=Idt,NObs=Idt@NObs,NIVar=Idt@NIVar,SelCrit=SelCrit,Hmcdt=FALSE,
           BestC=RepresHet[[bestHetMod]]@Conf,BestG=RepresHet[[bestHetMod]]@nG,
           logLiks=c(HomlogLiks,HetlogLiks),BICs=c(HomBICs,HetBICs),AICs=c(HomAICs,HetAICs),
           logLik=RepresHet[[bestHetMod]]@logLik,bic=RepresHet[[bestHetMod]]@bic,aic=RepresHet[[bestHetMod]]@aic,
           parameters=RepresHet[[bestHetMod]]@parameters,z=RepresHet[[bestHetMod]]@z,classification=RepresHet[[bestHetMod]]@classification,
           allres=list(RepresHet=RepresHet,RepresHet=RepresHet) )
         )
     }
  }  
)

EMControl <- function (nrep=0, maxiter=1000, convtol=0.01, protol=1e-3, seed=NULL, pertubfct=1, k2max=1e6, MaxVarGRt=1e6)  
{
    new("EMControl", nrep=nrep, maxiter=maxiter, convtol=convtol, protol=protol, seed=seed, pertubfct=pertubfct,
                     k2max=k2max, MaxVarGRt=MaxVarGRt)
}

SetCovZeros <- function(Cf,Sigma)
{
  p <- ncol(Sigma)
  q <- p/2
  if (Cf==2) { for (r in 2:p) for (c in 1:(r-1)) if (r!=c+q) Sigma[r,c] <- Sigma[c,r] <- 0. }
  else if (Cf==3) Sigma[1:q,(q+1):p] <- Sigma[(q+1):p,1:q] <- 0.
       else if (Cf==4) Sigma[row(Sigma)!=col(Sigma)] <- 0.
            else stop("Wrong value for Cf argument\n")  

  Sigma
}

