
MacroPCApredict <- function(Xnew,InitialMacroPCA,MacroPCApars=NULL){
  #
  # This function performs the MacroPCA prediction on Xnew.
  #
  # The inputs are:
  #
  # Xnew            : the new data (test data), which must be a 
  #                   matrix or a data frame. 
  #                   It must always be provided.
  # InitialMacroPCA : the output of the MacroPCA function on the 
  #                   initial (training) dataset. Must be provided.
  # MacroPCApars    : the input options to be used for the prediction.
  #                   By default the options of InitialMacroPCA
  #                   are used. For the complete list of options
  #                   see the function MacroPCA.
  #
  # The outputs are:
  #
  # MacroPCApars: the options used in the call.
  # scaleX      : scales of the columns from InitialMacroPCA.
  # k           : number of PC's, from InitialMacroPCA.  
  # loadings    : the loadings from InitialMacroPCA.
  # eigenvalues : the eigenvalues from InitialMacroPCA.
  # center      : center from InitialMacroPCA.
  # n.obs       : number of cases of Xnew.
  # It          : number of iteration steps.
  # diff        : convergence criterion.
  # X.NAimp     : Xnew with all NA's imputed by MacroPCA.
  # scores      : scores of X.NAimp 
  # OD          : orthogonal distances of the rows of X.NAimp 
  # cutoffOD    : cutoff value for the OD, from InitialMacroPCA.
  # SD          : score distances of the rows of X.NAimp 
  # cutoffSD    : cutoff value for the SD, from InitialMacroPCA.
  # indrows     : row numbers of rowwise outliers in Xnew.
  # residScale  : scale of the residuals, from InitialMacroPCA.
  # stdResid    : standardized residuals of Xnew. Note that these
  #               are NA for all missing values of Xnew.
  # indcells    : indices of cellwise outliers of Xnew.
  # NAimp       : various results for the NA-imputed Xnew.
  # Cellimp     : various results for the cell-imputed Xnew.
  # Fullimp     : various result for the fully imputed Xnew.
  # DDC         : result of DDCpredict which is the first step of
  #               MacroPCApredict. See the function DDCpredict.
  
  
  if (is.null(MacroPCApars)) { 
    MacroPCApars = InitialMacroPCA$MacroPCApars
  } else {
    if (!is.list(MacroPCApars)) {
      stop("MacroPCApars must be a list")
    }
    InitialMacroPCA$MacroPCApars[names(MacroPCApars)] <- MacroPCApars
    MacroPCApars <- InitialMacroPCA$MacroPCApars
  }
  
  if (!"distprob" %in% names(MacroPCApars)) {
    MacroPCApars$distprob <- 0.99
  }
  if (!"maxiter" %in% names(MacroPCApars)) {
    MacroPCApars$maxiter <- 20
  }
  if (!"tol" %in% names(MacroPCApars)) {
    MacroPCApars$tol <- 0.005
  }
  if (!"bigOutput" %in% names(MacroPCApars)) {
    MacroPCApars$bigOutput <- TRUE
  }
  
  Xnew <- as.matrix(Xnew)
  n    <- nrow(Xnew)
  d    <- ncol(Xnew)
  
  # Input arguments
  InitialDDC <- InitialMacroPCA$DDC
  DDCpars    <- MacroPCApars$DDCpars
  scaleX     <- InitialMacroPCA$scaleX
  distprob   <- MacroPCApars$distprob
  maxiter    <- MacroPCApars$maxiter
  tol        <- MacroPCApars$tol
  bigOutput <- MacroPCApars$bigOutput
  
  # Step1: DDCpredict
  resultDDC <- DDCpredict(Xnew,InitialDDC,DDCpars)
  DDCpars   <- resultDDC$DDCpars
  DDCimp    <- sweep(resultDDC$Ximp,2,scaleX,"/")
  Xnew      <- sweep(Xnew,2,scaleX,"/")
  XO        <- Xnew   # matrix with NA's
  rm(Xnew)            # to save space
  Xfi       <- XO  # initialize fully imputed data matrix  
  
  indcells <- resultDDC$indcells
  indNA    <- which(is.na(XO)) # positions of NA's
  indimp   <- unique(c(indcells,indNA)) 
  # positions of cells that may be imputed
  
  Xfi[indimp] <- DDCimp[indimp]   
  # XOimp <- Xfi  # store imputed rowoutliers
  Xind <- Xfi; Xind[indimp] <- NA; Xind <- is.na(Xind)
  
  # retrieve parameters
  k        <- InitialMacroPCA$k
  loadings <- InitialMacroPCA$loadings # existing loadings matrix
  center   <- InitialMacroPCA$center/scaleX # scaled center
  eigenvalues <- InitialMacroPCA$eigenvalues
  
  # Iterations for new imputations/estimations
  if (any(Xind) & maxiter>0){ 
    # Step 3: Iterate
    diff <- 2*tol # Inf # 100
    It <- 0
    while (It < maxiter & diff > tol) { # Iterate
      It <- It+1;
      
      Xfimis <- Xfi[Xind] # missing values matrix
      XfiC <- sweep(Xfi,2,center) # centered Xfi
      
      Tr <- XfiC%*%loadings # scores matrix
      Xfihat <- Tr%*%t(loadings) # fit to XfiC
      Xfihat <- sweep(Xfihat,2,center,"+") # fit to Xfi
      Xfi[Xind] <- Xfihat[Xind] # impute missings
      
      diff <- mean((Xfi[Xind]-Xfimis)^2)
    } # Iterate end
  } else {
    # don't iterate when there are no values to impute
    diff <- 0
    It <- 0
  }
  
  dimnames(loadings) <- list(colnames(XO),paste("PC", 
                                                seq_len(ncol(loadings)),sep=""))
  
  # Start the output list: res
  res <- list(MacroPCApars=MacroPCApars,scaleX=scaleX,k=k,
              loadings=loadings,eigenvalues=eigenvalues,
              center=center,It=It,diff=diff) 
  
  # NA-imputed data Xnai : scores and distances
  Xnai        <- XO
  Xnai[indNA] <- Xfi[indNA] # impute missings
  scoresnai   <- (Xnai - matrix(rep(center,times=n),nrow=n, 
                                byrow=TRUE)) %*% loadings
  res$scores  <- scoresnai
  NAimp       <- list(scoresnai=scoresnai)
  # res$h        <- n # needed for next line:
  cutoffOD    <- InitialMacroPCA$cutoffOD
  out <- pca.distancesNew(res, Xnai, scoresnai, ncol(Xnai), distprob,
                          cutOD = cutoffOD)
  res$OD       <- out$OD
  out$cutoffOD <- cutoffOD
  res$cutoffOD <- out$cutoffOD
  res$SD       <- out$SD
  out$cutoffSD <- InitialMacroPCA$cutoffSD
  res$cutoffSD <- out$cutoffSD
  out$indrowsnai <- which(out$OD > out$cutoffOD)
  res$indrows    <- out$indrowsnai
  NAimp <- c(NAimp, out); rm(out)  
  
  
  if (bigOutput) {
    # Fully imputed data Xfi : scores and distances
    scoresfi <- (Xfi - matrix(rep(center,times=n),nrow=n, 
                              byrow=TRUE)) %*% loadings
    Fullimp <- list(scoresfi=scoresfi)
    out <- pca.distancesNew(res,Xfi,scoresfi,ncol(Xfi),distprob,
                            cutOD=cutoffOD)
    out$cutoffOD  <- cutoffOD
    out$cutoffSD  <- InitialMacroPCA$cutoffSD  
    out$indrowsfi  <- which(out$OD > out$cutoffOD)
    Fullimp <- c(Fullimp,out); rm(out)
    
    
    # Cell-imputed matrix Xci : scores and distances
    Xci <- Xfi
    Xci[res$indrows] <- Xnai[res$indrows]
    scoresci <- (Xci - matrix(rep(center,times=n),nrow=n, 
                              byrow = TRUE)) %*% loadings
    Cellimp <- list(scoresci=scoresci)
    out <- pca.distancesNew(res,Xci,scoresci,ncol(Xci),distprob,
                            cutOD=cutoffOD)
    out$cutoffOD  <- cutoffOD
    out$cutoffSD  <- InitialMacroPCA$cutoffSD    
    out$indrowsci <- which(out$OD > out$cutoffOD)   
    Cellimp <- c(Cellimp,out); rm(out)
  }
  
  ## Unstandardization and residuals
  
  # Compute residuals
  # C is for Centering:  
  XOC   <- sweep(XO,2,center)    # original matrix X, with NA's
  if (bigOutput) {
    XnaiC <- sweep(Xnai,2,center)
    XciC  <- sweep(Xci,2,center)
    XfiC  <- sweep(Xfi,2,center) 
  }
  # Compute standardized residuals of XO (with NA's): 
  stdResid <- XOC - (scoresnai%*%t(loadings))
  res$residScale <- InitialMacroPCA$residScale
  res$stdResid <- sweep(stdResid,2,res$residScale,"/")
  res$indcells <- which(abs(res$stdResid)>
                          sqrt(qchisq(DDCpars$tolProb,1)))
  
  if (bigOutput) {
    # Compute standardized residuals of NA-imputed data:
    stdResidnai <- XnaiC - (scoresnai%*%t(loadings))
    NAimp$residScalenai <- InitialMacroPCA$residScale
    NAimp$stdResidnai <- sweep(stdResidnai,2,NAimp$residScalenai,"/")
    NAimp$indcellsnai <- which(abs(NAimp$stdResidnai)>
                                 sqrt(qchisq(DDCpars$tolProb,1)))
    
    # Compute standardized residuals of cell-imputed data:
    stdResidci           <- XciC - (scoresci%*%t(loadings))
    Cellimp$residScaleci <- InitialMacroPCA$Cellimp$residScaleci
    Cellimp$stdResidci   <- sweep(stdResidci,2,Cellimp$residScaleci,"/")
    Cellimp$indcellsci   <- which(abs(Cellimp$stdResidci)>
                                    sqrt(qchisq(DDCpars$tolProb,1)))
    
    # Compute standardized residuals of fully imputed data:
    stdResidfi           <- XfiC - (scoresfi%*%t(loadings))
    Fullimp$residScalefi <- InitialMacroPCA$Fullimp$residScalefi
    Fullimp$stdResidfi   <- sweep(stdResidfi,2,Fullimp$residScalefi,"/")
    Fullimp$indcellsfi   <- which(abs(Fullimp$stdResidfi)>
                                    sqrt(qchisq(DDCpars$tolProb,1)))
  }
  ## unstandardize:
  res$X.NAimp <- sweep(Xnai,2,scaleX,"*")
  if (bigOutput) {
    Cellimp$Xci <- sweep(Xci,2,scaleX,"*")  
    Fullimp$Xfi <- sweep(Xfi,2,scaleX,"*")
  }
  res$center  <- center * scaleX
  
  # add remainder of output:
  # res$h        <- NULL
  # res$scaleX   <- scaleX
  if (bigOutput) {
    res$NAimp    <- NAimp
    res$Cellimp  <- Cellimp   
    res$Fullimp  <- Fullimp
  }
  res$DDC      <- resultDDC 
  
  names(res$scaleX) <- colnames(Xnai)
  names(res$center) <- colnames(Xnai)
  names(res$residScale) <- colnames(Xnai)
  
  if (bigOutput) {
    names(Cellimp$residScaleci) <- colnames(Xnai)
    names(Fullimp$residScalefi) <- colnames(Xnai)
    names(NAimp$residScalenai) <- colnames(Xnai)
  }
  
  return(res)
} # ends MacroPCApredict
