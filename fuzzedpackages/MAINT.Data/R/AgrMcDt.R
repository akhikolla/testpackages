AgrMcDt <- function(MicDtDF,agrby,agrcrt="minmax")
{
  mcall <- match.call()$MicDtDF
  if (length(mcall) > 1) mcall <- "microdata"
  if (!(is.data.frame(MicDtDF))) stop("First argument of AgMicroData must be a data frame.\n")
  if (class(MicDtDF)[1]!="data.frame") MicDtDF <- as.data.frame(MicDtDF)
  if (any(!sapply(1:ncol(MicDtDF),function(ind) is.numeric(MicDtDF[,ind])))) {  
    stop(paste("Some of the columns of the",mcall,"data frame have non-numeric variables.\n"))
  }
  
  unvalidobs <- which(apply(MicDtDF,1,function(v) any(!is.finite(v))))
  nunvalid <- length(unvalidobs) 
  if (nunvalid>0) {
    MicDtDF <- MicDtDF[-unvalidobs,]
    agrby <- agrby[-unvalidobs]
    string2 <- paste("rows of the",mcall,"data frame were dropped because they included non-valid (non finite or missing values) observations.\n")
    if (nunvalid<=10) warning(paste("The",paste(row.names(MicDtDF)[unvalidobs],collapse=" "),string2))
    else warning(paste(nunvalid,string2,collapse=" "))
  }

  if (!is.factor(agrby)) stop("Argument agrby is not a factor\n")
  globaln <- nrow(MicDtDF)
  if (length(agrby)!=globaln) stop("Size of the agrby argument does not agree with the number of rows in the MicDtDF data frame.\n") 
  if ( agrcrt[1]!="minmax" && (class(agrcrt)[1]!="numeric" || length(agrcrt)!=2 || agrcrt[1]>=agrcrt[2] || agrcrt[1]<0. || agrcrt[2]>1.) )
    stop(paste("Wrong value for the agrcrt argument\n( it should be either the string minmax or a two-dim vector",
               "\nof a prob. value for the lower percentile, followed by the prob. value for the upper percentile - \nex:c(0.05,0.95) ).\n")) 
   
  if (length(unique(agrby))!=length(levels(agrby)))  agrby <- factor(agrby)
  grplvls <- levels(agrby)
  lbDF <- ubDF <- data.frame(MicDtDF[1,])
  bndsDF <- cbind.data.frame(lbDF,ubDF)

  ngrps <- length(grplvls)
  nvar <- ncol(MicDtDF)
  NbMicroUnits <- integer(ngrps)
  for (r in 1:ngrps) 
  { 
    grp <- grplvls[r]
    rind <- which(agrby==grp)
    NbMicroUnits[r] <- length(rind)
    for (c in 1:nvar) {
      if (agrcrt[1]=="minmax") {
        bndsDF[r,c] <- min(MicDtDF[rind,c])
        bndsDF[r,nvar+c] <- max(MicDtDF[rind,c])
      } else {
        bndsDF[r,c] <- quantile(MicDtDF[rind,c],probs=agrcrt[1])
        bndsDF[r,nvar+c] <- quantile(MicDtDF[rind,c],probs=agrcrt[2])
      }
    }  
  }
  res <- IData(bndsDF,Seq="AllLb_AllUb",VarNames=names(MicDtDF),ObsNames=grplvls)
  DegInT <- which(apply(res@LogR,1,function(v) any(!is.finite(v))))
  nDegInT <- length(DegInT)
  if (nDegInT>0) {
    if (nDegInT==res@NObs) {
      warning("No Idata object was created because all units had some degenerate intervals")
      return(NULL)
    }
    if (nDegInT<10) {
      if (nDegInT==1) {
        wmsg <- paste("Data unit",res@ObsNames[DegInT],"was eliminated because it lead to some degenerate intervals")
      } else {
        wmsg <- paste(
          "Data units",paste(res@ObsNames[DegInT],collapse=", "),"were eliminated because they lead to some degenerate intervals",sep="\n"
        )
      }  
    } else {
      wmsg <- paste(nDegInT,"data units were eliminated because they lead to some degenerate intervals")
    }
    warning(wmsg)
    res <- res[-DegInT,]
  }
  res@NbMicroUnits <- NbMicroUnits[-DegInT] 
  if (length(DegInT)>0) {
    res@NbMicroUnits <- NbMicroUnits[-DegInT]
  } else {
    res@NbMicroUnits <- NbMicroUnits
  }  
  names(res@NbMicroUnits) <- res@ObsNames
  
  res
}

