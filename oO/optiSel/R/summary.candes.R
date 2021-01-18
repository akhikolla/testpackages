
"summary.candes"<-function(object, tlim=range(object$phen$Born, na.rm=TRUE), histNe=NA, base=tlim[1], df=4, ...){
  ### Check arguments ###
  if(class(object)!="candes"){stop("Argument object must be created with function candes.\n")}
  if(!("Born" %in% colnames(object$phen))){     stop("Column 'Born' with Year-of-Birth is missing.\n")}
  if(!("I" %in% colnames(object$phen))){        stop("Column 'I' with the generation interval is missing.\n")}
  if(!("Offspring" %in% colnames(object$phen))){stop("Column 'Offspring' indicating the individuals with offspring.\n")}
  if(!is.numeric(object$phen$Born)){            stop("Column 'Born' is not numeric.\n")}  
  if(!is.numeric(object$phen$I)){               stop("Column 'I' is not numeric.\n")}  
  if(!is.logical(object$phen$Offspring)){       stop("Column 'Offspring' is not logical\n")}  
  if(all(is.na(object$phen$Born))){             stop("Column 'Born' contains only NA.\n")}  
  if(all(is.na(object$phen$I))){                stop("Column 'I' contains only NA.\n")}  
  if(any(is.na(object$phen$Offspring))){        stop("Column 'Offspring' contains NA values.\n")}  
  if(length(tlim)!=2){    stop("Argument tlim must be a numerc vector with length 2.\n")}
  if(!is.numeric(tlim)){  stop("Argument tlim must be a numerc vector with length 2.\n")}
  if(any(is.na(tlim))){   stop("Argument tlim contains NA.\n")}
  if(tlim[1]>tlim[2]){    stop("Ensure that tlim[1]<=tlim[2].\n")}
  if(tlim[1]<min(object$phen$Born, na.rm=TRUE)){stop(paste("Argument tlim[1]<min(Born)=", min(object$phen$Born, na.rm=TRUE), ".\n"))}
  if(tlim[2]>max(object$phen$Born, na.rm=TRUE)){
    tlim[2] <- max(object$phen$Born, na.rm=TRUE)
    cat(paste("Argument tlim[2] is set to max(Born)=", tlim[2], ".\n"))
    }
  if(length(base)!=1 | !is.numeric(base) | is.na(base)){stop("Argument base must be a numeric value.\n")}
  if(base>tlim[1]){stop(paste("Argument base>tlim[1]=", tlim[1], ".\n"))}
  storage.mode(histNe)<-"double"
  if(length(histNe)!=1 || !is.numeric(histNe)){stop("Argument histNe must be a numeric value.\n")}
  if(is.na(histNe) && (base < tlim[1])){stop("If base < tlim[1], then histNe cannot be NA.\n")}
  if(!is.na(histNe) && histNe<=0){stop("Ensure that histNe>0.\n")}
  if(length(df)!=1 || !is.numeric(df) || is.na(df) || df<=0){stop("Argument df must be a positive numeric value.\n")}
  
  ### Estimate generation interval I and sample sizes #####################
  ### Set I for individuals without offspring to onknown
  #object$phen$I[!object$phen$Offspring] <- NA
  ### Estimate I for founders with linear regression ######################
  Koef <- lm(I~Born,data=object$phen)$coef
  isFounder <- is.na(object$phen$Sire) & is.na(object$phen$Dam)
  object$phen$I[isFounder]<- Koef[1]+Koef[2]*object$phen$Born[isFounder]
  ### Estimate spline function ###
  isComplete <- !is.na(object$phen$I)&!is.na(object$phen$Born)
  x        <- object$phen$Born[isComplete]
  y        <- object$phen$I[isComplete]
  #w        <- 0.05+0.95*object$phen$Offspring[isComplete]
  w        <- rep(1, length(x))
  mySpline <- sm.spline(x, y, w, df=df)
  date     <- sort(unique(object$phen$Born))
  I        <- predict(mySpline, date)
  Param    <- data.frame(t=date, I=I, row.names=date)
  Param    <- Param[Param$t>=tlim[1] & Param$t<=tlim[2],]
  N        <- nrow(object$phen)
  
  object$phen$Born[is.na(object$phen$Born)]<- -1234567
  
  nameNatK <- character(0)
  for(i in seq_along(object$kinship)){
    if(class(object$kinship[[i]])%in% c("quadFun","ratioFun")){
      name <- object$kinship[[i]]$name
      Param[[name]] <- NA
      if(class(object$kinship[[i]])=="ratioFun"){
        nameNatK <- c(nameNatK, name)
        }
    }
  }
  
  for(i in seq_along(Param$t)){
    GI   <- round(Param$I[i])
    t    <- Param$t[i]
    Use  <- (object$phen$Born <= t) & (object$phen$Born > t-GI)
    Nobs <- sum(Use)*sum(Use)-sum(Use)
    Use  <- matrix(Use, N, N, byrow=TRUE) & matrix(Use, N, N, byrow=FALSE)
    diag(Use) <- FALSE
    
    for(k in seq_along(object$kinship)){
      if(class(object$kinship[[k]])=="quadFun"){
        name <- object$kinship[[k]]$name
        Param[i, name] <- sum((object$kinship[[k]]$Q)[Use])/Nobs
      }
      if(class(object$kinship[[k]])=="ratioFun"){
        name <- object$kinship[[k]]$name
        Param[i, name] <- sum((object$kinship[[k]]$Q1)[Use])/sum((object$kinship[[k]]$Q2)[Use])
      }
    }
  }
  
  lambda <- ifelse(base==tlim[1], 1, (1-1/(2*histNe))^((tlim[1]-base)/(Param$I[1])))
  
  for(i in nameNatK){
    nameNe  <- "Ne"
    nameNGE <- "NGE"
    if(length(nameNatK)>1){
      nameNe <- paste0(nameNe,   "for", i)
      nameNGE <- paste0(nameNGE, "for", i)
    }
    Param[[nameNe]] <- getNe(Param$t, 1-Param[[i]], df=df, I=Param$I)
    condGD <- 1-Param[[i]]
    condGD  <- lambda * condGD/condGD[1]
    Param[[nameNGE]] <- 1/(2*(1 - condGD))
  }
  cat("\n")
  print(format(round(Param,3), digits=3, scientific=FALSE))
  attributes(Param)$histNe <- histNe
  attributes(Param)$base   <- base
  attributes(Param)$tlim   <- tlim
  attributes(Param)$df     <- df
  
  invisible(Param)
}