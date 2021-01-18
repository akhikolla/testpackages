DACrossVal <- function(data,grouping,TrainAlg,EvalAlg=EvalClrule,Strfolds=TRUE,kfold=10,CVrep=20,
  prior="proportions",loo=FALSE,...)
{
  fold <- function(n,kfold,fi) {
    lb <- floor((fi-1)*n/kfold) + 1
    ub <- floor(fi*n/kfold)
    if (lb>ub) return(NULL)
    lb:ub  # return(lb:ub)
  }   

  if (!is.factor(grouping)) {
    grouping <- factor(grouping)
  }
  codes <- levels(grouping)
  nk <- table(grouping)
  k <- nrow(nk)
  if (k<2) stop("Factor grouping should have at least two different levels\n")
  nk <- as.numeric(nk)
  n <- sum(nk)
  if (loo) {
    kfold <- n; CVrep <- 1
  } else if (Strfolds) {
    permut <- vector("list",k)
  }
  trep <- kfold*CVrep
  EvalRes <- array(dim=c(trep,k,2),dimnames=list(NULL,codes,c("Nk","Clerr")))
  if (prior[1]=="proportions") { prior <- nk/n }
  for (i in 1:CVrep)
  {
    if (!loo) {
      if (Strfolds) {
        for (grp in 1:k) permut[[grp]] <- sort.int(runif(nk[grp]),index.return=TRUE)$ix
      } else {
       permut <- sort.int(runif(n),index.return=TRUE)$ix
      }
    }
    for (j in 1:kfold)
    {
      if (loo) {
        out <- j
      } else {
        if (Strfolds) {
          out <- which(grouping==codes[1])[permut[[1]]][fold(nk[1],kfold,j)]
          for (grp in 2:k) out <- c(out,which(grouping==codes[grp])[permut[[grp]]][fold(nk[grp],kfold,j)])
        } else {
          out <- permut[fold(n,kfold,j)]
        }
      }
      if (length(out)==0) next
      tres <- try(TrainAlg(data[-out,,drop=FALSE],grouping[-out],prior=prior,...))
      if (is.null(tres) || class(tres)[1] == "try-error")
      {
        EvalResij <- list(err=NA,Nk=as.numeric(table(grouping[-out])))
        warning("Non valid classification rule in fold ",j," of replication ",i,"\n")
      }  else {
        EvalResij <- EvalAlg(tres,data[out,],grouping[out],k=k)
        rep <- (i-1)*kfold+j
        EvalRes[rep,,"Clerr"] <- ifelse(EvalResij$Nk>0.,EvalResij$err,NA)
        EvalRes[rep,,"Nk"] <- EvalResij$Nk
      }
    }
  }

  glberr <- sum(EvalRes[,,"Nk"]*EvalRes[,,"Clerr"],na.rm=TRUE)/sum(EvalRes[,,"Nk"],na.rm=TRUE)
  errestimates <- c(colMeans(EvalRes[,,"Clerr"],na.rm=TRUE),Global=glberr)
  
  mcall <- match.call(expand.dots=FALSE)
  msg <- paste("Error rate estimates of algorithm",mcall[[4]])
  ndotargs <- length(mcall$...)
  if (ndotargs>0) {
    msg <- paste(msg,"with argument")
    if(ndotargs>1) msg <- paste(msg,"s",sep="")
    for (i in 1:ndotargs) {
      msg <- paste(msg," ",names(mcall$...)[i],"=",deparse(mcall$...[[i]]),sep="")
    }
  }

  cat(msg,"\n")
  print(errestimates)
  attr(EvalRes,"errestimates") <- errestimates
  invisible(EvalRes)  # returns EvalRes silently
}

EvalClrule <- function(darule,VData,Vgrp,k)
{
  grpcodes <- levels(Vgrp)
  errates <- array(dim=k)
  Nk <- array(dim=k)
  clres <- predict(darule,VData)$class
  for (grpInd in 1:k)
  {
    Nk[grpInd] <- length(Vgrp[Vgrp==grpcodes[grpInd]])
    thisgrpclres <- clres[Vgrp==grpcodes[grpInd]]
    levels(thisgrpclres) <- grpcodes
    thisgrperr <- thisgrpclres[grpcodes[grpInd]!=thisgrpclres]
    if (Nk[grpInd]>0) {
      errates[grpInd] <- length(thisgrperr)/Nk[grpInd]
    } else {
      errates[grpInd] <- 0
    }
  }
  list(err=errates,Nk=Nk)  #  return(list(err=errates,Nk=Nk))
}

