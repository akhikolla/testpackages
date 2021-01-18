`print.mdr` <- function(x,...){
  X <- list()
  tempFold <- x$fold
  fix <- x$fix
  precPos <- c(3,7,11,16)
  for(i in 1:tempFold){
    temp <- as.numeric(unlist(x$mdr[precPos[i]])[1:x$top])
    if(i==1){
      if(fix!=-1) temp <- temp[x$top]
      temp2 <- colnames(x$X)[unlist(x$mdr[precPos[i]+1])[1:x$top]]
      if(fix!=-1) temp2 <- temp2[1]
      
      if(length(x$res.sampled)>0){
        temp.cvc <- as.numeric(unlist(x$res.sampled[precPos[i]])[1:x$top])
        temp2.cvc <- list()
        for(cvcRun in 1:length(x$res.sampled)){
          if(fix!=-1) temp <- temp[x$top]
          temp2.cvc[[cvcRun]] <- colnames(x$X)[unlist(x$res.sampled[precPos[i]+1])[1:x$top]]
          if(fix!=-1) temp2 <- temp2[1]
        }
      }
    } else if(i==2){
      tempPos1 <- unlist(x$mdr[precPos[i]+1])[1:x$top]
      tempPos2 <- unlist(x$mdr[precPos[i]+2])[1:x$top]
      temp2 <- paste(colnames(x$X)[tempPos1],colnames(x$X)[tempPos2],sep=",")
      if(length(x$res.sampled)>0){
        temp2.cvc <- list()
        for(cvcRun in 1:length(x$res.sampled)){
          tempPos1 <- unlist(x$res.sampled[[cvcRun]][precPos[i]+1])[1:x$top]
          tempPos2 <- unlist(x$res.sampled[[cvcRun]][precPos[i]+2])[1:x$top]
          temp2.cvc[[cvcRun]] <- paste(colnames(x$X)[tempPos1],colnames(x$X)[tempPos2],sep=",")
        }
      }
    } else if(i==3){
      tempPos1 <- unlist(x$mdr[precPos[i]+1])[1:x$top]
      tempPos2 <- unlist(x$mdr[precPos[i]+2])[1:x$top]
      tempPos3 <- unlist(x$mdr[precPos[i]+3])[1:x$top]
      temp2 <- paste(colnames(x$X)[tempPos1],colnames(x$X)[tempPos2],colnames(x$X)[tempPos3],sep=",")
      if(length(x$res.sampled)>0){
        temp2.cvc <- list()
        for(cvcRun in 1:length(x$res.sampled)){
          tempPos1 <- unlist(x$res.sampled[[cvcRun]][precPos[i]+1])[1:x$top]
          tempPos2 <- unlist(x$res.sampled[[cvcRun]][precPos[i]+2])[1:x$top]
          tempPos3 <- unlist(x$res.sampled[[cvcRun]][precPos[i]+3])[1:x$top]
          temp2.cvc[[cvcRun]] <- paste(colnames(x$X)[tempPos1],colnames(x$X)[tempPos2],colnames(x$X)[tempPos3],sep=",")
        }
      }
    } else if(i==4){
      tempPos1 <- unlist(x$mdr[precPos[i]+1])[1:x$top]
      tempPos2 <- unlist(x$mdr[precPos[i]+2])[1:x$top]
      tempPos3 <- unlist(x$mdr[precPos[i]+3])[1:x$top]
      tempPos4 <- unlist(x$mdr[precPos[i]+4])[1:x$top]
      temp2 <- paste(colnames(x$X)[tempPos1],colnames(x$X)[tempPos2],colnames(x$X)[tempPos3],colnames(x$X)[tempPos4],sep=",")
      if(length(x$res.sampled)>0){
        temp2.cvc <- list()
        for(cvcRun in 1:length(x$res.sampled)){
          tempPos1 <- unlist(x$res.sampled[[cvcRun]][precPos[i]+1])[1:x$top]
          tempPos2 <- unlist(x$res.sampled[[cvcRun]][precPos[i]+2])[1:x$top]
          tempPos3 <- unlist(x$res.sampled[[cvcRun]][precPos[i]+3])[1:x$top]
          tempPos4 <- unlist(x$res.sampled[[cvcRun]][precPos[i]+4])[1:x$top]
          temp2.cvc[[cvcRun]] <- paste(colnames(x$X)[tempPos1],colnames(x$X)[tempPos2],colnames(x$X)[tempPos3],colnames(x$X)[tempPos4],sep=",")
        }
      }
    }
    if(length(x$res.sampled)>0){
      cvc <- rep(0,length(temp2))
      for(resRun in 1:length(temp2)){
        for(cvcRun in 1:length(x$res.sampled)){
          cvc[resRun] <- cvc[resRun] + sum(is.element(temp2[resRun], temp2.cvc[[cvcRun]]))
        }
      }
      cvc <- cvc/length(x$res.sampled)
      tempDF <- data.frame(SNP=temp2,Acc=temp, CVC=cvc)
    } else {
      tempDF <- data.frame(SNP=temp2,Acc=temp)      
    }

    X[[i]] <- tempDF[rev(rownames(tempDF)),]
    names(X)[i] <- paste("Fold",i,sep="")
    if((i==1) && (fix!=0)){
      
    }else {
      rownames(X[[i]]) <- 1:x$top
    }
  }
#  X$cvRes <- x$cvRes
  print(X,...)
} 
