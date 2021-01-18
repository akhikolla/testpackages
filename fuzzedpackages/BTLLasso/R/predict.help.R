predict.help <- function(coefs, object, newdata , type){

  q <- object$Y$q
  
  n.theta <- object$design$n.theta

  if(length(newdata)==0){
    des.mat <- object$design$design.repar
    if(type=="trait"){
      des.first <- design.BTLLasso(Y = object$Y, X = object$X, Z1 = object$Z1,
                                   Z2 = object$Z2, control = object$control,
                                   only.first = TRUE)$design.repar
      
      des.mat <- des.mat[seq(1,nrow(des.mat)-1,by=q),]
      des.mat[,1:n.theta] <- 0 
      
      des.first <- des.first[seq(1,nrow(des.first)-1,by=q),]
      des.first[,1:n.theta] <- 0 
    }
  }else{
    des.mat <- design.BTLLasso(Y = newdata$Y, X = newdata$X, Z1 = newdata$Z1,
                               Z2 = newdata$Z2, control = object$control, 
                               sd.X = object$design$sd.X, sd.Z1 = object$design$sd.Z1,
                               sd.Z2 = object$design$sd.Z2)$design.repar
    if(type=="trait"){
      des.first <- design.BTLLasso(Y = newdata$Y, X = newdata$X, Z1 = newdata$Z1,
                                   Z2 = newdata$Z2, control = object$control,
                                   sd.X = object$design$sd.X, sd.Z1 = object$design$sd.Z1,
                                   sd.Z2 = object$design$sd.Z2, only.first = TRUE)$design.repar
      
      des.mat <- des.mat[seq(1,nrow(des.mat)-1,by=q),]
      des.mat[,1:n.theta] <- 0 
      
      des.first <- des.first[seq(1,nrow(des.first)-1,by=q),]
      des.first[,1:n.theta] <- 0 
    }
  }
  
  
    ncoef <- nrow(coefs)

  
  ret.list <- list()
  for(l in 1:ncoef){
    coefs.l <- coefs[l,]

    if(type!="trait"){
      eta <- des.mat%*%coefs.l
      ret.mat <- matrix(eta, byrow=TRUE, ncol=object$Y$q)
      
      if(type=="response"){
        ret.mat <- exp(ret.mat)/(1+exp(ret.mat))
      }
      
    }else{
      
      
      eta.both <- des.mat%*%coefs.l
      eta.first <- des.first%*%coefs.l
      eta.second <- eta.first-eta.both
      
      ret.mat <- cbind(eta.first, eta.second)
      colnames(ret.mat) <- c("first.object", "second.object")
    }
    ret.list[[l]] <- ret.mat
  }
  
  ret.list
}