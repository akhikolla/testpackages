
grp_as_numeric<-function(group){
  if(!is.numeric(group)){
    group <-  as.factor(group)
    levels(group) <-1:length(unique(group))
    group <- as.numeric(levels(group)[group])
  }
  return(group)
}

ordering_data<-function(row_ordering,col_ordering,
                        X, z, group, weights){
  X_lu <- X[row_ordering,col_ordering,drop=F]
  group0 <- group; group <- grp_as_numeric(group)
  group0 <- group0[col_ordering]; 
  group <- group[col_ordering]
  z_lu <- z[row_ordering]
  if(!is.null(weights)){w_lu <- weights[row_ordering]}else{w_lu = NULL}
  
  return(list(X_lu =  X_lu,z_lu =  z_lu,w_lu =  w_lu,group = group,group0 = group0))
}


fitting_setup<-function(py1,lambda,lambdaMinRatio,nlambda,initial_coef,group,penalty,p){
  
  # Setup lambdaseq, initial_coef, gsize, pen, weiOption
  # Return user_lanbdaseq, lambdaseq, initial_coef, gsize, pen, weiOption
  if (is.null(lambda)) {
    if (lambdaMinRatio >= 1){stop("lambdaMinRatio should be less than 1")}
    if (nlambda < 1){stop("nlambda should be at least 1")}
    user_lambdaseq <- FALSE
    lambdaseq <- c(0.1,0.01) # will not be used
  } else {
    if (any(lambda < 0)){stop("lambdas should be non-negative")}
    user_lambdaseq <- TRUE
    lambdaseq <- sort(lambda, decreasing = TRUE)
  }
  
  if(is.null(initial_coef)){
    icoef <- rep(0,p+1); icoef[1] = log(py1/(1-py1))
    if(is.nan(icoef[1])){stop("not a valid py1=P(Y=1)")}
  }else{
    if(length(initial_coef)!=(p+1)){stop("length of initial_coef should be the same as ncol(X)+1")}
    icoef <- initial_coef
  }
  
  gsize <-  c(1,table(group))
  if(is.null(penalty)){
    pen <- c(0,rep(1,length(gsize)-1))*sqrt(gsize)
  } else{
    pen <- c(0, penalty)
  }
  
  return(list( user_lambdaseq = user_lambdaseq, lambdaseq = lambdaseq, 
               icoef = icoef, gsize = gsize, pen = pen))
}

opt_option_setup<-function(method, trace, stepSize, stepSizeAdjustment,samplingProbabilities){
  
  trace = switch(trace,"none"=0,"param"=1,"fVal"=2,"all"=3)
  
  if(is.null(stepSize)||is.null(samplingProbabilities)){use_Lipschitz_for_ss_or_sProb = TRUE}
  if(is.null(stepSize)){stepSize=0}
  if(is.null(samplingProbabilities)){samplingProbabilities=1}
  if(is.null(stepSizeAdjustment)){
    if(method %in% c("GD","SGD")){
      adj = 1
    }else if(method =="SVRG"){
      adj = 1/8
    }else{
      adj = 1/16
    }
  }else{adj = stepSizeAdjustment}
  
  return(list(trace = trace,
              use_Lipschitz_for_ss_or_sProb = use_Lipschitz_for_ss_or_sProb,
              stepSize = stepSize,
              stepSizeAdjustment= adj,
              samplingProbabilities = samplingProbabilities))
}

summary_cpp_results<-function(g,method,trace,colnames,group0){
  
  coef <-  g$coef
  colnames(coef) <-  paste("l",1:length(g$lambda),sep = "")
  rownames(coef) <- c("(Intercept)",colnames)
  
  std_coef <- g$std_coef
  colnames(std_coef) <-  paste("l",1:length(g$lambda),sep = "")
  rownames(std_coef) <- c("(Intercept)",paste("group",group0))
  
  if(method=="CD"){
    iters = g$nUpdates
  }else{iters=g$iters}
  
  p = nrow(coef)-1
  # when trace != "none"
  trace = switch(trace,"none"=0,"param"=1,"fVal"=2,"all"=3)
  std_coef_all = NULL
  fVals_all = NULL
  if(trace == 1){
    
    std_coef_all= list()
    for(k in 1:length(g$lambda)){
      std_coef_all[[k]] <- g$beta_all[((k-1)*(p+1)+1):(k*(p+1)),1:(iters[k]+1)]
      rownames(std_coef_all[[k]]) <- rownames(std_coef)
      colnames(std_coef_all[[k]]) <- 1:(iters[k]+1)
    }
    names(std_coef_all) <- paste("lambda",1:length(g$lambda),sep="")
    
  } else if(trace==2){
    fVals_all = g$fVals_all[1:max(iters),,drop=F]
    
  } else if(trace==3) {
    
    std_coef_all= list()
    for(k in 1:length(g$lambda)){
      std_coef_all[[k]] <- g$beta_all[((k-1)*(p+1)+1):(k*(p+1)),1:(iters[k]+1)]
      rownames(std_coef_all[[k]]) <- rownames(std_coef)
      colnames(std_coef_all[[k]]) <- 1:(iters[k]+1)
    }
    names(std_coef_all) <- paste("lambda",1:length(g$lambda),sep="")
    fVals_all = g$fVals_all[1:max(iters),,drop=F]
  }
  
  return(list(coef = coef, std_coef = std_coef, iters= iters, std_coef_all = std_coef_all, fVals_all = fVals_all))
}

input_check <- function(X,z,group,penalty,stepSize,samplingProbabilities,weights){
  if(is.null(dim(X))){stop("not a valid X")}
  if(nrow(X)!=length(z)){stop("nrow(X) must be same as length(z)")}
  if(!all(z%in%c(0,1))){stop("z should be 0 or 1")}
  if(length(group)!=ncol(X)){stop("lenght(group) should be the same as ncol(X)")}
  if(!is.null(penalty)){
    if(length(penalty)!=length(unique(group))){stop("length(penalty) should be the same as the group size")}
  }
  
  if(!is.null(stepSize)){
    if(stepSize <= 0){stop("step size should be > 0 ")}
  }
  if(!is.null(samplingProbabilities)){
    if(length(samplingProbabilities)!=nrow(X)){stop("length of sampling probability should be equal to the nrow(X)")}
    if(any(samplingProbabilities < 0)){stop("all sampling probabilities should >= 0")}
  }
  if(!is.null(weights)){
    if(any(weights < 0)){stop("all weights should be positive")}
    if(length(weights)!= nrow(X)){stop("length(w_lu) should be the same as nrow(X)")}}
}



createFolds <-function(row_ordering, nl,nu,nfolds){
  
  if(min(nl,nu)<nfolds){stop("min(nl,nu)<nfolds")}
  
  # Calculate each group size
  rmrdl=nl%%nfolds
  rmrdu=nu%%nfolds
  cvnl=floor(nl/nfolds)
  cvnu=floor(nu/nfolds)
  cvsizel<-c()
  cvsizeu<-c()
  for(i in 1:nfolds){
    cvsizel[i]<-ifelse(i<=rmrdl,cvnl+1,cvnl)
    cvsizeu[i]<-ifelse(i<=rmrdu,cvnu+1,cvnu)
  }
  
  labeled_samples_idx = row_ordering[1:nl]
  unlabeled_samples_idx = row_ordering[-(1:nl)]
  
  #idx: row_id in X matrix, assignments = test group
  assignMat=  rbind(data.frame(idx = labeled_samples_idx,assignments=sample(rep(1:nfolds,cvsizel))),
                    data.frame(idx = unlabeled_samples_idx, assignments=sample(rep(1:nfolds,cvsizeu))))
  
  assignMat= assignMat[order(assignMat$idx),]
  assignMat$assignments
}
