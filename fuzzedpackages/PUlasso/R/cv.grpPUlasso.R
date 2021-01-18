#' Cross-validation for PUlasso
#'
#' Do a n-fold cross-validation for PUlasso.
#'
#'@importFrom Rcpp evalCpp
#'@import methods
#'@import Matrix
#'@import parallel
#'@import doParallel
#'@import foreach
#'@importFrom stats sd
#'@useDynLib PUlasso
#'@param X Input matrix; each row is an observation. Can be a matrix or a sparse matrix.
#'@param z Response vector representing whether an observation is labeled or unlabeled.
#'@param py1 True prevalence Pr(Y=1)
#'@param initial_coef A vector representing an initial point where we start PUlasso algorithm from.
#'@param group A vector representing grouping of the coefficients. For the least ambiguity, it is recommended if group is provided in the form of vector of consecutive ascending integers.
#'@param penalty penalty to be applied to the model. Default is sqrt(group size) for each of the group.
#'@param lambda A user supplied sequence of lambda values. If unspecified, the function automatically generates its own lambda sequence based on nlambda and lambdaMinRatio.
#'@param nlambda The number of lambda values.
#'@param lambdaMinRatio Smallest value for lambda, as a fraction of lambda.max which leads to the intercept only model.
#'@param maxit Maximum number of iterations.
#'@param weights observation weights. Default is 1 for each observation.
#'@param eps Convergence threshold for the outer loop. The algorithm iterates until the maximum change in coefficients is less than eps in the outer loop.
#'@param inner_eps Convergence threshold for the inner loop. The algorithm iterates until the maximum change in coefficients is less than eps in the inner loop.
#'@param verbose A logical value. if TRUE, the function prints out the fitting process.
#'@param stepSize A step size for gradient-based optimization. if NULL, a step size is taken to be stepSizeAdj/mean(Li) where Li is a Lipschitz constant for ith sample
#'@param stepSizeAdjustment A step size adjustment. By default, adjustment is 1 for GD and SGD, 1/8 for SVRG and 1/16 for SAG.
#'@param batchSize A batch size. Default is 1.
#'@param updateFrequency An update frequency of full gradient for method =="SVRG"
#'@param samplingProbabilities sampling probabilities for each of samples for stochastic gradient-based optimization. if NULL, each sample is chosen proportionally to Li.
#'@param method Optimization method. Default is Coordinate Descent. CD for Coordinate Descent, GD for Gradient Descent, SGD for Stochastic Gradient Descent, SVRG for Stochastic Variance Reduction Gradient, SAG for Stochastic Averaging Gradient.
#'@param trace An option for saving intermediate quantities when fitting a full dataset.
#'@param nfolds Number of cross-validation folds to be created.
#'@param fitInd A vector of indices of cross-validation models which will be fitted. Default is to fit the model for each of the cross-validation fold.
#'@param nCores Number of threads to be used for parallel computing. If nCores=0, it is set to be (the number of processors available-1) . Default value is 1.
#'@return cvm Mean cross-validation error
#'@return cvsd Estimate of standard error of cvm
#'@return cvcoef Coefficients for each of the fitted CV models
#'@return cvstdcoef Coefficients in a standardized scale for each of the fitted CV models
#'@return lambda The actual sequence of lambda values used.
#'@return lambda.min Value of lambda that gives minimum cvm.
#'@return lambda.1se The largest value of lambda such that the error is within 1 standard error of the minimum cvm.
#'@return PUfit A fitted PUfit object for the full data
#'@examples
#'data("simulPU")
#'fit<-cv.grpPUlasso(X=simulPU$X,z=simulPU$z,py1=simulPU$truePY1)
#'@export
#'
cv.grpPUlasso <-function(X,z,py1,initial_coef=NULL,group=1:p,
                         penalty=NULL,lambda=NULL, nlambda = 100,
                         lambdaMinRatio=ifelse(N< p, 0.05, 0.005),
                         maxit=ifelse(method=="CD",1000,N*10),
                         weights = NULL,eps=1e-04,inner_eps = 1e-02,
                         verbose = FALSE, stepSize=NULL, stepSizeAdjustment = NULL, 
                         batchSize=1, updateFrequency = N,
                         samplingProbabilities=NULL, method=c("CD","GD","SGD","SVRG","SAG"),
                         nfolds=10,fitInd=1:nfolds,nCores=1,trace=c("none","param","fVal","all"))
{
  N = nrow(X); p = ncol(X)
  if(length(lambda)==1){stop("More than one lambda needed for cross-validation")}
  if(!any(fitInd%in%(1:nfolds))){stop("an element of fitInd should be an integer between 1 and nfolds")
  }else{nfits= length(fitInd)}
  if(nfits<2){stop("length(fitInd) should be at least two")}
  
  # Match arguments
  method = match.arg(method,choices=c("CD","GD","SGD","SVRG","SAG"))
  trace = match.arg(trace,choices=c("none","param","fVal","all"))
  
  if(verbose){cat("Run grpPUlasso on full dataset\n")}
  fit.full<- grpPUlasso(X = X,z = z,py1 = py1,initial_coef = initial_coef,
              group = group,penalty = penalty,lambda = lambda,
              nlambda = nlambda,lambdaMinRatio = lambdaMinRatio,
              maxit = maxit,weights = weights,eps = eps,
              inner_eps = inner_eps,verbose = verbose,
              stepSize = stepSize,stepSizeAdjustment = stepSizeAdjustment,
              batchSize = batchSize,updateFrequency = updateFrequency,
              samplingProbabilities = samplingProbabilities,
              method = method,
              trace = trace)
  
  lambdaseq = fit.full$lambda
  convFlag_f = fit.full$optResult$convergence
  ############################################################################
  ## CV
  ############################################################################
  if(verbose){cat("Start Cross-Validation\n")}
  
  # Reorder samples of X,z
  row_ordering= order(z,decreasing = T); 
  nl = sum(z); nu = length(z) -nl
  folds = createFolds(row_ordering,nl,nu,nfolds)
  
  if (nCores<1||nCores>1){isParallel=TRUE}else{isParallel=FALSE}

  if (isParallel) {

    if(nCores==0){nCores = max(1, detectCores() - 1)
    }else{
      nCores = min(nCores,detectCores())
    }
    cl <- makeCluster(nCores)
    registerDoParallel(cl)
    clusterCall(cl, function(x) .libPaths(x), .libPaths())
    if(verbose){cat('Cross-Validation with',getDoParWorkers(), 'workers\n')}

    g=foreach(l=1:nfits,
              .packages = "PUlasso",
              .combine = list,
              .multicombine = TRUE)  %dopar%
              {
                k = fitInd[l]
                train_X = X[folds!=k,,drop=F]
                train_z = z[folds!=k]
                train_weights = weights[folds!=k]
                
                fit_cv = grpPUlasso(X = train_X,z = train_z,py1 = py1,initial_coef = initial_coef,
                                     group = group,penalty = penalty,lambda = lambdaseq,
                                     nlambda = nlambda,lambdaMinRatio = lambdaMinRatio,
                                     maxit = maxit,weights = train_weights,eps = eps,
                                     inner_eps = inner_eps,verbose = verbose,
                                     stepSize = stepSize,stepSizeAdjustment = stepSizeAdjustment,
                                     batchSize = batchSize,updateFrequency = updateFrequency,
                                     samplingProbabilities = samplingProbabilities,
                                     method = method,
                                     trace = trace)
                fit_cv
              }

    cvdev = foreach(l=1:nfits,
                    .packages = "PUlasso",
                    .combine  = "cbind" )%dopar%
                    {
                      k = fitInd[l]
                      test_X = X[folds==k,,drop=F]
                      test_z = z[folds==k]
                      test_weights = weights[folds==k]
                      
                      cvdev<- deviances(X = test_X,z = test_z,
                                        py1 = pi,coefMat = g[[l]]$coef,
                                        weights = test_weights)
                      return(cvdev)
                    }
    stopCluster(cl)

  } else {
    g=list()
    cvdev = matrix(NA,ncol=nfits,nrow=length(lambdaseq))
    for(l in 1:nfits){
      k = fitInd[l]
      if(verbose){cat('Cross-Validation for dataset',k,'\n')}
      
      train_X = X[folds!=k,,drop=F]
      train_z = z[folds!=k]
      train_weights = weights[folds!=k]
      
      test_X = X[folds==k,,drop=F]
      test_z = z[folds==k]
      test_weights = weights[folds==k]

      g[[l]] = grpPUlasso(X = train_X,z = train_z,py1 = py1,initial_coef = initial_coef,
                           group = group,penalty = penalty,lambda = lambdaseq,
                           nlambda = nlambda,lambdaMinRatio = lambdaMinRatio,
                           maxit = maxit,weights = train_weights,eps = eps,
                           inner_eps = inner_eps,verbose = verbose,
                           stepSize = stepSize,stepSizeAdjustment = stepSizeAdjustment,
                           batchSize = batchSize,updateFrequency = updateFrequency,
                           samplingProbabilities = samplingProbabilities,
                           method = method,
                           trace = trace)
      
      cvdev[,l] <- deviances(X = test_X,z = test_z,py1 = pi,
                             coefMat = g[[l]]$coef,
                             weights = test_weights)
    }
  }# End of Fitting
  # Summary
  coefmat <- list()
  std_coefmat <- list()
  for (i in 1:min(nfolds,nfits)){
    coefmat[[i]] <- g[[i]]$coef
    std_coefmat[[i]] <- g[[i]]$std_coef
  }
  names(coefmat) <- paste("cv",fitInd,sep="")
  names(std_coefmat) <- paste("cv",fitInd,sep="")

  # cvdev=sapply(g,function(x){x$deviance})
  # rownames(cvdev)=paste("l",1:length(lambdaseq),sep = "")
  cvm=apply(cvdev,1,mean)
  cvsd <- apply(cvdev,1,sd)/sqrt(min(nfolds,nfits))
  names(cvm)=paste("l",1:length(lambdaseq),sep = "")
  names(cvsd)=paste("l",1:length(lambdaseq),sep = "")
  indmin <- min(which(cvm==min(cvm)))
  lambda.min <- lambdaseq[indmin]

  ind <-  intersect(which(cvm>=cvm[indmin]+cvsd[indmin]),(1:indmin))
  if(length(ind)==0){ind1se <-  indmin
  } else {
    ind1se <-  max(ind)
  }
  lambda.1se <- lambdaseq[ind1se]

  convFlagMat=sapply(g,function(x){x$optResult$convergence})

  # Warning
  if(method %in% c("CD","GD")){
    widx<-which(convFlag_f==1)
    if(length(widx)>0){
      for(i in 1:length(widx)){
        warning(paste("convergence failed at ",widx[i],"th lambda, ", fit.full$iters[widx[i]],"th iterations",sep=""))
      }
    }

    for(j in 1:min(nfolds,nfits)){
      widx<-which(convFlagMat[,j]==1)
      if(length(widx)>0){
        for(i in 1:length(widx)){
          warning(paste("cvset",fitInd[j]," convergence failed at ",widx[i],"th lambda",sep=""))}}
    }
  }else{
    if(verbose){
      widx<-which(convFlag_f==0)
      if(length(widx)>0){
        for(i in 1:length(widx)){
          cat('|param.diff| < eps at',widx[i],'th lambda,',fit.full$iters[widx[i]],'th iterations\n')
        }
      }

      for(j in 1:min(nfolds,nfits)){
        widx<-which(convFlagMat[,j]==0)
        if(length(widx)>0){
          for(i in 1:length(widx)){
            cat('cvset',j,'|param.diff| < eps at',widx[i],'th lambda\n')}}
      }
    }
  }

  result<-structure(list(cvm=cvm,cvsd=cvsd, 
                         cvcoef = coefmat, 
                         cvstdcoef = std_coefmat, 
                         lambda = lambdaseq, 
                         lambda.min= lambda.min,
                         lambda.1se= lambda.1se,
                         PUfit=fit.full,
                         cvfolds = folds,
                         call=match.call()),class="cvPUfit")

  return(result)
}
