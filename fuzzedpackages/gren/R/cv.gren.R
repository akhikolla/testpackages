# cross-validation performance of group-regularized elastic net
cv.gren <- function(x, y, m=rep(1, nrow(x)), unpenalized=NULL, partitions=NULL, 
                    alpha=0.5, lambda=NULL, intercept=TRUE, monotone=NULL, 
                    psel=TRUE, compare=TRUE, posterior=FALSE, nfolds=nrow(x), 
                    foldid=NULL, trace=TRUE,
                    control=list(epsilon=0.001, maxit=500, maxit.opt=1000, 
                                 maxit.vb=100),
                    keep.pred=TRUE, fix.lambda=FALSE, nfolds.out=nrow(x), 
                    foldid.out=NULL,
                    type.measure=c("auc", "deviance", "class.error")) {
  
  # save argument list
  argum <- formals(gren)
  fit.call <- match.call()
  
  # change some input for convenience
  if(is.data.frame(x)) {
    names.xr <- colnames(x)
    x <- as.matrix(x)
  } else {
    names.xr <- NULL
  }
  if(is.vector(partitions) & is.atomic(partitions)) {
    partitions <- list(partition1=partitions)
  }
  
  # save argument list and call
  argum <- formals(gren)
  fit.call <- match.call()
  
  # change some input for convenience
  if(is.data.frame(x)) {
    names.xr <- colnames(x)
    x <- as.matrix(x)
  } else {
    names.xr <- NULL
  }
  if(is.vector(partitions) & is.atomic(partitions)) {
    partitions <- list(partition1=partitions)
  }
  
  # check input
  if(!is.numeric(x) | !(is.numeric(y) | is.factor(y))) {
    stop("only numerical input data is supported at the moment")
  } else if(!is.null(ncol(y))) {
    if(ncol(y)!=2) {
      stop("y is either a vector or matrix with two columns")
    }
  } else if(ifelse(is.null(ncol(y)), length(y), nrow(y))!=length(m) | 
            nrow(x)!=length(m)) {
    stop("number of observations in y, m, and x not equal")
  } else if(!is.null(unpenalized) & !is.numeric(unpenalized) & 
            !is.data.frame(unpenalized)) {
    stop("only unpenalized data in a matrix or data.frame is supported")
  } else if(ifelse(is.null(ncol(unpenalized)), length(unpenalized), 
                   nrow(unpenalized))!=nrow(x) & !is.null(unpenalized)) {
    stop("number of observations in unpenalized not equal to number of 
         observations in y, x, and m")
  } else if(!is.list(partitions)) {
    stop("partitions should be either a list of partitions or one partition as 
         a numeric vector")
  } else if(any(lapply(partitions, length)!=ncol(x))) {
    stop("all partitions should be vectors of length ncol(x), containing the
         group identifiers of the features")
  } else if(!is.numeric(alpha) | length(alpha)!=1 | (alpha < 0) | (alpha > 1)) {
    stop("alpha should be of length one and a value between 0 and 1")
  } else if(!is.null(lambda) & !(is.vector(lambda) & is.atomic(lambda))) {
    stop("lambda should be either NULL or a numeric vector")
  } else if(!is.logical(intercept)) {
    stop("logical should be either TRUE or FALSE")
  } else if(!is.null(monotone) & 
            !(is.logical(monotone) & length(monotone)==length(partitions)) &
            !(is.list(monotone) & length(monotone)==2 & 
              is.logical(unlist(monotone)) & 
              length(unlist(monotone))==2*length(partitions))) {
    stop("monotone should be a vector of TRUEs and FALSEs of the same length
         as partitions or NULL")
  } else if(!is.null(psel) & !(is.vector(psel) & is.atomic(psel))) {
    stop("psel is either NULL or a vector of non-negative whole numbers")
  } else if(!is.logical(compare)) {
    stop("compare is either TRUE or FALSE")
  } else if(!is.logical(posterior)) {
    stop("posterior is either TRUE or FALSE")
  } else if(!is.numeric(nfolds) | length(nfolds)!=1) {
    stop("nfolds should be a whole number")
  } else if(!is.null(foldid) & !(is.numeric(foldid) & length(foldid==nrow(x)))) {
    stop("foldid is either NULL or a vector of length nrow(x)")
  } else if(!is.null(foldid) & !all(foldid %in% c(1:nrow(x)))) {
    stop("foldid must be a vector of length n, containing only whole numbers
         from 1 to nrow(x)")
  } else if(!is.logical(trace)) {
    stop("trace is either TRUE or FALSE")
  } else if(!is.list(control)) {
    stop("control must be a list")
  } else if(!is.numeric(control$epsilon) | length(control$epsilon)!=1) {
    stop("epsilon must be a non-negative number")
  } else if(!is.numeric(control$maxit) | length(control$maxit)!=1) {
    stop("maxit must be a non-negative whole number")
  } else if(!is.numeric(control$maxit.opt) | length(control$maxit.opt)!=1) {
    stop("maxit.opt must be a non-negative whole number")
  } else if(!is.numeric(control$maxit.vb) | length(control$maxit.vb)!=1) {
    stop("maxit.vb must be a non-negative whole number")
  } else if(!is.logical(keep.pred)) {
    stop("keep.pred should be either TRUE or FALSE")
  } else if(!is.logical(fix.lambda)) {
    stop("fix.lambda should be either TRUE or FALSE")
  } else if(!is.numeric(nfolds.out) | length(nfolds.out)!=1) {
    stop("nfolds.out should be a whole number")
  } else if(!is.null(foldid.out) & !(is.numeric(foldid.out) & 
                                     length(foldid.out==nrow(x)))) {
    stop("foldid.out is either NULL or a vector of length nrow(x)")
  } else if(!is.null(foldid.out) & !all(foldid.out %in% c(1:nrow(x)))) {
    stop("foldid.out must be a vector of length n, containing only whole numbers
         from 1 to nrow(x)")
  } else if(!is.character(type.measure)) {
    stop("type.measure should be one or more of: auc, deviance, class.error")
  } else if(!(all(type.measure %in% c("auc", "deviance", "class.error")))) {
    stop("type.measure should be one or more of: auc, deviance, class.error")
  } else if(any(m!=1) & any(type.measure %in% c("auc", "class.error"))) {
    stop("auc or misclassification error are currently only possible with 
         binary classification problems")
  }
  
  # extra check for data
  if(is.factor(y)) {
    y <- as.numeric(y) - 1
  }
  if(is.null(ncol(y))) {
    ymat <- cbind(m - y, y)
  } else {
    ymat <- y
    y <- y[, 2]
  }
  if(is.data.frame(unpenalized)) {
    unpenalized <- model.matrix(as.formula(paste("~", paste(
      colnames(unpenalized), collapse="+"))), data=unpenalized)[, -1]
  }
  u <- ifelse(is.null(unpenalized), 0, ncol(unpenalized))
  n <- nrow(x)
  r <- ncol(x)
  
  # determine the folds
  if(is.null(foldid.out)) {
    rest <- n %% nfolds.out
    foldsize <- c(rep(n %/% nfolds.out + as.numeric(rest!=0), times=rest),
                  rep(n %/% nfolds.out, times=nfolds.out - rest))
    foldid.out <- sample(rep(1:nfolds.out, times=foldsize))
  }
  nfolds.out <- length(unique(foldid.out))
  
  # if fix.lambda=TRUE, we estimate the global lambda only once
  if(fix.lambda) {
    if(trace) {cat("\r", "Estimating fixed global lambda by cross-validation", 
                   sep="")}
    srt <- proc.time()[3]
    cv.fit <- cv.glmnet(cbind(unpenalized, x), ymat, family="binomial", 
                        alpha=alpha, standardize=FALSE, intercept=intercept, 
                        foldid=foldid.out, grouped=FALSE,
                        penalty.factor=c(rep(0, u), rep(1, r)))
    lambda <- cv.fit$lambda.min
    cv.time <- proc.time()[3] - srt
    
    if(trace) {cat("\n", "Estimated fixed global lambda is ", round(lambda, 2), 
                   " in ", round(cv.time, 2), " seconds", sep="")
    }
  } else {
    lambda <- NULL
  }
  
  # set initial values to zero at beginning
  init <- cinit <- list(lambdag=NULL, mu=NULL, sigma=NULL, chi=NULL, ci=NULL)
  
  # create the predictions object based on the number of selected features
  if(!is.null(psel)) {
    if(is.numeric(psel)) {
      pred.groupreg <- matrix(NA, nrow=length(y), ncol=length(psel))
      if(compare) {
        pred.regular <- matrix(NA, nrow=length(y), ncol=length(psel))
      } else {
        pred.regular <- NULL
      }
    } else {
      pred.groupreg <- numeric(length(y))
      if(compare) {
        pred.regular <- numeric(length(y))
      } else {
        pred.regular <- NULL
      }
    }
  } else {
    pred.groupreg <- numeric(length(y))
    if(compare) {
      pred.regular <- numeric(length(y))
    } else {
      pred.regular <- NULL
    }
  }
  
  if(trace) {cat("\n", "Estimating performance by cross-validation", "\n", 
                 sep="")}
  srt <- proc.time()[3]
  cu <- u
  cpart <- partitions
  for(k in 1:nfolds.out) {
    cat("\r", "Fold ", k, sep="")
    # split the data into training and test data
    ytrain <- ymat[foldid.out!=k, ]
    xtrain <- matrix(x[foldid.out!=k, ], ncol=r)
    xtest <- matrix(x[foldid.out==k, ], ncol=r)
    
    # checking if any features are constant
    which.const1 <- apply(xtrain, 2, sd)==0
    if(any(which.const1)) {
      xtrain <- xtrain[, !which.const1]
      xtest <- xtest[, !which.const1]
      cpart <- lapply(partitions, function(l) {l[!which.const1]})
    } else {
      cpart <- partitions
    }
    if(is.null(unpenalized)) {
      utrain <- NULL
      utest <- NULL
      which.const2 <- NULL
    } else {
      utrain <- matrix(unpenalized[foldid.out!=k, ], ncol=u)
      utest <- matrix(unpenalized[foldid.out==k, ], ncol=u)
      which.const2 <- apply(utrain, 2, sd)==0
      if(any(which.const2)) {
        utrain <- utrain[, !which.const2]
        utest <- utest[, !which.const2]
      }
      cu <- ncol(utrain)
    }
    mtrain <- m[foldid.out!=k]
    mtest <- m[foldid.out==k]

    # fit the model on the training data
    if(k!=1) {
      if(intercept) {
        cinit$mu <- init$mu[!c(FALSE, which.const2, which.const1)]
        cinit$sigma <- init$sigma[!c(FALSE, which.const2, which.const1), 
                                  !c(FALSE, which.const2, which.const1)]
      } else {
        cinit$mu <- init$mu[!c(which.const2, which.const1)]
        cinit$sigma <- init$sigma[!c(which.const2, which.const1), 
                                  !c(which.const2, which.const1)]
      }
      cinit$chi <- init$chi[!which.const1]
      cinit$ci <- init$ci[foldid.out!=k]
    }
    fit.gren <- gren(xtrain, ytrain, mtrain, utrain, cpart, alpha, 
                     lambda, intercept, monotone, psel, compare, posterior, 
                     nfolds, foldid, FALSE, cinit, control)
    
    # predictions on the test data
    # if we set psel, we predict with the desired number of features
    if(!is.null(psel) & is.numeric(psel)) {
      for(csel in 1:length(psel)) {
        sel.lambda <- fit.gren$freq.model$groupreg$lambda[which.min(abs(
          psel[csel] - fit.gren$freq.model$groupreg$df - cu))[1]]
        pred.groupreg[foldid.out==k, ] <- predict(fit.gren$freq.model$groupreg,
                                                  cbind(utest, xtest),
                                                  s=sel.lambda, 
                                                  type="response")
          
      }
      if(compare) {
        sel.lambda <- fit.gren$freq.model$regular$lambda[which.min(abs(
          psel[csel] - fit.gren$freq.model$regular$df - cu))[1]]
        pred.regular[foldid.out==k, ] <- predict(fit.gren$freq.model$regular, 
                                                 cbind(utest, xtest),
                                                 s=sel.lambda, 
                                                 type="response")
      }
      pred.groupreg[foldid.out==k, ] <- predict(fit.gren$freq.model$groupreg, 
                                                cbind(utest, xtest),
                                                s=fit.gren$lambda, 
                                                type="response")
    } else {
      pred.groupreg[foldid.out==k] <- predict(fit.gren$freq.model$groupreg, 
                                              cbind(utest, xtest),
                                              s=fit.gren$lambda, 
                                              type="response")
      if(compare) {
        pred.regular[foldid.out==k] <- predict(fit.gren$freq.model$regular, 
                                               cbind(utest, xtest),
                                               s=fit.gren$lambda, 
                                               type="response")
      }
    }
    
    
    # update the initial values for the next iteration
    if(k==1) {
      init$mu <- rep(0, u + ncol(x) + intercept)
      init$sigma <- diag(u + ncol(x) + intercept)
      init$chi <- rep(1, ncol(x))
      init$ci <- rep(1, nrow(x))
    }
    if(intercept) {
      init$mu[!c(FALSE, which.const2, which.const1)] <- fit.gren$vb.post$mu
      init$sigma[!c(FALSE, which.const2, which.const1), 
                 !c(FALSE, which.const2, which.const1)] <- 
        fit.gren$vb.post$sigma
    } else {
      init$mu[!c(which.const2, which.const1)] <- fit.gren$vb.post$mu
      init$sigma[!c(which.const2, which.const1), 
                 !c(which.const2, which.const1)] <- fit.gren$vb.post$sigma
    }
    init$chi[!which.const1] <- fit.gren$vb.post$chi
    init$ci[foldid.out!=k] <- fit.gren$vb.post$ci
    init$lambdag <- fit.gren$lambdag
    
  }
  cv.time <- proc.time()[3] - srt
  if(trace) {cat("\n", "Cross-validated performance in ", round(cv.time, 2), 
                 " seconds", sep="")
  }
  
  # calculate requested performance measures
  out <- list(groupreg=list(pred=pred.groupreg), 
              regular=list(pred=pred.regular))
  if(any(type.measure=="auc")) {
    if(!is.null(psel) & is.numeric(psel)) {
      out$groupreg$auc <- apply(pred.groupreg, 2, function(pred) {
        pROC::roc(y, pred)$auc})
      if(compare) {
        out$regular$auc <- apply(pred.regular, 2, function(pred) {
          pROC::roc(y, pred)$auc})
      } 
    } else {
      out$groupreg$auc <- pROC::roc(y, pred.groupreg)$auc
      if(compare) {
        out$regular$auc <- pROC::roc(y, pred.regular)$auc
      }  
    }
  }
  if(any(type.measure=="deviance")) {
    if(!is.null(psel) & is.numeric(psel)) {
      out$groupreg$deviance <- colMeans((y==1)*log(pred.groupreg) + 
                                          (y==0)*log(1 - pred.groupreg))
      if(compare) {
        out$regular$deviance <- colMeans((y==1)*log(pred.regular) + 
                                           (y==0)*log(1 - pred.regular))
      }
    } else {
      out$groupreg$deviance <- mean((y==1)*log(pred.groupreg) + 
                                      (y==0)*log(1 - pred.groupreg))
      if(compare) {
        out$regular$deviance <- mean((y==1)*log(pred.regular) + 
                                       (y==0)*log(1 - pred.regular))
      }
    }
  }
  if(any(type.measure=="class.error")) {
    if(!is.null(psel) & is.numeric(psel)) {
      out$groupreg$class.error <- colMeans((pred.groupreg >= 0.5)==y)
      if(compare) {
        out$regular$class.error <- colMeans((pred.regular >= 0.5)==y)
      }
    } else {
      out$groupreg$class.error <- mean((pred.groupreg >= 0.5)==y)
      if(compare) {
        out$regular$class.error <- mean((pred.regular >= 0.5)==y)
      }
    }
  }
  
  return(out)
  
}