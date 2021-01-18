
fit.cv.BTLLasso <- function(response, design, penalty, q, m, 
  folds = 10, lambda, control = ctrl.BTLLasso(), cores = folds, 
  trace = TRUE, trace.cv = TRUE, cv.crit) {
  
  
  
  k <- q + 1
  n.design <- nrow(design)/q
  
  
  if (trace.cv) {
    cat("Full model", "\n")
  }
  m.all <- fit.BTLLasso(response = response, design = design, 
    penalty = penalty, lambda = lambda, k = k, m = m, control = control, 
    trace = trace)
  
  
  ### cross validation
  
  n.cv <- rep(floor(n.design/folds), folds)
  rest <- n.design%%folds
  if (rest > 0) {
    n.cv[1:rest] <- n.cv[1:rest] + 1
  }
  
  which.fold <- rep(1:folds, n.cv)
  
  id.fold <- rep(sample(which.fold, n.design, replace = FALSE), 
    each = q)

  
  cv.fun <- function(ff) {

      if (trace.cv) {
      cat("CV-fold:", ff, "out of", folds, "\n")
    }
    
    design.train <- design[which(id.fold != ff), , drop = FALSE]
    design.test <- design[which(id.fold == ff), , drop = FALSE]
    
    if(any(apply(design.train,2,var)==0)){
      stop("In cross-validation one of the parameters is not estimable, 
probably because all correponding observations were eliminated from the training data.
Please change your seed and/or increase the number of folds!")
    }

    response.train <- response[which(id.fold != ff)]
    response.test <- response[which(id.fold == ff)]
    
    
    fit.fold <- fit.BTLLasso(response.train, design.train, 
      penalty = penalty, lambda = lambda, k = k, m = m, 
      control = control, trace = trace)

    coef.fold <- fit.fold$coefs
    
      if (cv.crit == "Deviance") {
      y.test <- t(cbind(matrix(response.test, ncol = q, 
        byrow = TRUE), 1)) * (1:k)
      y.test[y.test == 0] <- k + 1
      y.test <- apply(y.test, 2, min)
      
      yhelp <- rep(y.test, each = k)
      yhelp <- as.numeric(yhelp == rep(1:k, length(y.test)))
      
      preds <- c()
      for (u in 1:length(lambda)) {
        preds <- cbind(preds, predict.BTLLasso(coef.fold[u, 
          ], q, design.test))
      }
      
      criterion <- -2 * colSums(yhelp * log(preds))
    } else {
      pi.test <- c()
      for (u in 1:length(lambda)) {
        eta.test <- design.test %*% coef.fold[u, ]
        pi.test <- cbind(pi.test, exp(eta.test)/(1 + 
          exp(eta.test)))
      }
      criterion <- colSums((pi.test - response.test)^2)
    }
    
    criterion
  }
  
  cat("Cross-Validation...", "\n")
  if (cores > 1) {
    cl <- makeCluster(cores, outfile = "")
    
    clusterExport(cl, varlist = c("response", "design", "id.fold", 
      "lambda", "control", "trace.cv", "trace", "k", "m", 
      "cv.crit"), envir = sys.frame(sys.nframe()))
    
    criterion <- rowSums(parSapply(cl, seq(folds), cv.fun))
    stopCluster(cl)
  } else {
    criterion <- rowSums(sapply(seq(folds), cv.fun))
  }
  
  
  
  ret.list <- list(coefs = m.all$coefs, criterion = criterion)
  
  
  return(ret.list)
}
