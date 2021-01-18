fit_cv_GPCMlasso <- function(model, design_list,
                             control, score_fun, loglik_fun, log_score_fun, Y){
  
  with(design_list,{

        ## create cv folds
    fold_list <- get_fold_list(n=n, folds = control$folds)
    
    cat("Cross-Validation...", "\n")

    cv_result <- sapply(seq(control$folds), cv_fun, fold_list, Y,
                      control,  model, design_list,
                      loglik_fun, score_fun, log_score_fun)

   cv_error <- rowSums(cv_result)/design_list$n

  return(cv_error)  
    
  })
}


cv_fun <- function(ff, fold_list, Y, control, model,
                   design_list, loglik_fun, score_fun, log_score_fun) {

  I <- design_list$I
  n <- design_list$n
  n_sigma <- design_list$n_sigma
  q <- design_list$q
  
  scale_fac <- 1-1/control$folds
  
  if (control$trace) {
    cat("CV-fold:", ff, "out of", control$folds, "\n")
  }

  index <- matrix(1:(n*sum(q)),nrow=sum(q))
  index.test <- c(index[,fold_list[[ff]]])
  index.train <- c(index[,-fold_list[[ff]]])
  
  design.test <- design_list$designX[index.test,]
  design_list$designX <- design_list$designX[index.train,]
  
  response.test <- design_list$response[index.test]
  design_list$response <- design_list$response[index.train]
  
  design_list$n <- length(index.train)/sum(q)
  
  ## get fitted parameters 
  fit.fold <- fit_GPCMlasso(model = model, loglik_fun = loglik_fun, score_fun = score_fun, 
                            log_score_fun = log_score_fun,
                            design_list = design_list, 
                            control = control, start = NULL, scale_fac = scale_fac)
  coef.fold <- fit.fold$coefficients
  px <- ncol(coef.fold)
  
  criterion <- c()
  
  ## get nodes and weights for Gauss-Hermite quadrature
  her_poly <- gauss.quad(control$Q, "hermite")
  GHnodes <- her_poly$nodes
  GHweights <- her_poly$weights * exp(GHnodes^2) * dnorm(GHnodes)
  

  
  for(e in 1:nrow(coef.fold)){
  criterion[e] <-  loglik_fun(coef.fold[e,],
                     response.test,
                     design_list$design,
                     design.test,
                     control$Q,
                     q, length(fold_list[[ff]]), I, px,
                     GHweights,
                     GHnodes,
                     design_list$acoefs,
                     0,
                     0,
                     control$cvalue,
                     control$cores,
                     rep(1,ncol(design_list$acoefs)),
                     n_sigma,1)
  }
  criterion
}



get_fold_list <- function(n, folds){
  
  n.cv <- rep(floor(n/folds), folds)
  rest <- n%%folds
  if (rest > 0) {
    n.cv[1:rest] <- n.cv[1:rest] + 1
  }
  
  remaining_persons <- 1:n
  fold_list <- list()
  
  for(o in 1:folds){
    sample.o <- sample(remaining_persons[!is.na(remaining_persons)], n.cv[o])
    fold_list[[o]] <- sample.o
    remaining_persons[sample.o] <- NA
  }
  
  fold_list
}

