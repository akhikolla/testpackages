#' Print function for GPCMlasso
#' 
#' Print function for a \code{GPCMlasso} object. Prints parameters estimates for all model
#' components for the optimal model chosen by a specific criterion (by default BIC). 
#' 
#' @param x \code{GPCMlasso} object
#' @param select Specifies which criterion to use for the optimal model, we recommend the 
#' default value "BIC". If cross-validation was performed, automatically the optimal
#' model according to cross-validation is used. Only the parameter estimates from
#' the chosen optimal model are printed.
#' @param ... Further print arguments.
#' @author Gunther Schauberger\cr \email{gunther.schauberger@@tum.de}
#' @references Schauberger, Gunther and Mair, Patrick (2019): A Regularization Approach for the Detection of Differential 
#' Item Functioning in Generalized Partial Credit Models, \emph{Behavior Research Methods}, \url{https://link.springer.com/article/10.3758/s13428-019-01224-2}
#' @seealso \code{\link{GPCMlasso}}
#' @examples
#' data(tenseness_small)
#' 
#' ## formula for simple model without covariates
#' form.0 <- as.formula(paste("cbind(",paste(colnames(tenseness_small)[1:5],collapse=","),")~0"))
#' 
#' ######
#' ## fit simple RSM where loglikelihood and score function are evaluated parallel on 2 cores
#' rsm.0 <- GPCMlasso(form.0, tenseness_small, model = "RSM", 
#' control= ctrl_GPCMlasso(cores=2))
#' rsm.0
#' 
#' \dontrun{
#' ## formula for model with covariates (and DIF detection)
#' form <- as.formula(paste("cbind(",paste(colnames(tenseness_small)[1:5],collapse=","),")~."))
#' 
#' ######
#' ## fit GPCM model with 10 different tuning parameters
#' gpcm <- GPCMlasso(form, tenseness_small, model = "GPCM", 
#'                   control = ctrl_GPCMlasso(l.lambda = 10))
#' gpcm
#' plot(gpcm)
#' pred.gpcm <- predict(gpcm)
#' trait.gpcm <- trait.posterior(gpcm)
#' 
#' ######
#' ## fit RSM, detect differential step functioning (DSF)
#' rsm.DSF <- GPCMlasso(form, tenseness_small, model = "RSM", DSF = TRUE, 
#'                      control = ctrl_GPCMlasso(l.lambda = 10))
#' rsm.DSF
#' plot(rsm.DSF)
#' 
#' ## create binary data set
#' tenseness_small_binary <- tenseness_small
#' tenseness_small_binary[,1:5][tenseness_small[,1:5]>1] <- 2
#' 
#' ######
#' ## fit and cross-validate Rasch model
#' set.seed(1860)
#' rm.cv <- GPCMlasso(form, tenseness_small_binary, model = "RM", cv = TRUE, 
#'                    control = ctrl_GPCMlasso(l.lambda = 10))
#' rm.cv
#' plot(rm.cv)
#' }
print.GPCMlasso <- function(x, select = c("BIC", "AIC", "cAIC", "cv"), ...){
  with(x$design_list,{

    select <- match.arg(select, c("BIC", "AIC", "cAIC", "cv"))

    if(select=="BIC"){
      criterion <- x$BIC
    }
    if(select=="AIC"){
      criterion <- x$AIC
    }
    if(select=="cAIC"){
      criterion <- x$cAIC
    }
    if(select=="cv" | !is.null(x$cv_error)){
      criterion <- x$cv_error
    }
    
      coefs <- x$coefficients[which.min(criterion),]
      loglik <- x$logLik[which.min(criterion)]


  if(RSM){
    delta <- coefs[1:I]
    names(delta) <- x$item.names
    alpha <- c(0,coefs[(I+1):(I+q[1]-1)])
    names(alpha) <- paste0("Catgr. ",1:(q[1]))
    start.gamma <- start.sigma <- I+q[1]
  }else{
    if(x$model %in% c("RM","2PL")){
      delta <- coefs[1:I]
      names(delta) <- x$item.names
      alpha <-NA
      start.gamma <- start.sigma <- I+1
    }else{
      delta <- coefs[1:sum(q)]
      if(sum(diff(q))==0){
        delta <- matrix(delta, byrow = TRUE, ncol = q[1])
        rownames(delta) <- x$item.names
        colnames(delta) <- paste0("Catgr. ",1:(q[1]))
      }else{
        delta <- split(delta,rep(1:I,q))
        names(delta) <- x$item.names
        for(ii in 1:I){
          names(delta[[ii]]) <- paste0("Catgr. ",1:(q[ii]))
        }
      }
      alpha <- NA
      start.gamma <- start.sigma <- sum(q)+1
    }
  }
  
  if(x$main.effects & m > 0){
    mains <- coefs[start.gamma:(start.gamma+m-1)]
    names(mains) <- x.names
    start.gamma <- start.sigma <- start.gamma + m
  }else{
    mains <- NA
  }    
      
  n.dif.par <- I*m
  if(x$DSF){
    n.dif.par <- sum(q)*m
  }

  if(m>0){
    if(x$DSF){
      if(sum(diff(q))==0){
        gamma <- matrix(coefs[start.gamma:(start.gamma + n.dif.par - 1)], byrow = TRUE, nrow = I)
        rownames(gamma) <- x$item.names
        colnames(gamma) <- paste(rep(x.names, q[1]), rep(1:q[1], each = length(x.names)), sep = ".")
        index.gamma <- c(t(matrix(1:(m*q[1]), nrow = m)))
        gamma <- gamma[, index.gamma]
      }else{
        gamma <- list()
        coefs.pen <- coefs[start.gamma:(start.gamma+n.dif.par-1)]
        start.coefs <- 1
        for(i in 1:I){
          coefs.i <- coefs.pen[start.coefs:(start.coefs-1+q[i]*m)]
          gamma.index <- matrix(1:(q[i]*m),nrow=m)
          gamma.i <- c()
          for(ii in 1:m){
            coefs.ii <- coefs.i[gamma.index[ii,]]
            names(coefs.ii) <- paste0(x.names[ii],".",1:(q[i]))
            gamma.i <- c(gamma.i,coefs.ii)
          }
          gamma[[i]] <- gamma.i
          start.coefs <- q[i]*m+start.coefs
        }
      names(gamma) <- x$item.names
      }
    }else{
      gamma <- matrix(coefs[start.gamma:(start.gamma+n.dif.par-1)],byrow=TRUE,nrow=I)
      rownames(gamma) <- x$item.names
      colnames(gamma) <- x.names
    }
    start.sigma <- start.gamma+n.dif.par
  }else{
    gamma <- NA
  }
  
  if(GPCM){
    sigma <- coefs[start.sigma:(start.sigma+I-1)]
    names(sigma) <- x$item.names
  }else{
    sigma <- coefs[start.sigma]
  }
  
  cat("Call: ")
  print(x$call,...)
  cat("\n")
  
  cat("Item parameters:","\n")
  print(delta, ...)
  
  
  
  suppressWarnings(if(!is.na(alpha[1])){
    cat("\n")
    cat("Category parameters:","\n")
    print(alpha, ...)
  cat("\n")  })
  
  suppressWarnings(if(!is.na(mains[1])){
    cat("\n")
    cat("Main covariate effects:","\n")
    print(mains, ...)
  cat("\n")  })
  
  
  suppressWarnings(if(!is.na(gamma[1])){
    cat("\n")
    cat("DIF parameters:","\n")
    print(gamma, ...)
  })
  cat("\n")
  
  cat("Discrimination parameter(s):","\n")
  cat(sigma,"\n")
  
  cat("\n")
  
  cat("Log-likelihood:","\n")
  cat(loglik)
  cat("\n")
  
  invisible(list(x = x, gamma = gamma))
  })

}