#' Predict function for GPCMlasso
#' 
#' Predict function for a \code{GPCMlasso} object. 
#' Predictions can be linear predictors or probabilities separately for each person and each item.
#' 
#' Results are lists of vectors with length equal to the number 
#' of response categories $k_i$ in case of
#' probabilities (\code{type="response"}) or 
#' $k_i-1$ in case of linear predictors (\code{type="link"}).
#' 
#' @param object \code{GPCMlasso} object
#' @param coefs Optional vector of coefficients, can be filled with a specific 
#' row from \code{object$coefficients}. If not specified, \code{coefs} are specififed to be
#' the BIC-optimal coefficients or, if cross-validation was performed, the optimal 
#' coefficients according to cross-validation.
#' @param newdata List possibly containing slots Y, X, Z1 and Z2 to use new data for prediction.
#' @param type Type "link" gives vectors of linear predictors for separate categories 
#' (of length $k_i-1$) and type "response" gives the respective probabilities (of length $k_i$).  
#' @param ... Further predict arguments.
#' @author Gunther Schauberger\cr \email{gunther@@stat.uni-muenchen.de}
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
predict.GPCMlasso <- function(object, coefs = NULL,
                              newdata = NULL,
                              type = c("link","response"),...){
  type <- match.arg(type)
  
  m <- object$design_list$m
  n <- object$design_list$n
  n_sigma <- object$design_list$n_sigma
  I <- object$design_list$I
  q <- object$design_list$q
  design <- object$design_list$design
  
  ## automatically get coefs if not provided
  if(is.null(coefs)){
    if(!is.null(object$cv_error)){
      cat("No coefs are provided, automatically cv-optimal model is chosen","\n")
      coefs <- object$coefficients[which.min(object$cv_error),]
    }else{
      cat("No coefs are provided, automatically BIC-optimal model is chosen","\n")
      coefs <- object$coefficients[which.min(object$BIC),]#
    }
  }
 
    

    ## extract covariates and response
  if(!is.null(newdata)){

  if(object$control$all.dummies){
    term.labels <- attr(terms(object$formula), "term.labels")
    X <- matrix(rep(1,nrow(newdata)))
    for(ij in 1:length(term.labels)){
      x.now <- newdata[,term.labels[ij], drop = FALSE]
      if(is.factor(x.now[,1])){
        if(nlevels(x.now[,1])==2){
          X <- cbind(X,model.matrix(~.,data=x.now)[,-1,drop=FALSE])
        }else{
          X <- cbind(X,model.matrix(~0+.,data=x.now))
        }
      }else{
        X <- cbind(X,model.matrix(~.,data=x.now)[,-1,drop=FALSE])
      }
    }
    X <- X[,-1,drop = FALSE]
  }else{
    X <- model.matrix(object$formula, data = newdata)
    if(ncol(X)>=1){
      if(colnames(X)[1]=="(Intercept)"){
        X <- X[,-1,drop = FALSE]
      }
    }
  }
    designX <- -get_designX(X, object$DSF, m, I, q, nrow(X))
  }else{
    designX <- object$design_list$designX
  }

    
    RSM <- GPCM <- FALSE
    if(object$model %in% c("GRSM","RSM")){
      RSM <- TRUE
    }
    if(object$model %in% c("GRSM","GPCM")){
      GPCM <- TRUE
    }
    
    des_total <- matrix(rep(t(design),n),byrow=TRUE,ncol=ncol(design))
    if(ncol(designX)>0){
      des_total <- cbind(des_total, designX)
    }
    
    coef_short <- head(coefs, length(coefs)-n_sigma)
    lin.pred <- des_total%*%coef_short
    
    theta <- trait.posterior(object, coefs, cores = object$control$cores)

    return_list <- rep(theta,each=sum(q))+lin.pred
    return_list <- split(return_list, rep(rep(1:I,q),n)+(rep(1:n,each=sum(q))-1)*I)
    
    if(type=="response"){
      for(i in 1:length(return_list)){
        pi.i <- responseFun(return_list[[i]])
        return_list[[i]] <- c(pi.i,1-sum(pi.i))
      }
    }
      
      names(return_list) <- paste(rep(paste0("Person_",1:n),each=I),rep(object$item.names,n),sep="_")
    return_list
    
  
}