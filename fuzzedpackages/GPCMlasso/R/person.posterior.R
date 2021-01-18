## from a vector of the linear predictor, get the respective
## probabilities for a adjacent categories model
responseFun <- function(eta) {
  q <- length(eta)
  eta[eta > 10] <- 10
  eta[eta < -10] <- -10
  eta.help <- matrix(rep(c(0, eta), each = q + 1), ncol = q + 
    1)
  eta.help[upper.tri(eta.help)] <- 0
  pi <- cumprod(c(1, exp(eta[-q])))/sum(apply(exp(eta.help), 
    1, prod))
  pi <- (pi - 0.5) * 0.99999 + 0.5
  pi
}


## Get likelihood contribution of one person, given an ability
## (theta), a vector conaining the linear predictor and the
## response vector.  For cross-validation, the test values are
## already eliminated.
likeli_i <- function(theta, lin.pred, sigma, yp, q) {
  
  lin.pred <- lin.pred + rep(sigma,q)*theta
  lin.pred <- split(lin.pred, rep(1:length(yp),q))
  
  likeli.i <- 1
  for (i in 1:length(lin.pred)) {
    mu_i <- responseFun(lin.pred[[i]])
    mu_i <- c(mu_i, 1 - sum(mu_i))
    likeli.i <- mu_i[yp[i]] * likeli.i
  }
  likeli.i
}

## Create vector of linear predictors for all observations of
## the model.  Linear predictor ALSO contains discrimination
## parameters.
lin_pred <- function(model, coef_short, sigma) {
  
  with(model$design_list, {
    ## create design matrix out ouf item-parameter part and
    ## covariate part

    des_total <- matrix(rep(t(design),n),byrow=TRUE,ncol=ncol(design))
    if (ncol(designX) > 0) {
      des_total <- cbind(des_total, designX)
    }
    
    ## multiply design matrix by discrimination parameters, so
    ## that discrimination parameters are included in linear
    ## predictor.
      des_total <- rep(rep(sigma,q), n) * des_total

    lin.pred <- des_total %*% coef_short
    
    return(lin.pred)
  })
  
}


## vectorized version (with respect to x, which is ability) of
## the function to be integrated for person abilities.
## Vectorization is favorable for integrate.
estim_i <- function(x, yp, lin.pred, sigma, q) {

  dens.norm <- dnorm(x)
  lik.i <- sapply(x, likeli_i, lin.pred, sigma, yp, q)
  ret <- dens.norm * lik.i * x
  return(ret)
}

## equivalent to estim_i except that x is not multiplied,
## necessary for scaling
factor_i <- function(x, yp, lin.pred, sigma, q) {
  
  dens.norm <- dnorm(x)
  lik.i <- sapply(x, likeli_i, lin.pred, sigma, yp, q)
  ret <- dens.norm * lik.i
  
  return(ret)
}

## Function to estimate person ability parameter.  
## integrate and pcubature are almost equally fast
person_i <- function(person, person.index, all.lin.preds, sigma, 
  Y, tol, limit, q) {

  yp <- Y[person, ]
  sigma <- sigma[!is.na(yp)]
  q <- q[!is.na(yp)]
  yp <- yp[!is.na(yp)]
  
  lin.pred <- all.lin.preds[person.index == person]
  
  est_i <- pcubature(estim_i, lowerLimit = -limit, upperLimit = limit, 
                     yp = yp, lin.pred = lin.pred, sigma = sigma, q = q, 
                     tol = tol)
  
  est_i <- est_i$integral/pcubature(factor_i, lowerLimit = -limit, 
                                 upperLimit = limit, yp = yp, lin.pred = lin.pred, sigma = sigma, 
                                 q = q, tol = tol)$integral
  
  
  est_i
}

#' Calculate Posterior Estimates for Trait Parameters
#' 
#' Calculates posterior estimates for trait/person parameters using the assumption 
#' of Gaussian distributed parameters.
#' 
#' 
#' @param model Object of class \code{GPCMlasso}. 
#' @param coefs Vector of coefficients to be used for prediction. If \code{coefs = NULL}, 
#' the parameters from  the BIC-optimal model will be used. 
#' If cross-validation was performed, automatically the parameters from the optimal
#' model according to cross-validation are used. 
#' @param cores Number of cores to be used in parallelized computation.
#' @param tol The maximum tolerance for numerical integration, 
#' for more details see \code{\link{pcubature}}.
#' @return 
#' \item{}{Vector containing all estimates of trait/person parameters.} 
#' @author Gunther Schauberger\cr \email{gunther.schauberger@@tum.de}
#' @references Schauberger, Gunther and Mair, Patrick (2019): A Regularization Approach for the Detection of Differential 
#' Item Functioning in Generalized Partial Credit Models, \emph{Behavior Research Methods}, \url{https://link.springer.com/article/10.3758/s13428-019-01224-2}
#' @seealso \code{\link{GPCMlasso}} \code{\link{GPCMlasso-package}}
#' @keywords GPCMlasso
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
trait.posterior <- function(model, coefs = NULL, cores = 25, tol = 1e-4) {
  
  n <- model$design_list$n
  I <- model$design_list$I
  q <- model$design_list$q
  n_sigma <- model$design_list$n_sigma
  
  ## automatically get coefs if not provided
  if(is.null(coefs)){
    if(!is.null(model$cv_error)){
      cat("No coefs are provided, automatically cv-optimal model is chosen","\n")
      coefs <- model$coefficients[which.min(model$cv_error),]
    }else{
      cat("No coefs are provided, automatically BIC-optimal model is chosen","\n")
      coefs <- model$coefficients[which.min(model$BIC),]#
    }
  }
  

  ## separate sigma and other coefficients
  coef_short <- head(coefs, length(coefs)-n_sigma)
  sigma <- tail(coefs, n_sigma)
  if(n_sigma==1){
    sigma <- rep(sigma, I)
  }
    
  ## create person index
  person.index <- rep(1:n, each = sum(q))
  
  Y <- matrix(as.numeric(as.matrix(model$Y)), ncol = ncol(model$Y))
  
  ## create linear predictor without theta
  all.lin.preds <- lin_pred(model, coef_short, sigma)  

  estimates <- person.fit(n, I, q, Y, sigma, person.index, 
    all.lin.preds, cores, tol)
  
  names(estimates) <- rownames(model$data)
  
  estimates
}


person.fit <- function(n, I, q, Y, sigma, person.index, 
  all.lin.preds, cores = 1, tol = 1e-3) {
  
  
  limit <- qnorm(0.999)
  
  if (cores > 1) {
    cl <- makeCluster(cores, outfile = "")
    
    clusterExport(cl, varlist = c("Y", "factor_i", "person.index", 
      "all.lin.preds", "sigma", "q", "estim_i", "likeli_i", 
       "tol", "limit"), 
      envir = sys.frame(sys.nframe()))
    
    estimates <- parSapply(cl, seq(n), person_i, person.index = person.index, 
      all.lin.preds = all.lin.preds, sigma = sigma, Y = Y, 
      tol = tol, limit = limit, q = q)
    stopCluster(cl)
  } else {
    estimates <- sapply(seq(n), person_i, person.index = person.index, 
      all.lin.preds = all.lin.preds, sigma = sigma, Y = Y, 
      tol = tol, limit = limit, q = q)
  }
  

  estimates
}