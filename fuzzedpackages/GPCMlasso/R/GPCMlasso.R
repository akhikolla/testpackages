#' GPCMlasso
#' 
#' Performs GPCMlasso, a method to identify differential item functioning (DIF) in Generalized Partial Credit Models. 
#' A joint parametric model is set up based on an IRT model chosen by the user. Several variables can be considered simultaneously.
#' For each pair between variable and item, a parametric DIF effect is introduced which
#' indicates DIF if the respective parameter is selected (estimated to be unequal zero). 
#' Parameter selection is done using a lasso-type penalization term.
#' 
#' @param formula Formula to indicate which items are considered and which covariates 
#' should be used to find DIF. Items are considered to be the response and are 
#' concatenated by \code{cbind()}. If the RHS of the formula is \code{~0}, simply the 
#' model specified in \code{model} is calulated.
#' @param data Data frame containing the ordinal item response data 
#' (as ordered factors) and all covariates.
#' @param DSF Should Differential Step Functioning (DSF) be considered? If \code{DSF = TRUE}, one parameter per step between two response categories is introduced. For binary items, DSF and DIF conincide. 
#' @param model Specify the underlying basic model. Currently, you can choose 
#' between the partial credit model and the rating scale model and
#' the respective generalized versions of both models called 'PCM', 'RSM', 
#' 'GPCM' and 'GRSM'. Generalized models allow for different discrimination 
#' parameters between items. 
#' @param control Control argument to specify further arguments for the algorithm
#' and numerical optimization, specified by \code{\link{ctrl_GPCMlasso}}.
#' @param cv Should cross-validation be performed? Cross-validation can be used as an alternative to BIC to select the optimal tuning parameter.
#' @param main.effects Should also main effects of the variables be included in the model? Default is \code{TRUE}. Here, positive parameter estimates correspond 
#' to an increase of the respective trait if the variable increases.
#' @return 
#' \item{coefficients}{Matrix containing all parameters for the GPCMlasso model, one row
#' per tuning parameter lambda. Due to the penalty the parameters are scaled and, therefore,
#' are comparable with respect to their size.} 
#' \item{logLik}{Vector of log-likelihoods, one value per tuning parameter lambda.}
#' \item{call}{The function call of \code{GPCMlasso}}
#' \item{cv_error}{Vector of cv_errors, one per tuning parameter. Only relevant if \code{cv = TRUE}.}
#' \item{model}{Basic IRT model chosen by user.}
#' \item{data}{Data from call.}
#' \item{control}{Control list.}
#' \item{DSF}{DSF from call.}
#' \item{formula}{Formula from call.}
#' \item{item.names}{Item names.}
#' \item{Y}{Matrix containing item responses.}
#' \item{design_list}{List containing several helpful objects for internal use.}
#' \item{AIC}{Vector of AIC values, one per tuning parameter.}
#' \item{BIC}{Vector of BIC values, one per tuning parameter.}
#' \item{cAIC}{Vector of corrected AIC values, one per tuning parameter.}
#' \item{df}{Vector of degrees of freedom, one per tuning parameter.}
#' \item{coef.rescal}{Matrix containing all rescaled parameters for the GPCMlasso model, one row
#' per tuning parameter lambda. In contrast to \code{coefficients}, all parameters are rescaled back
#' to their original scales.}
#' @author Gunther Schauberger\cr \email{gunther.schauberger@@tum.de}
#' @references Schauberger, Gunther and Mair, Patrick (2019): A Regularization Approach for the Detection of Differential 
#' Item Functioning in Generalized Partial Credit Models, \emph{Behavior Research Methods}, \url{https://link.springer.com/article/10.3758/s13428-019-01224-2}
#' @seealso \code{\link{GPCMlasso-package}}, \code{\link{ctrl_GPCMlasso}}, \code{\link{print.GPCMlasso}}, 
#' \code{\link{plot.GPCMlasso}}, \code{\link{predict.GPCMlasso}}, \code{\link{trait.posterior}}
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
GPCMlasso <- function(formula, data, DSF = FALSE, 
                      model = c("PCM","RSM","GPCM","GRSM","RM","2PL"),
                      control = ctrl_GPCMlasso(), cv = FALSE, main.effects = TRUE){

  ### save old contrasts to set them back to default later
  # old.contrasts <- options()$contrasts
  # options(contrasts = control$contrasts)
  # on.exit(options(contrasts = old.contrasts))
  
  ## check if data is data.frame
  if(!is.data.frame(data))
    stop("data has to be a data.frame!")
  
  ## apply model response and get Y
  Y <- model.response(model.frame(formula, data = data))
  ## factorize, also for binary Y is now 1,2,...
  for(i in 1:ncol(Y)){
    Y[,i] <- as.factor(Y[,i])
  }
  ## get item names
  item.names <- colnames(Y)
  
  ## get model, define respective dummies
  model <- match.arg(model)
  ## RSM makes rsm paramterization
  ## GPCM uses different instead of equal discrim parameters
  RSM <- GPCM <- FALSE
  if(model %in% c("GRSM","RSM")){
    RSM <- TRUE
  }
  if(model %in% c("GRSM","GPCM","2PL")){
    GPCM <- TRUE
  }

  ## for binary, DSF or cross-validation with deviance does not make sense
  k <- apply(Y, 2, function(x){length(levels(as.factor(x)))})
  
  if(all(k==2) | model %in% c("RM","2PL")){
    DSF <- FALSE
    control$cv.crit <- "deviance"
  }
  
  ## check if for RSM == TRUE all numbers of categories are equal
  if(length(unique(k))>1 & RSM){
    stop("Models RSM and GRSM cannot be used for unequal numbers of response categories!")
  }
  

  ## get design matrices and all further parameters  
  design_list <- design_GPCMlasso(formula = formula, data = data, 
                                  Y = Y, RSM = RSM, GPCM = GPCM, 
                                  DSF = DSF, all.dummies = control$all.dummies,
                                  main.effects = main.effects)

  
  ## if no covariates are used, we don't need lambda or adaptive
  if(design_list$m==0){
    control$lambda <- 0
    control$adaptive <- FALSE
    cv <- FALSE
  }
  if(length(control$lambda)==1){
    cv <- FALSE
  }
  
  ## which loglik and which score function is used,
  loglik_fun <- loglikPCMlasso
  score_fun <- scorePCMlasso
  log_score_fun <- loglikscorePCMlasso2
  if(all(k==2)){
    loglik_fun <- loglikDIFlasso
    score_fun <- scoreDIFlasso
    log_score_fun <- loglikscoreDIFlasso2
  }
  
  
  ## get fitted parameters and logliks
  fit <- fit_GPCMlasso(model = model,loglik_fun = loglik_fun, score_fun = score_fun, log_score_fun,
                       design_list = design_list, 
                       control = control,  start = NULL, scale_fac = 1, main.effects = main.effects)
  
  ## update lambda if lambda was not specified in advance
  if(is.null(control$lambda)){
    control$lambda <- fit$lambda
  }
  
  ## extract coefficients
  coefficients <-  fit$coefficients
  
  coef.rescal <- fit$coef.rescal
  
  ##extract vector of logliks and degrees of freedom
  logLik <- fit$logLik
  df <- fit$df
  
  BIC <- -2*logLik + log(design_list$n)*df
  AIC <- -2*logLik + 2*df
  cAIC <- -2*logLik + 2*df +2*df*(df-1)/(design_list$n-df-1)
  

  ## apply cross validation if wanted
  if(cv){
    cv_error <- fit_cv_GPCMlasso(model, design_list, 
                                  control, score_fun,
                                  loglik_fun, log_score_fun, Y)
  }else{
    cv_error <- NULL
  }

  ## create return list
  ret.list <- list(coefficients = coefficients, logLik = logLik, 
                   cv_error = cv_error, call = match.call(), 
                   model = model, data = data, control = control,
                   DSF = DSF, formula = formula,
                   item.names = item.names, Y = Y, design_list = design_list,
                   AIC = AIC, BIC = BIC, cAIC = cAIC, df = df, 
                   coef.rescal = coef.rescal, main.effects = main.effects)
  
  class(ret.list) <- "GPCMlasso"
  
  return(ret.list)
}
