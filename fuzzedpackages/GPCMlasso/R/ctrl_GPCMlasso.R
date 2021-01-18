#' Control function for GPCMlasso
#' 
#' Control parameters for penalty terms and for tuning the fitting algorithm.
#' 
#' @param log.lambda Should the grid of tuning parameters be created on a log scale?
#' @param lambda Optional argument to specify a vector of tuning parameters. If \code{lambda = NULL}, a vector of length \code{l.lambda} is created automatically. 
#' @param l.lambda Specifies the length of the grid of tuning parameters.
#' @param lambda.min Minimal value used if the grid of tuning parameters is created automatically. 
#' @param adaptive Should adaptive lasso be used? Default is \code{TRUE}.
#' @param weight.penalties Should penalties be weightes accoreding to the number of penalty term and the number of parameters
#' corresponding to one pair between item and covariate. Only relevant if both \code{DSF = TRUE} and the number of 
#' response categories differs across items (because only then these values can differ).  
#' @param ada.lambda Size of tuning parameter for Ridge-regularized estimation of parameters used for adaptive weights.
#' @param ada.power By default, 1st power of absolute values of Ridge-regularized estimates are used. Could be changed to squared values by \code{ada-power = 2}.
#' @param Q Number of nodes to be used in Gauss-Hermite quadrature.
#' @param lambda2 Tuning parameter for ridge penalty on all coefficients except 
#' sigma/slope parameters.
#' Should be small, only used to stabilize results.
#' @param cvalue Internal parameter for the quadratic approximation of the L1
#' penalty. Should be sufficiently small. For details see
#' \code{\link[gvcm.cat]{cat_control}}.
#' @param trace Should the trace of the progress (current tuning parameter) be printed?
#' @param folds Number of folds for cross-validation. Only relevant if \code{cv = TRUE} in \code{\link{GPCMlasso}}.
#' @param cores Number of cores to be used parallel when fitting the model.
#' @param null_thresh Threshold which is used to distinguih between values equal and unequal to zero.
#' @param gradtol Parameter to tune optimization accuracy, for details see \code{\link{nlm}}.
#' @param steptol Parameter to tune optimization accuracy, for details see \code{\link{nlm}}.
#' @param iterlim Parameter to tune optimization accuracy, for details see \code{\link{nlm}}.
#' @param precision Number of decimal places used to round coefficient estimates. 
#' @param all.dummies Should (in case of factors with more than 2 categories) the dummy variables for all categories be included in the design matrix? If \code{all.dummies = TRUE}, the dependence on the reference category is eliminated for multi-categorical covariates.
#' @author Gunther Schauberger\cr \email{gunther.schauberger@@tum.de}
#' @references Schauberger, Gunther and Mair, Patrick (2019): A Regularization Approach for the Detection of Differential 
#' Item Functioning in Generalized Partial Credit Models, \emph{Behavior Research Methods}, \url{https://link.springer.com/article/10.3758/s13428-019-01224-2}
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
ctrl_GPCMlasso <- function(log.lambda = TRUE, lambda = NULL, l.lambda = 50,
                           lambda.min = 0.1,
                           adaptive = TRUE, weight.penalties = TRUE, ada.lambda = 1e-4, ada.power = 1,
                           Q = 15, lambda2 = 1e-4, cvalue = 1e-5, 
                           trace = TRUE, folds = 10, cores = 25, 
                           null_thresh = 0.01, gradtol = 1e-6,
                           steptol = 1e-6, iterlim = 500, 
                            precision = 3,
                           all.dummies = FALSE){
  
  RET <- list(log.lambda = log.lambda, lambda = lambda, l.lambda = l.lambda, 
              cores = cores, adaptive = adaptive, weight.penalties = weight.penalties, ada.lambda = ada.lambda, Q = Q, 
              lambda2 = lambda2, cvalue = cvalue, trace = trace, ada.power = ada.power,
              folds = folds, null_thresh = null_thresh,
              gradtol = gradtol, steptol = steptol, iterlim = iterlim, 
              lambda.min = lambda.min, precision = precision, 
              all.dummies = all.dummies)
  RET
}
