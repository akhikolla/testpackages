#' Control function for BTLLasso
#' 
#' Control parameters for different penalty terms and for tuning the fitting algorithm.
#' 
#' @param l.lambda Number of tuning parameters. Applies only if \code{lambda = NULL} in the 
#' main function. 
#' @param log.lambda Should the grid of tuning parameters be created on a logarithmic scale 
#' rather than equidistant. Applies only if \code{lambda = NULL} in the main function. 
#' @param lambda.min Minimal value for tuning parameter. Applies only if \code{lambda = NULL} in the 
#' main function. 
#' @param adaptive Should adaptive lasso be used? Default is TRUE.
#' @param scale Should the covariates be scaled so that they are on comparable scales? Default is TRUE.
#' Variables will be properly scaled AND centered. Please note that results will refer to scaled covariates.
#' If \code{adaptive = TRUE} scaling is not necessary to keep penalties comparable.
#' @param norm Specifies the norm used in the penalty term. Currently, only
#' 'L1' and 'L2' are possible. Default is to 'L1', only 'L1' allows for
#' clustering and variable selection.
#' @param epsilon Threshold value for convergence of the algorithm.
#' @param lambda2 Tuning parameter for ridge penalty on all coefficients.
#' Should be small, only used to stabilize results.
#' @param c Internal parameter for the quadratic approximation of the L1
#' penalty. Should be sufficiently small. For details see
#' \code{\link[gvcm.cat]{cat_control}}.
#' @param precision Precision for final parameter estimates, specifies number of decimals.
#' @param weight.penalties Should the penalties across the different model components 
#' (i.e. intercepts, order effects, X, Z1, Z2) be weighted according to the number of
#' penalties included? Default is \code{TRUE} to minimize the risk of selection bias
#' across different model components.  
#' @param include.intercepts Should intercepts be included in the model?
#' @param order.effect Should a global order effect (corresponding to home effect
#' in sports applications) be included in the model?
#' @param object.order.effect Should object-specific order effects (corresponding to home effects
#' in sports applications) be included in the model?
#' @param order.center Should (in case of object-specific order effects) the order effects be centered
#' in the design matrix? Centering is equivalent to the coding scheme of effect coding instead of 
#' dummy coding.
#' @param name.order How should the order effect(s) be called in plots or prints.
#' @param penalize.intercepts Should intercepts be penalized? If \code{TRUE},
#' all pairwise differences between intercepts are penalized.
#' @param penalize.X Should effects from X matrix be penalized? If \code{TRUE},
#' all pairwise differences corresponding to one covariate are penalized. Can also be used with
#' a character vector as input. Then, the character vector contains the names of the variables
#' from X whose parameters should be penalized.
#' @param penalize.Z2 Should absolute values of effects from Z2 matrix be
#' penalized? Can also be used with
#' a character vector as input. Then, the character vector contains the names of the variables
#' from Z2 whose parameters should be penalized.
#' @param penalize.Z1.absolute Should absolute values of effects from Z1 matrix
#' be penalized? Can also be used with
#' a character vector as input. Then, the character vector contains the names of the variables
#' from Z1 whose parameters should be penalized.
#' @param penalize.Z1.diffs Should differences of effects from Z1 matrix be
#' penalized? If \code{TRUE}, all pairwise differences corresponding to one
#' covariate are penalized. Can also be used with
#' a character vector as input. Then, the character vector contains the names of the variables
#' from Z1 whose parameters should be penalized.
#' @param penalize.order.effect.absolute Should absolute values of order effect(s) be penalized?
#' Only relevant if either \code{object.order.effect = TRUE} or \code{order.effect = TRUE}.
#' @param penalize.order.effect.diffs Should differences of order effects be
#' penalized? If \code{TRUE}, all pairwise differences are penalized. Only relevant if 
#' \code{object.order.effect = TRUE}
#' @author Gunther Schauberger\cr \email{gunther.schauberger@@tum.de}
#' @seealso \code{\link{BTLLasso}}, \code{\link{cv.BTLLasso}}
#' @references Schauberger, Gunther and Tutz, Gerhard (2019): BTLLasso - A Common Framework and Software 
#' Package for the Inclusion  and Selection of Covariates in Bradley-Terry Models, \emph{Journal of 
#' Statistical Software}, 88(9), 1-29, \url{https://doi.org/10.18637/jss.v088.i09}
#' 
#' Schauberger, Gunther and Tutz, Gerhard (2017): Subject-specific modelling 
#' of paired comparison data: A lasso-type penalty approach, \emph{Statistical Modelling},
#' 17(3), 223 - 243
#' 
#' Schauberger, Gunther, Groll Andreas and Tutz, Gerhard (2018): 
#' Analysis of the importance of on-field covariates in the German Bundesliga, 
#' \emph{Journal of Applied Statistics}, 45(9), 1561 - 1578
#' @keywords BTLLasso control
#' @examples
#' 
#' \dontrun{
#' op <- par(no.readonly = TRUE)
#' 
#' ##############################
#' ##### Example with simulated data set containing X, Z1 and Z2
#' ##############################
#' data(SimData)
#' 
#' ## Specify control argument
#' ## -> allow for object-specific order effects and penalize intercepts
#' ctrl <- ctrl.BTLLasso(penalize.intercepts = TRUE, object.order.effect = TRUE,
#'                       penalize.order.effect.diffs = TRUE)
#' 
#' ## Simple BTLLasso model for tuning parameters lambda
#' m.sim <- BTLLasso(Y = SimData$Y, X = SimData$X, Z1 = SimData$Z1,
#'                   Z2 = SimData$Z2, control = ctrl)
#' m.sim
#' 
#' par(xpd = TRUE)
#' plot(m.sim)
#' 
#' 
#' ## Cross-validate BTLLasso model for tuning parameters lambda
#' set.seed(1860)
#' m.sim.cv <- cv.BTLLasso(Y = SimData$Y, X = SimData$X, Z1 = SimData$Z1,
#'                         Z2 = SimData$Z2, control = ctrl)
#' m.sim.cv
#' coef(m.sim.cv)
#' logLik(m.sim.cv)
#' 
#' head(predict(m.sim.cv, type="response"))
#' head(predict(m.sim.cv, type="trait"))
#' 
#' plot(m.sim.cv, plots_per_page = 4)
#' 
#' 
#' ## Example for bootstrap intervals for illustration only
#' ## Don't calculate bootstrap intervals with B = 20!!!!
#' set.seed(1860)
#' m.sim.boot <- boot.BTLLasso(m.sim.cv, B = 20, cores = 20)
#' m.sim.boot
#' plot(m.sim.boot, plots_per_page = 4)
#' 
#' 
#' ##############################
#' ##### Example with small version from GLES data set
#' ##############################
#' data(GLESsmall)
#' 
#' ## extract data and center covariates for better interpretability
#' Y <- GLESsmall$Y
#' X <- scale(GLESsmall$X, scale = FALSE)
#' Z1 <- scale(GLESsmall$Z1, scale = FALSE)
#' 
#' ## vector of subtitles, containing the coding of the X covariates
#' subs.X <- c('', 'female (1); male (0)')
#' 
#' ## Cross-validate BTLLasso model
#' m.gles.cv <- cv.BTLLasso(Y = Y, X = X, Z1 = Z1)
#' m.gles.cv
#' 
#' coef(m.gles.cv)
#' logLik(m.gles.cv)
#' 
#' head(predict(m.gles.cv, type="response"))
#' head(predict(m.gles.cv, type="trait"))
#' 
#' par(xpd = TRUE, mar = c(5,4,4,6))
#' plot(m.gles.cv, subs.X = subs.X, plots_per_page = 4, which = 2:5)
#' paths(m.gles.cv, y.axis = 'L2')
#' 
#' 
#' ##############################
#' ##### Example with Bundesliga data set
#' ##############################
#' data(Buli1516)
#' 
#' Y <- Buli1516$Y5
#' 
#' Z1 <- scale(Buli1516$Z1, scale = FALSE)
#' 
#' ctrl.buli <- ctrl.BTLLasso(object.order.effect = TRUE, 
#'                            name.order = "Home", 
#'                            penalize.order.effect.diffs = TRUE, 
#'                            penalize.order.effect.absolute = FALSE,
#'                            order.center = TRUE, lambda2 = 1e-2)
#' 
#' set.seed(1860)
#' m.buli <- cv.BTLLasso(Y = Y, Z1 = Z1, control = ctrl.buli)
#' m.buli
#' 
#' par(xpd = TRUE, mar = c(5,4,4,6))
#' plot(m.buli)
#' 
#' 
#' ##############################
#' ##### Example with Topmodel data set
#' ##############################
#' data("Topmodel2007", package = "psychotree")
#' 
#' Y.models <- response.BTLLasso(Topmodel2007$preference)
#' X.models <- scale(model.matrix(preference~., data = Topmodel2007)[,-1])
#' rownames(X.models) <- paste0("Subject",1:nrow(X.models))
#' colnames(X.models) <- c("Gender","Age","KnowShow","WatchShow","WatchFinal")
#' 
#' set.seed(5)
#' m.models <- cv.BTLLasso(Y = Y.models, X = X.models)
#' plot(m.models, plots_per_page = 6)
#' 
#' par(op)
#' }
ctrl.BTLLasso <- function(l.lambda = 30, log.lambda = TRUE, lambda.min = 0.05,
                          adaptive = TRUE, scale = TRUE, norm = c("L1", 
  "L2"), epsilon = 1e-04, lambda2 = 1e-04, c = 1e-09, precision = 3, 
  weight.penalties = TRUE, include.intercepts = TRUE, order.effect = FALSE, 
  object.order.effect = FALSE, order.center = FALSE, name.order = "Order", 
  penalize.intercepts = FALSE, penalize.X = TRUE, penalize.Z2 = FALSE, 
  penalize.Z1.absolute = TRUE, penalize.Z1.diffs = TRUE, penalize.order.effect.absolute = TRUE, 
  penalize.order.effect.diffs = FALSE) {
  norm <- match.arg(norm)
  
  RET <- list(l.lambda = l.lambda, log.lambda = log.lambda, lambda.min = lambda.min, 
              adaptive = adaptive, scale = scale, norm = norm, 
    epsilon = epsilon, lambda2 = lambda2, c = c, penalize.X = penalize.X, 
    penalize.Z1.diffs = penalize.Z1.diffs, penalize.Z2 = penalize.Z2, 
    penalize.Z1.absolute = penalize.Z1.absolute, penalize.intercepts = penalize.intercepts, 
    include.intercepts = include.intercepts, order.effect = order.effect, 
    object.order.effect = object.order.effect, order.center = order.center, 
    penalize.order.effect.diffs = penalize.order.effect.diffs, 
    penalize.order.effect.absolute = penalize.order.effect.absolute, 
    name.order = name.order, precision = precision, weight.penalties = weight.penalties)
  RET
}




#' Function to perform BTLLasso
#' 
#' Performs BTLLasso, a method to model heterogeneity in paired comparison
#' data. Different types of covariates are allowd to have an influence on the
#' attractivity/strength of the objects. Covariates can be subject-specific, 
#' object-specific or subject-object-specific. L1 penalties are used to reduce the 
#' complexiy of the model by enforcing clusters of equal effects or by elimination of irrelevant
#' covariates.  
#' 
#' 
#' @param Y A \code{response.BTLLasso} object created by
#' \code{\link{response.BTLLasso}}.
#' @param X Matrix containing all \bold{subject-specific covariates} that are
#' to be included with \bold{object-specific effects}. One row represents one
#' subject, one column represents one covariate. X has to be standardized.
#' @param Z1 Matrix containing all \bold{object-subject-specific covariates}
#' that are to be included with \bold{object-specific effects}. One row
#' represents one subject, one column represents one combination between
#' covariate and object. Column names have to follow the scheme
#' 'firstvar.object1',...,'firstvar.objectm',...,'lastvar.objectm'. The object
#' names 'object1',...,'objectm' have to be identical to the object names used
#' in the \code{response.BTLLasso} object \code{Y}. The variable names and the
#' object names have to be separated by '.'.  The rownames of the matrix',
#' Z.name, 'have to be equal to the subjects specified in the response object.
#' Z1 has to be standardized.
#' @param Z2 Matrix containing all \bold{object-subject-specific covariates or
#' object-specific covariates} that are to be included with \bold{global
#' effects}. One row represents one subject, one column represents one
#' combination between covariate and object. Column names have to follow the
#' scheme 'firstvar.object1',...,'firstvar.objectm',...,'lastvar.objectm'. The
#' object names 'object1',...,'objectm' have to be identical to the object
#' names used in the \code{response.BTLLasso} object \code{Y}. The variable
#' names and the object names have to be separated by '.'.  The rownames of the
#' matrix', Z.name, 'have to be equal to the subjects specified in the response
#' object. Z2 has to be standardized.
#' @param lambda Vector of tuning parameters. If \code{NULL}, automatically a grid
#' of tuning parameters is created. 
#' @param control Function for control arguments, mostly for internal use. See
#' also \code{\link{ctrl.BTLLasso}}.
#' @param trace Should the trace of the BTLLasso algorithm be printed?
#' @return 
#' \item{coefs}{Matrix containing all (original) coefficients, one row
#' per tuning parameter, one column per coefficient.} 
#' \item{coefs.repar}{Matrix
#' containing all reparameterized (for symmetric side constraint) coefficients,
#' one row per tuning parameter, one column per coefficient.}
#' \item{logLik}{Vector of log-likelihoods, one value per tuning parameter.}
#' \item{design}{List containing design matrix and several additional information like, 
#' e.g., number and names of covariates.} 
#' \item{Y}{Response object.} 
#' \item{penalty}{List containing all penalty matrices and some further information on penalties.} 
#' \item{response}{Vector containing 0-1 coded
#' response.} 
#' \item{X}{X matrix containing subject-specific covariates.} 
#' \item{Z1}{Z1 matrix containing subject-object-specific covariates.} 
#' \item{Z2}{Z2 matrix containing (subject)-object-specific covariates.} 
#' \item{lambda}{Vector of tuning parameters.} 
#' \item{control}{Control argument, specified by \code{\link{ctrl.BTLLasso}}.}
#' \item{df}{Vector containing degrees of freedom for all models along the grid 
#' of tuning parameters.}
#' @author Gunther Schauberger\cr \email{gunther.schauberger@@tum.de}
#' @seealso \code{\link{cv.BTLLasso}}, \code{\link{boot.BTLLasso}}, \code{\link{ctrl.BTLLasso}},
#' \code{\link{plot.BTLLasso}}, \code{\link{paths}}, \code{\link{print.BTLLasso}}, 
#' \code{\link{predict.BTLLasso}}, \code{\link{coef}}
#' @references Schauberger, Gunther and Tutz, Gerhard (2019): BTLLasso - A Common Framework and Software 
#' Package for the Inclusion  and Selection of Covariates in Bradley-Terry Models, \emph{Journal of 
#' Statistical Software}, to appear
#' 
#' Schauberger, Gunther and Tutz, Gerhard (2017): Subject-specific modelling 
#' of paired comparison data: A lasso-type penalty approach, \emph{Statistical Modelling},
#' 17(3), 223 - 243
#' 
#' Schauberger, Gunther, Groll Andreas and Tutz, Gerhard (2018): 
#' Analysis of the importance of on-field covariates in the German Bundesliga, 
#' \emph{Journal of Applied Statistics}, 45(9), 1561 - 1578
#' @keywords BTLLasso
#' @examples
#' 
#' \dontrun{
#' op <- par(no.readonly = TRUE)
#' 
#' ##############################
#' ##### Example with simulated data set containing X, Z1 and Z2
#' ##############################
#' data(SimData)
#' 
#' ## Specify control argument
#' ## -> allow for object-specific order effects and penalize intercepts
#' ctrl <- ctrl.BTLLasso(penalize.intercepts = TRUE, object.order.effect = TRUE,
#'                       penalize.order.effect.diffs = TRUE)
#' 
#' ## Simple BTLLasso model for tuning parameters lambda
#' m.sim <- BTLLasso(Y = SimData$Y, X = SimData$X, Z1 = SimData$Z1,
#'                   Z2 = SimData$Z2, control = ctrl)
#' m.sim
#' 
#' par(xpd = TRUE)
#' plot(m.sim)
#' 
#' 
#' ## Cross-validate BTLLasso model for tuning parameters lambda
#' set.seed(1860)
#' m.sim.cv <- cv.BTLLasso(Y = SimData$Y, X = SimData$X, Z1 = SimData$Z1,
#'                         Z2 = SimData$Z2, control = ctrl)
#' m.sim.cv
#' coef(m.sim.cv)
#' logLik(m.sim.cv)
#' 
#' head(predict(m.sim.cv, type="response"))
#' head(predict(m.sim.cv, type="trait"))
#' 
#' plot(m.sim.cv, plots_per_page = 4)
#' 
#' 
#' ## Example for bootstrap intervals for illustration only
#' ## Don't calculate bootstrap intervals with B = 20!!!!
#' set.seed(1860)
#' m.sim.boot <- boot.BTLLasso(m.sim.cv, B = 20, cores = 20)
#' m.sim.boot
#' plot(m.sim.boot, plots_per_page = 4)
#' 
#' 
#' ##############################
#' ##### Example with small version from GLES data set
#' ##############################
#' data(GLESsmall)
#' 
#' ## extract data and center covariates for better interpretability
#' Y <- GLESsmall$Y
#' X <- scale(GLESsmall$X, scale = FALSE)
#' Z1 <- scale(GLESsmall$Z1, scale = FALSE)
#' 
#' ## vector of subtitles, containing the coding of the X covariates
#' subs.X <- c('', 'female (1); male (0)')
#' 
#' ## Cross-validate BTLLasso model
#' m.gles.cv <- cv.BTLLasso(Y = Y, X = X, Z1 = Z1)
#' m.gles.cv
#' 
#' coef(m.gles.cv)
#' logLik(m.gles.cv)
#' 
#' head(predict(m.gles.cv, type="response"))
#' head(predict(m.gles.cv, type="trait"))
#' 
#' par(xpd = TRUE, mar = c(5,4,4,6))
#' plot(m.gles.cv, subs.X = subs.X, plots_per_page = 4, which = 2:5)
#' paths(m.gles.cv, y.axis = 'L2')
#' 
#' 
#' ##############################
#' ##### Example with Bundesliga data set
#' ##############################
#' data(Buli1516)
#' 
#' Y <- Buli1516$Y5
#' 
#' Z1 <- scale(Buli1516$Z1, scale = FALSE)
#' 
#' ctrl.buli <- ctrl.BTLLasso(object.order.effect = TRUE, 
#'                            name.order = "Home", 
#'                            penalize.order.effect.diffs = TRUE, 
#'                            penalize.order.effect.absolute = FALSE,
#'                            order.center = TRUE, lambda2 = 1e-2)
#' 
#' set.seed(1860)
#' m.buli <- cv.BTLLasso(Y = Y, Z1 = Z1, control = ctrl.buli)
#' m.buli
#' 
#' par(xpd = TRUE, mar = c(5,4,4,6))
#' plot(m.buli)
#' 
#' 
#' ##############################
#' ##### Example with Topmodel data set
#' ##############################
#' data("Topmodel2007", package = "psychotree")
#' 
#' Y.models <- response.BTLLasso(Topmodel2007$preference)
#' X.models <- scale(model.matrix(preference~., data = Topmodel2007)[,-1])
#' rownames(X.models) <- paste0("Subject",1:nrow(X.models))
#' colnames(X.models) <- c("Gender","Age","KnowShow","WatchShow","WatchFinal")
#' 
#' set.seed(5)
#' m.models <- cv.BTLLasso(Y = Y.models, X = X.models)
#' plot(m.models, plots_per_page = 6)
#' 
#' par(op)
#' }
BTLLasso <- function(Y, X = NULL, Z1 = NULL, Z2 = NULL, lambda = NULL, 
  control = ctrl.BTLLasso(), trace = TRUE) {
  

  ## create design matrix
  get.design <- design.BTLLasso(Y = Y, X = X, Z1 = Z1, Z2 = Z2, 
    control = control)
  
  ## exclude missing values
  na.response <- is.na(Y$response)
  na.design <- colSums(matrix(is.na(rowSums(get.design$design)), 
    nrow = Y$q)) != 0
  na.total <- na.response | na.design
  Y$response <- Y$response[!na.total]
  Y$first.object <- Y$first.object[!na.total]
  Y$second.object <- Y$second.object[!na.total]
  Y$subject <- Y$subject[!na.total]
  Y$subject.names <- levels(as.factor(Y$subject))
  Y$n <- length(Y$subject.names)
  
  get.design$design <- get.design$design[!rep(na.total, each = Y$q), 
    ]
  
  get.design$design.repar <- get.design$design.repar[!rep(na.total, each = Y$q), 
                                         ]

  ## create response vector
  if(identical(levels(Y$response),c("0","1"))){
      response <- as.numeric(Y$response) -1
  } else {
    response <- cumul.response(Y)
  }

  ## create penalty matrix
  get.penalties <- penalties.BTLLasso(Y = Y, X = X, Z1 = Z1, 
    Z2 = Z2, control = control, get.design = get.design)
  
  ## create sequence of tuning parameters if not pre-specified
  if(is.null(lambda)){
    lambda <- find.lambda(response = response, 
                          design = get.design$design, penalties = get.penalties, 
                           k = Y$k, m = Y$m, control = control, trace = trace)
  }
  
  ## fit BTLLasso model, with coefficients and degrees of
  ## freedom
  fit <- fit.BTLLasso(response, get.design$design, get.penalties, 
    lambda, Y$k, Y$m, control, trace)
  coefs <- fit$coefs
  
  
  ## reparameterize coefficients, from reference object to
  ## symmetric side constraint
  coefs.repar <- round(expand.coefs(coefs, get.design, Y, name.order = control$name.order), 
                       control$precision)
  
  
  ## calculate log likelihood
  logLik <- c()
  for (j in 1:nrow(coefs)) {
    logLik[j] <- loglik(coefs[j, ], Y$response, get.design$design, 
      Y$k)
  }
  
  coefs <- round(coefs, control$precision)

  df <- df.BTLLasso(coefs.repar, get.design, Y$m)
  
  ## return stuff
  ret.list <- list(coefs = coefs, coefs.repar = coefs.repar, 
    logLik = logLik, design = get.design, Y = Y, penalty = get.penalties, 
    response = response, X = X, Z1 = Z1, Z2 = Z2, lambda = lambda, 
    control = control, df = df)
  
  class(ret.list) <- "BTLLasso"
  
  return(ret.list)
  
}
