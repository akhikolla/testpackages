#' Bootstrap function for BTLLasso
#' 
#' Performs bootstrap for BTLLasso to get bootstrap intervals. Main
#' input argument is a \code{cv.BTLLasso} object. The bootstrap is (recommended to be) 
#' performed on level of the cross-validation. Therefore, within every bootstrap iteration
#' the complete cross-validation procedure from the \code{cv.BTLLasso} object
#' is performed. A \code{\link[=plot.boot.BTLLasso]{plot}} function can be applied 
#' to the resulting \code{boot.BTLLasso} object to plot bootstrap intervals.
#' 
#' The method can be highly time-consuming, for high numbers of tuning
#' parameters, high numbers of folds in the cross-validation and high number of
#' bootstrap iterations B.  The number of tuning parameters can be reduced by
#' specifying \code{lambda} in the \code{boot.BTLLasso} function. You can control if
#' the range of prespecified tuning parameters was to small by looking at the
#' output values \code{lambda.max.alert} and \code{lambda.min.alert}. They are
#' set \code{TRUE} if the smallest or largest of the specifed lambda values was
#' chosen in at least one bootstrap iteration.
#' 
#' @param model A \code{cv.BTLLasso} object.
#' @param B Number of bootstrap iterations.
#' @param lambda Vector of tuning parameters. If not specified (default),
#' tuning parameters from \code{cv.BTLLasso} object are used. See also details.
#' @param cores Number of cores for (parallelized) computation.
#' @param trace Should the trace of the BTLLasso algorithm be printed?
#' @param trace.cv Should the trace fo the cross-validation be printed? If
#' parallelized, the trace is not working on Windows machines.
#' @param with.cv Should cross-validation be performed separately on every 
#' bootstrap sample? If \code{FALSE}, the tuning parameter is fixed to the value chosen 
#' in the \code{cv.BTLLasso} object. 
#' @return \item{cv.model}{\code{cv.BTLLasso} object} \item{estimatesB}{Matrix
#' containing all B estimates for original parameters. For internal use.}
#' \item{estimatesBrepar}{Matrix containing all B estimates for reparameterized
#' (symmetric side constraints) parameters.} \item{lambdaB}{vector of used
#' tuning parameters} \item{lambda.max.alert}{Was the largest value of lambda chosen
#' in at least one bootstrap iteration?} \item{lambda.min.alert}{Was the
#' smallest value of lambda chosen in at least one bootstrap iteration?} \item{number.na}{Total number
#' of failed bootstrap iterations.}
#' @author Gunther Schauberger\cr \email{gunther.schauberger@@tum.de}
#' @seealso \code{\link{BTLLasso}}, \code{\link{cv.BTLLasso}},
#' \code{\link{plot.boot.BTLLasso}}
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
#' @keywords BTLLasso interval bootstrap
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
boot.BTLLasso <- function(model, B = 500, lambda = NULL, cores = 1, 
                          trace = TRUE, trace.cv = TRUE, with.cv = TRUE){
  
  cv.crit <- model$cv.crit
  
  design <- model$design$design
  response <- model$response
  penalty <- model$penalty
  control <- model$control
  
  q <- model$Y$q
  k <- model$Y$k
  n.design <- nrow(design)/q
  m <- model$Y$m
  folds <- model$folds
  
  if (is.null(lambda)) {
    lambda <- model$lambda
  }

  boot.fun <- function(b) {
    cat("Bootstrap sample:", b, "out of", B, "\n")
    
    sample.b <- sample(x = 1:n.design, size = n.design, replace = TRUE)
    
    id.vec <- c(t(matrix(1:nrow(design), ncol = q, byrow = TRUE)[sample.b, 
      ]))
    
    design.b <- design[id.vec, ]
    response.b <- response[id.vec]
    
    if(with.cv){
      model.b <- try(fit.cv.BTLLasso(response = response.b,
        design = design.b, penalty = penalty, q = q, m = m,
        folds = folds, lambda = lambda, control = control,
        cores = 1, trace = trace, trace.cv = trace.cv, cv.crit = cv.crit))
      
        if (inherits(model.b, "try-error")) {
          coef.b <- NA
          lambda.b <- NA
        } else {
          coef.b <- model.b$coefs[which.min(model.b$criterion),
            ]
          lambda.b <- lambda[which.min(model.b$criterion)]
        }
        }else{
          lambda.b <- model$lambda[which.min(model$criterion)]
          model.b <- try(fit.BTLLasso(response = response.b, 
                                      design = design.b, penalty = penalty, 
                                      lambda = lambda.b, k = k, m = m,
                                      control = control, trace = trace))

            if (inherits(model.b, "try-error")) {
              coef.b <- NA
              lambda.b <- NA
            } else {
              coef.b <- c(model.b$coefs)
            }
        }

    return(list(coef.b = coef.b, lambda.b = lambda.b))
  }
  
  
  if (cores > 1) {
    cl <- makeCluster(cores, outfile = "")
    clusterSetRNGStream(cl, NULL)
    clusterExport(cl, varlist = c("design", "response", "penalty", 
      "q", "m", "control", "folds", "lambda", "cores", 
      "B", "n.design", "trace", "trace.cv", "cv.crit"), 
      envir = sys.frame(sys.nframe()))
    
    outputB <- parLapply(cl, seq(B), boot.fun)
    stopCluster(cl)
  } else {
    outputB <- lapply(seq(B), boot.fun)
  }
  
  
  estimatesB <- matrix(0, ncol = ncol(design), nrow = B)
  lambdaB <- c()
  for (b in 1:B) {
    if (any(is.na(outputB[[b]]$coef.b))) {
      estimatesB[b, ] <- rep(NA, ncol(estimatesB))
      lambdaB[b] <- NA
    } else {
      estimatesB[b, ] <- outputB[[b]]$coef.b
      lambdaB[b] <- outputB[[b]]$lambda.b
    }
  }
  
  number.na <- sum(rowSums(is.na(estimatesB))>0)
  if(number.na>0){
    warning(number.na, " out of ", B, " bootstrap samples did not converge!")
  }

  estimatesBrepar <- round(expand.coefs(estimatesB, model$design, 
    model$Y, symmetric = TRUE, model$control$name.order), model$control$precision)
  
  estimatesB <- round(estimatesB, model$control$precision)
  
  # conf.ints <- apply(estimatesB, 2, quantile, probs = quantiles, 
  #   type = 1, na.rm = TRUE)
  # conf.ints.repar <- apply(estimatesBrepar, 2, quantile, probs = quantiles, 
  #   type = 1, na.rm = TRUE)
  
  lambda.min.alert <- any(lambdaB == min(lambda))
  lambda.max.alert <- any(lambdaB == max(lambda))
  
  returns <- list(cv.model = model, estimatesB = estimatesB, 
    estimatesBrepar = estimatesBrepar, lambdaB = lambdaB, 
    lambda.max.alert = lambda.max.alert, lambda.min.alert = lambda.min.alert,
    number.na = number.na)
  
  class(returns) <- "boot.BTLLasso"
  
  return(returns)
  
}
