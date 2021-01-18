#' Print function for cv.BTLLasso objects
#' 
#' Prints the most important output of \code{cv.BTLLasso} objects.
#' 
#' @method print cv.BTLLasso
#' @param x \code{cv.BTLLasso} object
#' @param rescale Should the parameter estimates be rescaled for plotting? Only 
#' applies if \code{scale = TRUE} was specified in \code{BTLLasso} or \code{cv.BTLLasso}.
#' @param \dots possible further arguments for print command
#' @author Gunther Schauberger\cr \email{gunther.schauberger@@tum.de}
#' @seealso \code{\link{cv.BTLLasso}}
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
print.cv.BTLLasso <- function(x, rescale = FALSE, ...) {
  
  m <- x$Y$m
  n <- x$Y$n
  k <- x$Y$q + 1
  n.theta <- x$design$n.theta
  n.intercepts <- x$design$n.intercepts
  if (n.intercepts != 0) {
    n.intercepts <- n.intercepts + 1
  }
  n.order <- x$design$n.order
  p.X <- x$design$p.X
  p.Z1 <- x$design$p.Z1
  p.Z2 <- x$design$p.Z2
  lambda <- x$lambda
  
  vars.X <- x$design$vars.X
  vars.Z1 <- x$design$vars.Z1
  vars.Z2 <- x$design$vars.Z2
  
  labels <- x$Y$object.names
  
  cv.crit <- x$cv.crit
  
  
  cat("Output of BTLLasso estimation:", "\n")
  
  cat("---", "\n")
  
  cat("Setting:")
  cat("\n", n, "subjects")
  cat("\n", m, "objects")
  cat("\n", k, "response categories")
  cat("\n", p.X, "subject-specific covariate(s)")
  cat("\n", p.Z1, "subject-object-specific covariate(s) with object-specific effects")
  cat("\n", p.Z2, "(subject-)object-specific covariate(s) with global effects")
  if (n.order == m) {
    cat("\n", n.order, "subject-specific order effects")
  }
  if (n.order == 1) {
    cat("\n", "Global order effect")
  }
  if (n.order == 0) {
    cat("\n", "No order effect")
  }
  cat("\n", length(lambda), "different tuning parameters", 
    "\n")
  cat("\n Cross-validation criterion:", cv.crit, "\n")
  
  cat("---", "\n")
  
  cat("Parameter estimates after", x$folds, "-", "fold cross-validation", 
    "\n")
  
  cat("\n")
  coefs <- x$coefs.repar[which.min(x$criterion), ]
  
  theta <- intercepts <- order.effects <- gamma.X <- gamma.Z1 <- gamma.Z2 <- c()
  
  if (n.theta > 0) {
    cat("Thresholds:", "\n")
    theta <- coefs[1:n.theta]
    names(theta) <- paste0("theta", 1:n.theta)
    print(theta, ...)
    cat("\n")
  }
  
  if (n.order > 0) {
    cat(paste0(x$control$name.order, ":"), "\n")
    orders <- coefs[(n.theta + 1):(n.theta + n.order)]
    if (n.order == m) {
      names(orders) <- labels
    }
    if (n.order == 1) {
      names(orders) <- NULL
    }
    print(orders, ...)
    cat("\n")
  }else{
    orders <- NULL
  }
  
  if (n.intercepts > 0) {
    cat("Intercepts:", "\n")
    intercepts <- coefs[(n.theta + n.order + 1):(n.theta + 
      n.order + n.intercepts)]
    names(intercepts) <- labels
    print(intercepts, ...)
    cat("\n")
  }
  
  if (p.X > 0) {
    cat("Object-specific effects for subject-specific covariate(s):", 
      "\n")
    gamma.X <- matrix(coefs[(n.theta + n.order + n.intercepts + 
      1):(n.theta + n.order + n.intercepts + p.X * m)], 
      nrow = p.X, byrow = TRUE)
    if (rescale) {
      gamma.X <- t(t(gamma.X)/rep(x$design$sd.X, each = m))
    }
    colnames(gamma.X) <- labels
    rownames(gamma.X) <- vars.X
    print(gamma.X, ...)
    cat("\n")
  }
  
  if (p.Z1 > 0) {
    cat("Object-specific effects for subject-object-specific covariate(s):", 
      "\n")
    gamma.Z1 <- matrix(coefs[(n.theta + n.order + n.intercepts + 
      p.X * m + 1):(n.theta + n.order + n.intercepts + 
      p.X * m + p.Z1 * m)], nrow = p.Z1, byrow = TRUE)
    if (rescale) {
      gamma.Z1 <- t(t(gamma.Z1)/rep(x$design$sd.Z1, each = m))
    }
    colnames(gamma.Z1) <- labels
    rownames(gamma.Z1) <- vars.Z1
    print(gamma.Z1, ...)
    cat("\n")
  }
  
  if (p.Z2 > 0) {
    
    cat("Global effects for (subject-)object-specific covariate(s):", 
      "\n")
    gamma.Z2 <- coefs[(n.theta + n.order + n.intercepts + 
      p.X * m + p.Z1 * m + 1):(n.theta + n.order + n.intercepts + 
      p.X * m + p.Z1 * m + p.Z2)]
    if (rescale) {
      gamma.Z2 <- t(t(gamma.Z2)/x$design$sd.Z2)
    }
    names(gamma.Z2) <- vars.Z2
    print(gamma.Z2, ...)
    cat("\n")
  }
  
  cat("---", "\n")
  
  cat("\n")
  
  cat("Optimal lambda:", x$lambda[which.min(x$criterion)], 
    "\n")
  
  cat("\n")
  
  cat("Log likelihood:", x$logLik[which.min(x$criterion)], 
    "\n")
  
  
  coef.opt <- list(theta = theta, intercepts = intercepts, 
    order.effects = orders, gamma.X = gamma.X, gamma.Z1 = gamma.Z1, 
    gamma.Z2 = gamma.Z2)
  
  invisible(coef.opt)
  
  
}

