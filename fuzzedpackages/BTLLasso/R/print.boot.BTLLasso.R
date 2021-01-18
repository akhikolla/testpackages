#' Print function for boot.BTLLasso objects
#' 
#' Prints the most important output of \code{boot.BTLLasso} objects.
#' 
#' @method print boot.BTLLasso
#' @param x \code{boot.BTLLasso} object
#' @param quantiles Which empirical quantiles of the bootstrap estimates should be printed?
#' @param rescale Should the parameter estimates be rescaled for plotting? Only 
#' applies if \code{scale = TRUE} was specified in \code{BTLLasso} or \code{cv.BTLLasso}.
#' @param \dots possible further arguments for print command
#' @author Gunther Schauberger\cr \email{gunther.schauberger@@tum.de}
#' @seealso \code{\link{boot.BTLLasso}}
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
print.boot.BTLLasso <- function(x, quantiles = c(0.025, 0.975), rescale = FALSE, ...) {
  
  model <- x$cv.model
  epsilon <- model$control$epsilon
  accuracy <- -log10(epsilon)
  covariates <- c(model$design$vars.X, model$design$vars.Z1, 
    model$design$vars.Z2)
  
  # conf.ints <- apply(x$estimatesB, 2, quantile, probs = quantiles,
  #   type = 1, na.rm = TRUE)
  conf.ints <- apply(x$estimatesBrepar, 2, quantile, probs = quantiles,
    type = 1, na.rm = TRUE)

  
  m <- model$Y$m
  labels <- model$Y$x.names
  
  n.theta <- model$design$n.theta
  n.order <- model$design$n.order
  n.intercepts <- model$design$n.intercepts
  if (n.intercepts > 0) {
    n.intercepts <- n.intercepts + 1
  }
  p.X <- model$design$p.X
  p.Z1 <- model$design$p.Z1
  p.Z2 <- model$design$p.Z2
  
  estimates <- model$coefs.repar[which.min(model$criterion), 
    ]
  estimates <- round(estimates, accuracy)
  conf.ints <- round(conf.ints, accuracy)
  
  gamma.total <- c()
  
  start <- 1
  end <- n.theta
  
  cat("Bootstrap intervals:\n")
  
  cat("---", "\n")
  
  if (n.theta > 0) {
    cat("Thresholds:", "\n")
    gamma <- matrix(NA, nrow = 3, ncol = n.theta)
    rownames(gamma)[c(1, 3)] <- rownames(conf.ints)
    rownames(gamma)[2] <- "estimate"
    gamma[c(1, 3), ] <- conf.ints[, start:end]
    gamma[2, ] <- estimates[start:end]
    colnames(gamma) <- names(estimates[start:end])
    print(gamma, ...)
    cat("\n")
    gamma.total <- cbind(gamma.total, gamma)
  }
  
  start <- n.theta + 1
  
  if (n.order > 0) {
    end <- start + n.order - 1
    gamma <- matrix(NA, nrow = 3, ncol = n.order)
    rownames(gamma)[c(1, 3)] <- rownames(conf.ints)
    rownames(gamma)[2] <- "estimate"
    cat(paste0(model$control$name.order, ":"), "\n")
    gamma[c(1, 3), ] <- conf.ints[, start:end]
    gamma[2, ] <- estimates[start:end]
    colnames(gamma) <- names(estimates[start:end])
    print(gamma, ...)
    cat("\n")
    gamma.total <- cbind(gamma.total, gamma)
  }
  
  start <- n.theta + n.order + 1
  
  if (n.intercepts > 0) {
    end <- start + n.intercepts - 1
    gamma <- matrix(NA, nrow = 3, ncol = n.intercepts)
    rownames(gamma)[c(1, 3)] <- rownames(conf.ints)
    rownames(gamma)[2] <- "estimate"
    cat("Intercepts:\n")
    gamma[c(1, 3), ] <- conf.ints[, start:end]
    gamma[2, ] <- estimates[start:end]
    colnames(gamma) <- names(estimates[start:end])
    print(gamma, ...)
    cat("\n")
    gamma.total <- cbind(gamma.total, gamma)
  }
  
  start <- n.theta + n.order + n.intercepts + 1
  
  if (p.X > 0) {
    end <- start + p.X * m - 1
    if (rescale) {
      est <- estimates[start:end]/rep(model$design$sd.X, 
        each = m)
      est.ci <- t(t(conf.ints[, start:end])/rep(model$design$sd.X, 
        each = m))
    } else {
      est <- estimates[start:end]
      est.ci <- conf.ints[, start:end]
    }
    
    gamma <- matrix(NA, nrow = 3, ncol = p.X * m)
    rownames(gamma)[c(1, 3)] <- rownames(conf.ints)
    rownames(gamma)[2] <- "estimate"
    cat("Object-specific effects for subject-specific covariate(s):", 
      "\n")
    gamma[c(1, 3), ] <- est.ci
    gamma[2, ] <- est
    colnames(gamma) <- names(estimates[start:end])
    print(gamma, ...)
    cat("\n")
    gamma.total <- cbind(gamma.total, gamma)
  }
  
  start <- n.theta + n.order + n.intercepts + p.X * m + 1
  
  
  
  if (p.Z1 > 0) {
    end <- start + p.Z1 * m - 1
    if (rescale) {
      est <- estimates[start:end]/rep(model$design$sd.Z1, 
        each = m)
      est.ci <- t(t(conf.ints[, start:end])/rep(model$design$sd.Z1, 
        each = m))
    } else {
      est <- estimates[start:end]
      est.ci <- conf.ints[, start:end]
    }
    gamma <- matrix(NA, nrow = 3, ncol = p.Z1 * m)
    rownames(gamma)[c(1, 3)] <- rownames(conf.ints)
    rownames(gamma)[2] <- "estimate"
    cat("Object-specific effects for subject-object-specific covariate(s):", 
      "\n")
    gamma[c(1, 3), ] <- est.ci
    gamma[2, ] <- est
    colnames(gamma) <- names(estimates[start:end])
    print(gamma, ...)
    cat("\n")
    gamma.total <- cbind(gamma.total, gamma)
  }
  
  
  start <- n.theta + n.order + n.intercepts + p.X * m + p.Z1 * 
    m + 1
  
  if (p.Z2 > 0) {
    end <- start + p.Z2 - 1
    if (rescale) {
      est <- estimates[start:end]/model$design$sd.Z2
      est.ci <- t(t(conf.ints[, start:end])/model$design$sd.Z2)
    } else {
      est <- estimates[start:end]
      est.ci <- conf.ints[, start:end]
    }
    gamma <- matrix(NA, nrow = 3, ncol = p.Z2)
    rownames(gamma)[c(1, 3)] <- rownames(conf.ints)
    rownames(gamma)[2] <- "estimate"
    cat("Global effects for (subject-)object-specific covariate(s):", 
      "\n")
    gamma[c(1, 3), ] <- est.ci
    gamma[2, ] <- est
    colnames(gamma) <- names(estimates[start:end])
    print(gamma, ...)
    gamma.total <- cbind(gamma.total, gamma)
  }
  
  
  invisible(gamma.total)
  
}

