#' Plot covariate paths for BTLLasso
#' 
#' Plots paths for every covariate of a BTLLasso object or a cv.BTLLasso
#' object. In contrast to \code{\link{plot.BTLLasso}}, only one plot is created,
#' every covariate is illustrated by one path. For \code{cv.BTLLasso} objects, the
#' optimal model according to the cross-validation is marked by a vertical
#' dashed line.
#' 
#' @param model \code{BTLLasso} or \code{cv.BTLLasso} object
#' @param y.axis Two possible values for the y-axis. Variables can either be plotted
#' with regard to their contribution to the total penalty term (\code{y.axis='penalty'}) or
#' with regard to the $L_2$ norm of the corresponding parameter vector (\code{y.axis='L2'}).
#' @param x.axis Should the paths be plotted against log(lambda+1) or against lambda?' 
#' @author Gunther Schauberger\cr \email{gunther.schauberger@@tum.de}
#' @seealso \code{\link{BTLLasso}}, \code{\link{cv.BTLLasso}},
#' \code{\link{plot.BTLLasso}}
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
#' @keywords BTLLasso covariate paths
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
paths <- function(model, y.axis = c("penalty", "L2"), x.axis = c("loglambda", 
  "lambda")) {

  
  x.axis <- match.arg(x.axis)
  y.axis <- match.arg(y.axis)
  
  if (y.axis == "penalty") {
    y.text <- "penalty size"
  }
  if (y.axis == "L2") {
    y.text <- "L2 norm"
  }
  
  coefs <- model$coefs
  covar <- c(model$design$vars.X, model$design$vars.Z1, model$design$vars.Z2)
  

  if (x.axis == "lambda") {
    norm <- model$lambda
    norm.range <- rev(range(norm))
    x.axis.name <- expression(lambda)
  }
  
  if (x.axis == "loglambda") {
    norm <- log(model$lambda + 1)
    norm.range <- rev(range(norm))
    x.axis.name <- expression(log(lambda + 1))
  }
  
  
  
  m <- model$Y$m
  n.theta <- model$design$n.theta
  n.order <- model$design$n.order
  n.intercepts <- model$design$n.intercepts
  acoefs <- model$penalty$acoefs
  
  p.X <- model$design$p.X
  p.Z1 <- model$design$p.Z1
  p.Z2 <- model$design$p.Z2
  
  
  numpen.order <- model$penalty$numpen.order
  numpen.intercepts <- model$penalty$numpen.intercepts
  numpen.X <- model$penalty$numpen.X
  numpen.Z1 <- model$penalty$numpen.Z1
  numpen.Z2 <- model$penalty$numpen.Z2
  
  labels <- model$Y$object.names
  
  criterion <- model$criterion
  
  order.effects <- intercepts <- gamma.X <- gamma.Z1 <- gamma.Z2 <- c()
  
  index.cols.X <- index.cols.Z1 <- index.cols.Z2 <- c()
  index.rows.X <- index.rows.Z1 <- index.rows.Z2 <- c()
  
  if (n.order > 0) {
    order.effects <- coefs[, (n.theta + 1):(n.theta + n.order)]
  }
  
  if (n.intercepts > 0) {
    intercepts <- coefs[, (n.theta + n.order + 1):(n.theta + 
      n.order + n.intercepts)]
  }
  
  p <- p.X + p.Z1 + p.Z2
  
  paths <- c()
  
  start.row <- n.theta + n.intercepts + n.order
  if (p.X > 0) {
    index <- rep((1:p.X), each = m - 1)
    for (i in 1:p.X) {
      if (y.axis == "penalty") {
        paths <- cbind(paths, rowSums(abs(coefs[, start.row + 
          which(index == i), drop = FALSE] %*% acoefs[start.row + 
          which(index == i), , drop = FALSE])))
      } else {
        paths <- cbind(paths, sqrt(rowSums(coefs[, start.row + 
          which(index == i), drop = FALSE]^2)))
      }
      
    }
    start.row <- start.row + length(index)
  }
  
  if (p.Z1 > 0) {
    index <- rep(1:p.Z1, each = m)
    for (i in 1:p.Z1) {
      if (y.axis == "penalty") {
        paths <- cbind(paths, rowSums(abs(coefs[, start.row + 
          which(index == i), drop = FALSE] %*% acoefs[start.row + 
          which(index == i), , drop = FALSE])))
      } else {
        paths <- cbind(paths, sqrt(rowSums(coefs[, start.row + 
          which(index == i), drop = FALSE]^2)))
      }
      
    }
    start.row <- start.row + length(index)
  }
  
  if (p.Z2 > 0) {
    index <- 1:p.Z2
    for (i in 1:p.Z2) {
      if (y.axis == "penalty") {
        paths <- cbind(paths, rowSums(abs(coefs[, start.row + 
          which(index == i), drop = FALSE] %*% acoefs[start.row + 
          which(index == i), , drop = FALSE])))
      } else {
        paths <- cbind(paths, sqrt(rowSums(coefs[, start.row + 
          which(index == i), drop = FALSE]^2)))
      }
      
    }
  }
  
  
  
  if (!is.null(criterion)) {
    x.axis.min <- norm[which.min(criterion)]
  }
  
  
  
  if (numpen.intercepts > 0) {
    if (y.axis == "penalty") {
      paths <- cbind(rowSums(abs(intercepts %*% acoefs[(n.theta + 
        n.order + 1):(n.theta + n.order + n.intercepts), 
        (numpen.order + 1):(numpen.order + numpen.intercepts)])), 
        paths)
    } else {
      paths <- cbind(sqrt(rowSums(intercepts^2)), paths)
    }
    
    p <- p + 1
    covar <- c("Intercept", covar)
  }
  
  if (numpen.order > 0) {
    if (y.axis == "penalty") {
      paths <- cbind(rowSums(abs(order.effects %*% acoefs[(n.theta + 
        1):(n.theta + n.order), 1:numpen.order])), paths)
    } else {
      paths <- cbind(sqrt(rowSums(order.effects^2)), paths)
    }
    p <- p + 1
    covar <- c(model$control$name.order, covar)
  }
  
  
  plot(norm, paths[, 1], type = "l", ylim = range(paths), ylab = y.text, 
    xlab = x.axis.name, xlim = norm.range, las = 1,lwd=par()$lwd,
    frame.plot = FALSE)
  for (o in 2:p) {
    lines(norm, paths[, o],lwd=par()$lwd)
  }
  
  if (!is.null(criterion)) {
    segments( x.axis.min, min(paths),
              x.axis.min, max(paths) ,
              col=2,lty=2,lwd=par()$lwd)
    
  }
  
  
  x.lab1 <- norm[length(norm)]-abs(diff(range(norm)))*0.02
  x.lab2 <- norm[length(norm)]-abs(diff(range(norm)))*0.005
  y.lab1 <- paths[nrow(paths), ]
  y.lab2 <- spread.labs(y.lab1, 1.2*strheight("A"))
  
  text( x.lab1, y.lab2, covar ,pos=4)
  segments( x.lab2, y.lab1,
            x.lab1, y.lab2 ,col="gray")
  
  
  
}
