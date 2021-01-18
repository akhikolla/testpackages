#' Plot bootstrap intervals for BTLLasso
#' 
#' Plots bootstrap intervals for every single coefficient based on bootstrap estimates 
#' calculated by \code{\link{boot.BTLLasso}}. Bootstrap
#' intervals are separated by covariates, every covariate is plotted
#' separately. 
#' 
#' @param x boot.BTLLasso object
#' @param quantiles Which empirical quantiles of the bootstrap estimates should be plotted?
#' @param plots_per_page Number of plots per page, internally specified by \code{par(mfrow=...)}.
#' @param ask_new If TRUE, the user is asked before each plot.
#' @param rescale Should the parameter estimates be rescaled for plotting? Only 
#' applies if \code{scale = TRUE} was specified in \code{cv.BTLLasso}.
#' @param which Integer vector to specify which parameters/variables to plot. 
#' @param include.zero Should all plots contain zero?
#' @param rows Optional argument for the number of rows in the plot. 
#' Only applies if \code{plots_per_page>1}.
#' @param subs.X Optional vector of subtitles for variables in \code{X}. Can be used
#' to note the encoding of the single covariates, especially for dummy
#' variables.
#' @param subs.Z1 Optional vector of subtitles for variables in \code{Z1}. Can be used
#' to note the encoding of the single covariates, especially for dummy
#' variables.
#' @param main.Z2 Optional character containg main for plot
#' containing intervals for Z2 parameters. 
#' @param ... other parameters to be passed through to plot function.
#' @author Gunther Schauberger\cr \email{gunther.schauberger@@tum.de}
#' @seealso \code{\link{boot.BTLLasso}}, \code{\link{BTLLasso}},
#' \code{\link{cv.BTLLasso}}
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
plot.boot.BTLLasso <- function(x, quantiles = c(0.025, 0.975), 
                               plots_per_page = 1,  ask_new = TRUE, rescale = FALSE, 
                               which = "all", include.zero = TRUE, rows = NULL, 
                               subs.X = NULL, subs.Z1 = NULL,
                               main.Z2 = "Obj-spec. Covariates", ...){

  op <- par(no.readonly = TRUE)
  
  ## extract important things from the cv.BTLLasso object
  model <- x$cv.model
  epsilon <- model$control$epsilon
  accuracy <- -log10(epsilon)
  covariates <- c(model$design$vars.X, model$design$vars.Z1, 
    model$design$vars.Z2)
  
  ## running index for current plot
  index.num <- 1
  
  ## create matrix containing bootstrap intervals
  conf.ints <- apply(x$estimatesBrepar, 2, quantile, probs = quantiles,
                     type = 1, na.rm = TRUE)
  
  ## some more important parameters
  m <- model$Y$m
  labels <- model$Y$object.names
  
  n.theta <- model$design$n.theta
  n.order <- model$design$n.order
  n.intercepts <- model$design$n.intercepts
  if (n.intercepts > 0) {
    n.intercepts <- n.intercepts + 1
  }
  p.X <- model$design$p.X
  p.Z1 <- model$design$p.Z1
  p.Z2 <- model$design$p.Z2
  
  estimates <- model$coefs.repar[which.min(model$criterion),]
  estimates <- round(estimates, accuracy)
  conf.ints <- round(conf.ints, accuracy)
  
  p <- p.global <- 0
  gamma <- gamma.ci <- global <- global.ci <- c()
  covar <- covar.global <- c()
  
  start <- n.theta + 1
  
  all.subs <- c()
  
  if (n.order > 0) {
    end <- start + n.order - 1
    if (n.order == 1) {
      global <- c(global, estimates[start:end])
      global.ci <- cbind(global.ci, conf.ints[, start:end])
      covar.global <- c(covar.global, model$control$name.order)
      p.global <- p.global + 1
    }
    if (n.order > 1) {
      gamma <- c(gamma, estimates[start:end])
      gamma.ci <- cbind(gamma.ci, conf.ints[, start:end])
      covar <- c(covar, model$control$name.order)
      p <- p + 1
      all.subs <- c(all.subs, "")
    }
  }
  
  start <- n.theta + n.order + 1
  
  if (n.intercepts > 0) {
    end <- start + n.intercepts - 1
    covar <- c(covar, "Intercept")
    gamma <- c(gamma, estimates[start:end])
    gamma.ci <- cbind(gamma.ci, conf.ints[, start:end])
    p <- 1 + p
    all.subs <- c(all.subs, "")
  }
  
  start <- n.theta + n.order + n.intercepts + 1
  
  if (p.X > 0) {
    end <- start + p.X * m - 1
    covar <- c(covar, model$design$vars.X)
    if (rescale) {
      est <- estimates[start:end]/rep(model$design$sd.X, 
        each = m)
      est.ci <- t(t(conf.ints[, start:end])/rep(model$design$sd.X, 
        each = m))
    } else {
      est <- estimates[start:end]
      est.ci <- conf.ints[, start:end]
    }
    gamma <- c(gamma, est)
    gamma.ci <- cbind(gamma.ci, est.ci)
    p <- p + p.X
    if (is.null(subs.X)) {
      subs.X <- rep("", p.X)
    }
    all.subs <- c(all.subs, subs.X)
  }
  
  start <- n.theta + n.order + n.intercepts + p.X * m + 1
  
  if (p.Z1 > 0) {
    end <- start + p.Z1 * m - 1
    covar <- c(covar, model$design$vars.Z1)
    if (rescale) {
      est <- estimates[start:end]/rep(model$design$sd.Z1, 
        each = m)
      est.ci <- t(t(conf.ints[, start:end])/rep(model$design$sd.Z1, 
        each = m))
    } else {
      est <- estimates[start:end]
      est.ci <- conf.ints[, start:end]
    }
    gamma <- c(gamma, est)
    gamma.ci <- cbind(gamma.ci, est.ci)
    p <- p + p.Z1
    if (is.null(subs.Z1)) {
      subs.Z1 <- rep("", p.Z1)
    }
    all.subs <- c(all.subs, subs.Z1)
  }
  
  
  start <- n.theta + n.order + n.intercepts + p.X * m + p.Z1 * 
    m + 1
  
  if (p.Z2 > 0) {
    end <- start + p.Z2 - 1
    covar.global <- c(covar.global, model$design$vars.Z2)
    if (rescale) {
      est <- estimates[start:end]/model$design$sd.Z2
      est.ci <- t(t(conf.ints[, start:end])/model$design$sd.Z2)
    } else {
      est <- estimates[start:end]
      est.ci <- conf.ints[, start:end]
    }
    global <- c(global, est)
    global.ci <- cbind(global.ci, est.ci)
    p.global <- p.global + p.Z2
  }
  
  p.tot <- p
  if (p.global > 0) {
    p.tot <- p.tot + 1
  }
  
  
  suppressWarnings(if(which=="all"){
    which <- 1:p.tot
  })
  pages <- ceiling(length(which)/plots_per_page)
  
  if (is.null(rows)) {
    rows <- floor(sqrt(plots_per_page))
  } 
  
  cols <- ceiling(plots_per_page/rows)
  
  plots_on_page <- 0
  pages_done <- 0
  par(mfrow=c(rows, cols))

  index <- 1
  for (i in 1:p) {
    if(i %in% which){
    xlim <- range(gamma.ci[, index:(index + m - 1)])
    if (include.zero) {
      xlim <- range(0, xlim)
    }
    plot(gamma[index:(index + m - 1)], 1:m, xlim = xlim, 
      pch = 16, yaxt = "n", xlab = "", ylab = "", main = "", ...)
    
    segments(y0 = 1:m, x0 = gamma.ci[1, index:(index + m - 
      1)], x1 = gamma.ci[2, index:(index + m - 1)])
    axis(2, at = 1:m, labels = labels, las = 2)
    title(covar[i], line = 1.2)
    mtext(all.subs[i], side = 3, line = 0.2, cex = par()$cex)
   
    segments( 0, 1, 0, m , col="lightgray",lty=2,lwd=par()$lwd)
      
    plots_on_page <- plots_on_page+1
    if(plots_on_page==plots_per_page & pages_done<(pages-1)){
      plots_on_page <- 0
      pages_done <- pages_done+1
      if(interactive() & ask_new)
      {readline("Press enter for next plot!")}
      par(mfrow=c(rows, cols))
    }
    
    }
    index <- index + m
    
  }
  
  if (p.global > 0 & (p+1) %in% which) {
    xlim <- range(global.ci)
    if (include.zero) {
      xlim <- range(0, xlim)
    }
    plot(global, 1:p.global, xlim = xlim, pch = 16, yaxt = "n", 
      xlab = "", ylab = "", main = "Global Parameters", ...)
    segments(y0 = 1:p.global, x0 = global.ci[1, ], x1 = global.ci[2,])
    axis(2, at = 1:p.global, labels = covar.global, las = 2)
    segments( 0, 1, 0, p.global , col="lightgray",lty=2,lwd=par()$lwd)
  }
  
  par(op)
  invisible(conf.ints)
}

