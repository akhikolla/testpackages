#' Plot parameter paths for BTLLasso
#' 
#' Plots single paths for every parameter of a \code{BTLLasso} object or a \code{cv.BTLLasso}
#' object. In contrast, to \code{\link{paths}}, one plot per covariate is
#' created, every single parameter is illustrated by one path. For \code{cv.BTLLasso}
#' objects, the optimal model according to the cross-validation is marked by a
#' vertical dashed line.
#' 
#' @param x BTLLasso or cv.BTLLasso object
#' @param plots_per_page Number of plots per page, internally specified by \code{par(mfrow=...)}.
#' @param ask_new If TRUE, the user is asked before each plot.
#' @param rescale Should the parameter estimates be rescaled for plotting? Only 
#' applies if \code{scale = TRUE} was specified in \code{BTLLasso} or \code{cv.BTLLasso}.
#' @param which Integer vector to specify which parameters/variables to plot. 
#' @param equal.ranges Should all single plots (for different covariates) have
#' equal ranges on the y-axes. FALSE by default.
#' @param x.axis Should the paths be plotted against log(lambda+1) or against lambda?
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
#' @param ... Further plot arguments.
#' @author Gunther Schauberger\cr \email{gunther.schauberger@@tum.de}
#' @seealso \code{\link{BTLLasso}}, \code{\link{cv.BTLLasso}},
#' \code{\link{paths}}
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
#' @keywords BTLLasso paths parameter paths
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
plot.BTLLasso <- function(x, plots_per_page = 1, ask_new = TRUE, 
  rescale = FALSE, which = "all",
  equal.ranges = FALSE, x.axis = c("loglambda", "lambda"), 
  rows = NULL, subs.X = NULL, subs.Z1 = NULL,
  main.Z2 = "Obj-spec. Covariates", ...) {
  
  op <- par(no.readonly = TRUE)
  
  if(length(x$lambda)==1){
    stop("Only one tuning parameter, nothing to plot!")
  }
  
  ## get correct x axis
  x.axis <- match.arg(x.axis)
  
  
  if (x.axis == "lambda") {
    norm <- x$lambda
    norm.range <- rev(range(norm))
    
    x.axis.name <- expression(lambda)
  }
  
  if (x.axis == "loglambda") {
    norm <- log(x$lambda + 1)
    norm.range <- rev(range(norm))
    
    x.axis.name <- expression(log(lambda + 1))
  }
  
  ## get basic parameters
  m <- x$Y$m
  n.theta <- x$design$n.theta
  n.order <- x$design$n.order
  n.intercepts <- x$design$n.intercepts
  if (n.intercepts > 0) {
    n.intercepts <- n.intercepts + 1
  }
  
  p.X <- x$design$p.X
  p.Z1 <- x$design$p.Z1
  p.Z2 <- x$design$p.Z2
  
  ## get object labels
  labels <- x$Y$object.names
  
  ## get all coefficients
  coefs <- x$coefs.repar
  
  ## get range of all coefs if necessary
  y.range <- NA
  if(equal.ranges){
    y.range <- range(coefs)
  }
  
  ## initialize matrix of all parameters to plot
  coef.plot <- c()
  ## index.plots  will contain number indicating rows for equal plots
  index.plots <- c()
  ## running index for current plot
  index.num <- 1
  ## will contain all labels for right-hand labels
  all.labs <- c()
  ## will contain all mains
  all.mains <- c()
  all.subs <- c()
  
  ## go through all possible model components 
  if (n.order > 0) {
    order.effects <- coefs[, (n.theta + 1):(n.theta + n.order)]
    coef.plot <- cbind(coef.plot, order.effects)
    index.plots <- c(index.plots, rep(index.num,n.order))
    index.num <- index.num+1
    if(n.order>1){
      all.labs <- c(all.labs,labels)
    }else{
      all.labs <- c(all.labs,"")
    }
    all.mains <- c(all.mains, x$control$name.order)
    all.subs <- c(all.subs, "")
  }
  
  if (n.intercepts > 0) {
    intercepts <- coefs[, (n.theta + n.order + 1):(n.theta + 
      n.order + n.intercepts), drop = FALSE]
    coef.plot <- cbind(coef.plot, intercepts)
    index.plots <- c(index.plots, rep(index.num, n.intercepts))
    index.num <- index.num+1
    all.labs <- c(all.labs,labels)
    all.mains <- c(all.mains, "Intercepts")
    all.subs <- c(all.subs, "")
  }
  
  if (p.X > 0) {
    gamma.X <- coefs[, (n.theta + n.order + n.intercepts + 
      1):(n.theta + n.order + n.intercepts + p.X * m), 
      drop = FALSE]
    if (rescale) {
      gamma.X <- t(t(gamma.X)/rep(x$design$sd.X, each = m))
    }
    coef.plot <- cbind(coef.plot, gamma.X)
    index.plots <- c(index.plots, rep(index.num:(index.num+p.X-1), each=m))
    index.num <- index.num+p.X
    all.labs <- c(all.labs,rep(labels,p.X))
    all.mains <- c(all.mains, x$design$vars.X)
    all.subs <- c(all.subs, subs.X)
  }
  
  if (p.Z1 > 0) {
    gamma.Z1 <- coefs[, (n.theta + n.order + n.intercepts + 
      p.X * m + 1):(n.theta + n.order + n.intercepts + 
      p.X * m + p.Z1 * m), drop = FALSE]
    if (rescale) {
      gamma.Z1 <- t(t(gamma.Z1)/rep(x$design$sd.Z1, 
        each = m))
    }
    coef.plot <- cbind(coef.plot, gamma.Z1)
    index.plots <- c(index.plots, rep(index.num:(index.num+p.Z1-1), each=m))
    index.num <- index.num+p.Z1
    all.labs <- c(all.labs,rep(labels,p.Z1))
    all.mains <- c(all.mains, x$design$vars.Z1)
    all.subs <- c(all.subs, subs.Z1)
  }
  
  if (p.Z2 > 0) {
    gamma.Z2 <- coefs[, (n.theta + n.order + n.intercepts + 
      p.X * m + p.Z1 * m + 1):(n.theta + n.order + n.intercepts + 
      p.X * m + p.Z1 * m + p.Z2), drop = FALSE]
    if (rescale) {
      gamma.Z2 <- t(t(gamma.Z2)/x$design$sd.Z2)
    }
    coef.plot <- cbind(coef.plot, gamma.Z2)
    index.plots <- c(index.plots, rep(index.num, p.Z2))
    index.num <- index.num+1
    all.labs <- c(all.labs,x$design$vars.Z2)
    all.mains <- c(all.mains, main.Z2)
    all.subs <- c(all.subs, "")
  }
  
  n.plots <- index.num-1
  
  suppressWarnings(if(which=="all"){
    which <- 1:n.plots
  })
  pages <- ceiling(length(which)/plots_per_page)
  
  if (is.null(rows)) {
    rows <- floor(sqrt(plots_per_page))
  } 
  
  cols <- ceiling(plots_per_page/rows)
  
  plots_on_page <- 0
  pages_done <- 0
  par(mfrow=c(rows, cols))
  
  for(u in 1:n.plots){
    if(u %in% which){
      
      plot.comp(u, norm, coef.plot, index.plots,
                all.labs, all.mains, all.subs, y.range,
                x.axis.name, x$criterion, norm.range, ...)
        
      
    plots_on_page <- plots_on_page+1
      if(plots_on_page==plots_per_page & pages_done<(pages-1)){
        plots_on_page <- 0
        pages_done <- pages_done+1
        if(interactive() & ask_new)
        {readline("Press enter for next plot!")}
        par(mfrow=c(rows, cols))
      }
    }
  }
  
  
  par(op)
}


plot.comp <- function(u, norm, coef.plot, index.plots,
                      all.labs, all.mains, all.subs, y.range,
                      x.axis.name, criterion, norm.range, ...){
  
  
  if (!is.null(criterion)) {
    x.axis.min <- norm[which.min(criterion)]
  }
  
  index.u <- which(index.plots==u)
  l.u <- length(index.u)

  
  cur.coef <- coef.plot[,index.u,drop=FALSE]
  final.u <- cur.coef[nrow(cur.coef),]
  
  if(is.na(y.range)){
    y.range.u <- range(cur.coef)
  }else{
    y.range.u <- y.range
  }

      plot(norm, cur.coef[,1], ylim = y.range.u, type = "l", 
       main = "", ylab = "estimates", xlab = x.axis.name, 
       xlim = norm.range,  frame.plot = FALSE,
       lwd=par()$lwd, ...)
    
  if(l.u>1){
    for (uu in 2:l.u) {
      lines(norm, cur.coef[, uu], lwd=par()$lwd)
    }
  }
  
  title(main = all.mains[u], line = 1.2)
  mtext(all.subs[u], side = 3, line = 0.2, cex = par()$cex)
  
  if (!is.null(criterion)) {
    segments( x.axis.min, min(y.range.u),
              x.axis.min, max(y.range.u) ,
              col=2,lty=2,lwd=par()$lwd)
    
    }

  
  x.lab1 <- norm[length(norm)]-abs(diff(range(norm)))*0.02
  x.lab2 <- norm[length(norm)]-abs(diff(range(norm)))*0.005
  y.lab1 <- final.u
  y.lab2 <- spread.labs(y.lab1, 1.2*strheight("A"))
  
  text( x.lab1, y.lab2, all.labs[index.u],pos=4)
  segments( x.lab2, y.lab1,
            x.lab1, y.lab2 ,col="gray")
  
} 