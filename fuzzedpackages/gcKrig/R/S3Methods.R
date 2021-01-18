#### S3 methods that can be applied to the results of three classes:
#### simgc (plot only), mlegc (summary, print, plot, vcov, profile) and predgc (summary and plot)


#### S3 Method on Simulated Datasets --------------------------------------------------------------------
#### ----------------------------------------------------------------------------------------------------

plot.simgc <- function(
    x, index = 1, plottype = 'Text', xlab = "xloc", ylab = "yloc", xlim = NULL, ylim = NULL,
    pch = 20, textcex = 0.8, plotcex = 1, angle = 60, col = 4, col.regions = gray(90:0/100),...
){
  X <- x
  rm(x) # to avoid name conflicts
  if(is.null(ncol(X$data))){
    Data <- X$data
  }else{
  Data <- X$data[index,]
  }
  r1 <- apply(X$locs, 2, range)
  if(is.null(xlim)) xlim <- r1[,1]
  if(is.null(ylim)) ylim <- r1[,2]

  if (requireNamespace("latticeExtra", quietly = TRUE)) {
  ## Plot1
   plot1 <- lattice::levelplot(Data ~  X$locs[,1] + X$locs[,2], col.regions = col.regions,
                            xlab = xlab, ylab = ylab, cex = plotcex,
                            panel = function(x,y,z,...) {
                              lattice::panel.levelplot(x,y,z,...)
                              lattice::panel.text(X$locs[,1], X$locs[,2],  Data, col = col, cex = textcex)
                            })

  ## Plot2
  plot2 <- lattice::levelplot(Data ~ X$locs[,1] + X$locs[,2], panel = latticeExtra::panel.levelplot.points,
                     col.regions = col.regions, xlab = xlab, ylab = ylab, cex = plotcex)
  }else{
    stop("Unable to Provide level Plot. Please install {latticeExtra} first!")
}

  ## Plot3
  if(plottype == '3D'){
  if (requireNamespace("scatterplot3d", quietly = TRUE)) {
    scatterplot3d::scatterplot3d(x = X$locs[,1], y = X$locs[,2], z = Data, highlight.3d = TRUE,
                                 col.grid = "lightblue",  type='h', xlim = xlim, ylim = ylim,
                                 col.axis = "blue", xlab = xlab, ylab = ylab, angle = angle,  pch = pch,
                                 scale.y=1.2, y.margin.add = 0.2)
  }else{
    stop("Unable to Provide Three Dimensional Plot.
         Please install {scatterplot3d} first!")
  }
 }
  if(plottype == 'Text'){
    print(plot1)
  }
  if(plottype == 'Dot'){
    print(plot2)
  }
}



#### S3 Method on Maximum Likelihood Estimate Results ---------------------------------------------------
#### ----------------------------------------------------------------------------------------------------

plot.mlegc <- function(x, plotdata = "Observed", plottype = "2D",
                       xlab = "xloc", ylab = "yloc", xlim = NULL, ylim = NULL,
                       pch = 20, textcex = 0.8, plotcex = 1, angle = 60, col = 4,
                       col.regions = gray(90:0/100),...)
{
  fitted <- x
  rm(x)
  Fitmean <- fitted$args$marginal$fm$linkinv(fitted$MLE[1:ncol(fitted$x)]%*%
                                  t(fitted$x)*fitted$args$effort)

  Data <- fitted$args$y


  if (requireNamespace("latticeExtra", quietly = TRUE)) {
    ####  Plot1 Originl Count
    plot1 <- lattice::levelplot(Data ~ fitted$args$locs[,1] + fitted$args$locs[,2],
                             col.regions = col.regions, xlab = xlab, ylab = ylab, cex = plotcex,
                             main = NULL, panel = function(x, y, z, ...) {
                               lattice::panel.levelplot(x, y, z, ...)
                               lattice::panel.text(fitted$args$locs[,1], fitted$args$locs[,2],
                                          Data, col = col, cex = textcex)
                             })

    #### Plot2 Fitted Means
    plot2 <- lattice::levelplot(Fitmean ~ fitted$args$locs[,1] + fitted$args$locs[,2],
                             col.regions = col.regions, xlab = xlab, ylab = ylab, cex = plotcex,
                             main = NULL)

    #                         panel = function(...) {
    #                           lattice::panel.levelplot(...)
    #                         }))
  }else{
    stop("Unable to Provide level Plot. Please install {latticeExtra} first!")
  }

  if(plotdata == "Observed" & plottype == "2D") print(plot1)
  if(plotdata == "Fitted" & plottype == "2D") print(plot2)


   #### Plot3 3d Plot Original

  if (requireNamespace("scatterplot3d", quietly = TRUE)) {


    if(plotdata == "Observed" & plottype == "3D"){

    scatterplot3d::scatterplot3d(x = fitted$args$locs[,1], y = fitted$args$locs[,2],
                                 z = Data, color = col, col.grid = "lightblue",
                                 type='h', xlim = xlim, ylim = ylim, main = NULL,
                                 col.axis = "blue", xlab = xlab, ylab = ylab, angle = angle,  pch = pch,
                                 scale.y=1.2, y.margin.add = 0.2)

    }
    #### Plot4 3d Plot Fitted

    if(plotdata == "Fitted" & plottype == "3D"){

    scatterplot3d::scatterplot3d(x = fitted$args$locs[,1], y = fitted$args$locs[,2],
                                 z = Fitmean, color = col, col.grid = "lightblue",
                                 type='h', xlim = xlim, ylim = ylim, main = NULL,
                                 col.axis = "blue", xlab = xlab, ylab = ylab, angle = angle,  pch = pch,
                                 scale.y=1.2, y.margin.add = 0.2)
    }

  }else{
    stop("Unable to Provide Three Dimensional Plot.
         Please install {scatterplot3d} first!")
  }
}






print.mlegc <- function(x, digits = max(3, getOption("digits") - 3),...)
{
  fitted <- x
  rm(x)
  cat("\nCall:", deparse(fitted$call, width.cutoff = floor(getOption("width")*0.6)), "", sep = "\n")
  cat("Parameter Estimates:\n")
  print.default(format(fitted$MLE , digits = digits), print.gap = 2.5, quote = FALSE)
  cat("\n")
  cat("Estimated Standard Deviation:\n")
  vcov <- try(solve(fitted$hessian), silent = TRUE)
  if(inherits(vcov, "try-error")){
    if (requireNamespace("Matrix", quietly = TRUE)){
      vcov <- solve(Matrix::nearPD(fitted$hessian)$mat)
    }else{
      stop("Please Install {Matrix} first!")
   }
  }else{
  vmat <- diag(vcov)
  print.default(format(sqrt(vmat)) , digits = digits, print.gap = 2.5, quote = FALSE)
  }
  cat("\n")
  print(paste0("AIC: ", signif(fitted$AIC,4) ))
  print(paste0("Log-likelihood: ", signif(fitted$log.lik, 4) ))
}



vcov.mlegc <- function(object, digits = max(3, getOption("digits") - 3), ...)
{
  vcov <- try(solve(object$hessian), silent = TRUE)

  if(inherits(vcov, "try-error")){
    if (requireNamespace("Matrix", quietly = TRUE)){
      vcov <- solve(Matrix::nearPD(object$hessian)$mat)
      warning("Initial  Variance-Covariance Matrix is NOT Positive Definite, Approximation is used.")
    }else{
      stop("Please install {Matrix} first!")
    }
  }
    cat("Estimated Variance-Covariance Matrix is:\n")
    print(vcov)
    cat("\n")
    cat("Estimated Correlation Matrix is:\n")
    print(cov2cor(vcov))
}






profile.mlegc <- function(fitted, par.index, alpha = 0.05, start.point = NULL,
                            method = 'GQT', nrep = 1000, seed = 12345,...){

  eps1 <- .Machine$double.eps^(0.2)
  maxiter <- 10
  if (!method %in% c("GQT", "GHK"))
    stop("'method' must be GQT or GHK.")

  if(is.null(start.point)){
    vcov <- try(solve(fitted$hessian), silent = TRUE)
    if(inherits(vcov, "try-error")){
      stop("Please specify the starting points (left and right) for computing profile likelihoods!")
    }
    se <- ifelse(diag(vcov) > 0, sqrt(diag(vcov)), 0)
    if(par.index <= fitted$nreg){
    s01 = fitted$MLE[par.index] -se[par.index]*qnorm(1-alpha/2)
    s02 = fitted$MLE[par.index] +se[par.index]*qnorm(1-alpha/2)
    }else if(fitted$nug == 1 & par.index == fitted$par.df){
      s01 = eps1; s02 = fitted$MLE[par.index] + 2*eps1
    }else{
    s01 <- max(fitted$MLE[par.index] - eps1, fitted$optlb[par.index])
    s02 <- min(fitted$MLE[par.index] + eps1, fitted$optub[par.index])
    }
    start.point <- c(s01,s02)
  }
  ## begin iteration, with the lower bound first
  l <- rep(0,maxiter); u <- rep(0,maxiter)

  if(method == 'GQT'){
    if(fitted$MLE[par.index] <= fitted$optlb[par.index]){
      lb = fitted$optlb[par.index]
    }else{
    l.vec <- mleProfilelikGQT(fitted = fitted, par.index = par.index, fixvalue = start.point[1], single = FALSE,
                              start = fitted$MLE[-par.index], alpha = alpha)
    l[1] <- mleProfilelikGQT(fitted = fitted, par.index = par.index, fixvalue = l.vec, single = TRUE,
                             start = start.point[1], alpha = alpha)
    for(i in 1:maxiter){
      l.vec <- mleProfilelikGQT(fitted = fitted, par.index = par.index, fixvalue = l[i], single = FALSE,
                                start = fitted$MLE[-par.index], alpha = alpha)
      l[i+1] <- mleProfilelikGQT(fitted = fitted, par.index = par.index, fixvalue = l.vec, single = TRUE,
                                 start = l[i], alpha = alpha)
      if(abs(l[i+1]-l[i]) < abs(0.05*l[i])) break;
      if(i == 10) warning("Number of Iterations in profiling may not enough!")
    }
    lb = l[i+1]
  }

    #### Begin iteration, the upper bound
    if(fitted$MLE[par.index] == fitted$optub[par.index]){
      ub = fitted$optub[par.index]
    }else{
    u.vec <- mleProfilelikGQT(fitted = fitted, par.index = par.index, fixvalue = start.point[2], single = FALSE,
                              start = fitted$MLE[-par.index], alpha = alpha)
    u[1] <- mleProfilelikGQT(fitted = fitted, par.index = par.index, fixvalue = u.vec, single = TRUE,
                             start = start.point[2], alpha = alpha)
    for(j in 1:maxiter){
      u.vec = mleProfilelikGQT(fitted = fitted, par.index = par.index, fixvalue = u[j], single = FALSE,
                               start = fitted$MLE[-par.index], alpha = alpha)
      u[j+1] = mleProfilelikGQT(fitted = fitted, par.index = par.index, fixvalue = u.vec, single = TRUE,
                                start = u[j], alpha = alpha)
      if(abs(u[j+1]-u[j]) < abs(0.05*u[j])) break;
      if(j == 10) warning("Number of Iterations in profiling may not enough!")
     }
    ub = u[j+1]
   }
  }else if(method == 'GHK'){

    if(fitted$MLE[par.index] <= fitted$optlb[par.index]){
      lb = fitted$optlb[par.index]
    }else{
    l.vec <- mleProfilelikGHK(fitted = fitted, par.index = par.index, fixvalue = start.point[1], single = FALSE,
                              start = fitted$MLE[-par.index], alpha = alpha, nrep = nrep, seed = seed)
    l[1] <- mleProfilelikGHK(fitted = fitted, par.index = par.index, fixvalue = l.vec, single = TRUE,
                             start = start.point[1], alpha = alpha, nrep = nrep, seed = seed)
    for(i in 1:maxiter){
      l.vec <- mleProfilelikGHK(fitted = fitted, par.index = par.index, fixvalue = l[i], single = FALSE,
                                start = fitted$MLE[-par.index], alpha = alpha, nrep = nrep, seed = seed)
      l[i+1] <- mleProfilelikGHK(fitted = fitted, par.index = par.index, fixvalue = l.vec, single = TRUE,
                                 start = l[i], alpha = alpha, nrep = nrep, seed = seed)
      if(abs(l[i+1]-l[i]) < abs(0.05*(l[i]+.Machine$double.eps))) break;
      if(i == 10) warning("Number of Iterations in profiling may not enough!")
    }
    lb = l[i+1]
    }
    #### Begin iteration, the upper bound

    if(fitted$MLE[par.index] == fitted$optub[par.index]){
      ub = fitted$optub[par.index]
    }else{

    u.vec <- mleProfilelikGHK(fitted = fitted, par.index = par.index, fixvalue = start.point[2],
                              single = FALSE, start = fitted$MLE[-par.index],
                              alpha = alpha, nrep = nrep, seed = seed)
    u[1] <- mleProfilelikGHK(fitted = fitted, par.index = par.index, fixvalue = u.vec,
                             single = TRUE, start = start.point[2], alpha = alpha,
                             nrep = nrep, seed = seed)
    for(j in 1:maxiter){
      u.vec = mleProfilelikGHK(fitted = fitted, par.index = par.index, fixvalue = u[j], single = FALSE,
                               start = fitted$MLE[-par.index], alpha = alpha, nrep = nrep, seed = seed)
      u[j+1] = mleProfilelikGHK(fitted = fitted, par.index = par.index, fixvalue = u.vec, single = TRUE,
                                start = u[j], alpha = alpha, nrep = nrep, seed = seed)
      if(abs(u[j+1]-u[j]) < abs(0.05*u[j])) break;
      if(j == 10) warning("Number of Iterations in profiling may not enough!")
      }
     ub = u[j+1]
    }
  }
  percent <- 100*(1-alpha)
  return(cat("An approximate", percent, "percent profile likelihood-based confidence interval is (",
             format(round(lb,3)), ",", format(round(ub,3)),")."))
  return(c(lb,ub))
}


summary.mlegc <- function(object, ...)
{
#### The Summary Format is same as the Standard one in package "gcmr" and "betareg".
  cf <- object$MLE
  vcov <- try(solve(object$hessian), silent = TRUE)

  if(inherits(vcov, "try-error")){
    if (requireNamespace("Matrix", quietly = TRUE)){
      vcov <- solve(Matrix::nearPD(object$hessian)$mat)
    }else{
      stop("Please install {Matrix} first!")
    }
  }
  se <- suppressWarnings(ifelse(diag(vcov) > 0, sqrt(diag(vcov)), 0))
  cf <- cbind(cf, se, cf/se, 2 * pnorm(-abs(cf/se)))
  ## set to NA se, z and p-value for parameters with fixed values
  colnames(cf) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  cf <- list(parest = cf[seq.int(length.out = object$par.df), , drop = FALSE])
  rownames(cf$parest) <- names(object$MLE)[seq.int(length.out = object$par.df)]
  object$coefficients <- cf

  ## delete some slots
  object$hessian <- object$par.df <- object$N <- object$D <- object$optlb <- NULL
  object$MLE <- object$optub <- object$args <- NULL
  class(object) <- "summary.mlegc"
  object
}

print.summary.mlegc <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
  fitted <- x
  rm(x)
  cat("\nCall:", deparse(fitted$call, width.cutoff = floor(getOption("width")*0.6)), "", sep = "\n")
      cat("\nCoefficients of the model:\n")
      printCoefmat(fitted$coefficients$parest, digits = digits, signif.legend = TRUE)
     cat("\nlog likelihood = ", formatC(fitted$log.lik, digits = max(5L, digits + 1L)),
        ",  AIC = ", format(fitted$AIC, digits = max(4L, digits + 1L)),
        ",  BIC = ", format(fitted$BIC, digits = max(4L, digits + 1L)),
        ",  AICc = ", format(fitted$AICc, digits = max(4L, digits + 1L)),
         "\n", sep = "")
  invisible(fitted)
}


#### S3 Method on prediction results --------------------------------------------------------------------
#### ----------------------------------------------------------------------------------------------------
plot.predgc <- function(x, plottype = "2D", xlab = "xloc", ylab = "yloc",
                        xlim = NULL, ylim = NULL, pch = 20, textcex = 0.6, plotcex = 1,
                        angle = 60, col = c(2, 4), col.regions = gray(90:0/100),...)
{
# 2D; Predicted Counts; Predicted Means; Predicted Variance; 3D

  X <- x
  rm(x)
  locall <- rbind(X$obs.locs, X$pred.locs)
  Data <- c(X$obs.y, X$predCount)
  r1 <- apply(locall, 2, range)
  if(is.null(xlim)) xlim <- r1[,1]
  if(is.null(ylim)) ylim <- r1[,2]

  if (requireNamespace("latticeExtra", quietly = TRUE)) {
 ####  Plot1

  if(plottype == "2D"){
  print(lattice::levelplot(Data ~ locall[,1] + locall[,2], col.regions = col.regions,
                           xlab = xlab, ylab = ylab, cex = plotcex, main = NULL,
                           panel = function(...) {
                             lattice::panel.levelplot(...)
                             lattice::panel.text(X$obs.locs[,1], X$obs.locs[,2],
                                        X$obs.y, col = col[1], cex = textcex)
                             lattice::panel.text(X$pred.locs[,1], X$pred.locs[,2],
                                        X$predCount, col = col[2], cex = textcex)
                           }))
  }

  #### Plot2 Predicted Counts
  if(plottype == "Predicted Counts"){
    print(lattice::levelplot(X$predCount ~ X$pred.locs[,1] + X$pred.locs[,2],
                               panel = latticeExtra::panel.levelplot.points,
                               col.regions = col.regions, xlab = xlab, ylab = ylab, cex = plotcex,
                               main = NULL))
    }


  #### Plot3 Predicted Mean

  if(plottype == "Predicted Mean"){
  print(lattice::levelplot(X$predValue ~ X$pred.locs[,1] + X$pred.locs[,2],
                           panel = latticeExtra::panel.levelplot.points,
                           col.regions = col.regions, xlab = xlab, ylab = ylab, cex = plotcex,
                           main = NULL))
  }

  #### Plot4 Predicted Variance
  if(plottype == "Predicted Variance"){
  print(lattice::levelplot(X$predVar ~ X$pred.locs[,1] + X$pred.locs[,2],
                           panel = latticeExtra::panel.levelplot.points,
                           col.regions = col.regions, xlab = xlab, ylab = ylab, cex = plotcex,
                           main = NULL))
  }
  }else{
    stop("Unable to Provide level Plot. Please install {latticeExtra} first!")
  }

  ####
  if(plottype == "3D"){
  if (requireNamespace("scatterplot3d", quietly = TRUE)) {
  #### Plot 5

  K1 <- length(X$obs.y)
  K2 <- nrow(X$pred.locs)
  scatterplot3d::scatterplot3d(x = locall[,1], y = locall[,2], z = Data,
                               color = c(rep(col[1], K1), rep(col[2], K2)),
                               col.grid = "lightblue",  type='h', xlim = xlim, ylim = ylim,
                               col.axis = "blue", xlab = xlab, ylab = ylab, angle = angle,  pch = pch,
                               scale.y=1.2, y.margin.add = 0.2)
  }else{
    stop("Unable to Provide Three Dimensional Plot.
         Please install {scatterplot3d} first!")
  }
 }
}




summary.predgc <- function(object,...)
{
  if(!is.null(object$ConfidenceLevel)){
    answer1 <- as.data.frame(cbind(object$pred.locs,
                                   predMean = object$predValue,
                                   predCount = object$predCount,
                                   predVar = object$predVar))
    answer2 <- data.frame(ci1 = object$predInterval.EqualTail[,1],
                          ci2 = object$predInterval.EqualTail[,2],
                          ci3 = object$predInterval.Shortest[,1],
                          ci4 = object$predInterval.Shortest[,2])

    names(answer2) <- c(paste("Lower",object$ConfidenceLevel*100,"EqualTail"),
                        paste("Upper",object$ConfidenceLevel*100,"EqualTail"),
                        paste("Lower",object$ConfidenceLevel*100,"Shortest"),
                        paste("Upper",object$ConfidenceLevel*100,"Shortest"))
    answer <- cbind(answer1, answer2)
  } else {
    answer <- as.data.frame(cbind(object$pred.locs,
                                  predMean = object$predValue,
                                  predCount = object$predCount,
                                  predVar = object$predVar))
  }
  return(answer)
}

