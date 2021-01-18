#######################    plotJoint     ########################
#' Plot Joint CDF after Nonparametric Analysis
#'
#' @description
#' This function plots the joint cdf for a \verb{bivrecNP} object.
#'
#' @param object A \verb{bivrecNP} object.
#' @param type  Optional vector of strings to label Type I and Type II gap times. Default is c("Type I", "Type II").
#' @importFrom stats ftable
#' @export
#'

#@return A contour plot of joint cdf.

plotJoint <- function(object, type) {

  x = object

  if (!inherits(x, "bivrecNP")) stop("Must be a bivrecNP object.")
  if (missing(type)) {type = c("Type I", "Type II")}

  forplot <- x$joint_cdf

  #####OLD MAY RE-USE LATTER: Wald CI and plot
  # rgl::plot3d(forplot[,1], forplot[,2], forplot[,3], col = "black", xlab = "x",
  #        main = "Joint cdf", ylab ="y", zlab = expression(P(X^0 <= x, Y^0 <= y)),  expand = 1.1, cex.lab = 1.5)
  # for (i in 1:nrow(forplot)) {
  #   rgl::rgl.lines(forplot[i,1], forplot[i,2], as.matrix(forplot[i,5:6]), col="red")
  # }

  forplot <- forplot[1:3]
  colnames(forplot) <- c("X", "Y", "Cumm.Prob")
  myx <- as.factor(forplot$X)
  myy <- as.factor(forplot$Y)
  lx <- length(levels(myx))
  forplot2 <- as.matrix(ftable(forplot, row.vars = 1, col.vars = 2))
  rownames(forplot2) <- levels(myx)
  colnames(forplot2) <- levels(myy)
  for (i in 1:lx) {
    index <- which(forplot$X==as.numeric(rownames(forplot2)[i]))
    forplot2[i,] = forplot$Cumm.Prob[index]
  }

  graphics::filled.contour(x=as.numeric(levels(myx)), y= as.numeric(levels(myy)),
                           forplot2, color.palette = grDevices::heat.colors, cex.main=1.5,
                           xlab=type[1], ylab=type[2],
                           main = expression(P(X^0 <= x, Y^0 <= y)))

}

########################    plotMarg     ########################
#' Plot Marginal Survival after Nonparametric Analysis
#'
#' @description
#' This function plots the marginal survival for a \verb{bivrecNP} object.
#'
#' @param object An object of \verb{bivrecNP} class.
#' @param type Optional string to label the Type I gap time. Default is "Type I".
#' @export

#@return A plot of marginal survival vs. first gap time with confidence interval.

plotMarg <- function(object, type) {
  x <- object

  if (!inherits(x, "bivrecNP")) stop("Must be a bivrecNP object.")
  if (missing(type)) {type = "Type I"}

  xij <- x$xij
  forplot <- x$marginal_survival

  #formula <- bivrec.nonparam.result$formula

  ## variables <- all.vars(formula)
  mx <- 10 * (max(xij) %/% 10 + 1)
  forplot <- rbind(c(0, 1, 0, 1, 1), forplot, c(mx, 0, forplot[nrow(forplot),3], 0, 0))

  ##### Wald CI and plot
  index <- which(forplot$`Lower .99`<0)
  forplot[index, -1] <- forplot[index[1]-1, -1]

  plot(forplot[,1], forplot[,2], type = "l", xlab = type,
       ylab = "Marginal Survival", yaxp  = c(0, 1, 10),
       xaxp  = round(c(0, mx, 10), digits=1), main = expression(1 - P(X^0 <= x))
  )
  graphics::lines(forplot[,1], forplot[,4], lty = 2)
  graphics::lines(forplot[,1], forplot[,5], lty = 2)

}

########################    plotCond     ########################
#' Plot Conditional CDF after Nonparametric Analysis
#'
#' @description
#' This function plots the conditional cdf for a \verb{bivrecNP} object.
#'
#' @param object An object of \verb{bivrecNP} class where the analysis has specified \verb{conditional = TRUE}.
#' @param type Optional string to label the Type II gap time. Default is "Type II".
#'
#' @importFrom stats ftable
#' @export

#@return A plot of conditional cdf in the given interval.

plotCond <- function(object, type) {
  x=object

  cond <- x$conditional_cdf$conditional

  if (missing(type)) {type = "Type II"}

  plot(cond$Time, cond[,5], type="l", lty = 2, xlab = type,
       ylab = "Conditional Probability", xlim=c(0, round(max(x$conditional_cdf$ygrid), digits=1)),
       ylim=c(0, round(max(cond[,5]), digits=1)),
       main=substitute(paste("P(", Y^0 <= y, "|", X^0 %in% "[", gi1, ",", gi2, "])"),
                       list(gi1 = x$given.interval[1], gi2 = x$given.interval[2]))
  )
  graphics::lines(cond$Time, cond[,4], lty = 2)
  graphics::lines(cond$Time, cond[,2], lty = 1)
}

########################    plot.bivrecNP     ########################
#' Plot Results of Nonparametric Analysis
#'
#' @description
#' This function plots all the estimated functions (joint cdf, marginal survival and conditional cdf if \verb{conditional=TRUE} during analysis) from a \verb{bivrecNP} object in one step.
#'
#' @import graphics
#' @importFrom utils tail
#' @importFrom stats model.frame
#' @importFrom stats na.omit
#' @importFrom stats model.matrix
#'
#' @param x A \verb{bivrecNP} object.
#' @param y  Either empty or NULL.
#' @param type Optional vector of strings to label Type I and Type II gap times. Default is c("Type I", "Type II").
#' @param ... Additional arguments to be passed to graphical methods if needed.
#'
#' @export
#'

plot.bivrecNP <-function(x, y=NULL, type=NULL, ...){

  if (!inherits(x, "bivrecNP")) stop("Must be a bivrecNP object.")

  cond=x$conditional #boolean saying if conditional is in bivrecNP object

  if (cond==FALSE){
    par(mar=c(5,4,4,2)+0.1)
    if (is.null(type)) {plotJoint(x)} else {plotJoint(x, type)}
    par(mar=c(5,4,4,2)+0.1)
    if (is.null(type)) {plotMarg(x)} else {plotMarg(x, type[1])}
  } else {
    par(mar=c(5,4,4,2)+0.1)
    if (is.null(type)) {plotJoint(x)} else {plotJoint(x, type)}
    par(mar=c(5,4,4,2)+0.1, mfrow=c(1,2))
    if (is.null(type)) {plotMarg(x)} else {plotMarg(x, type[1])}
    if (is.null(type)) {plotCond(x)} else {plotCond(x, type[2])}
    par(mfrow=c(1, 1))
  }
}
