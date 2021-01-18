
#' Conduct goodness of fit diagnostics
#'
#'
#' @param object the object to evaluate
#' @param ... additional parameters
#'
#'
#' @details
#' see \code{\link{gofit.lolog}}
gofit <- function(object, ...) {
  UseMethod("gofit")
}


#' Goodness of Fit Diagnostics for a LOLOG fit
#'
#'
#' @param object the object to evaluate
#' @param formula A formula specifying the statistics on which to evaluate the fit
#' @param nsim The number of simulated statistics
#' @param ... additional parameters
#'
#' @examples
#' library(network)
#' data(ukFaculty)
#' 
#' # Delete vertices missing group
#' delete.vertices(ukFaculty, which(is.na(ukFaculty %v% "Group")))
#' 
#' # A dyad independent model
#' fitind <- lolog(ukFaculty ~ edges() + nodeMatch("GroupC") + nodeCov("GroupC"))
#' summary(fitind)
#'
#' # Check gof on degree distribution (bad!)
#' gind <- gofit(fitind, ukFaculty ~ degree(0:50))
#' gind
#' plot(gind)
#'
#' #check gof on esp distribution (bad!)
#' gind <- gofit(fitind, ukFaculty ~ esp(0:25))
#' gind
#' plot(gind)
#'
#' \dontrun{
#'
#' #include triangles and 2-stars (in and out)
#' fitdep <- lolog(ukFaculty ~ edges() + nodeMatch("GroupC") + nodeCov("GroupC") + 
#'                 triangles + star(2, direction="in") + star(2, direction="out"), nsamp=1500)
#' summary(fitdep)
#'
#' # Check gof on (in + out) degree distribution (good!)
#' gdep <- gofit(fitdep, ukFaculty ~ degree(0:50))
#' gdep
#' plot(gdep)
#'
#' #check gof on esp distribution (good!)
#' gdep <- gofit(fitdep, ukFaculty ~ esp(0:25))
#' gdep
#' plot(gdep)
#'
#' }
#'
#'
#' @method gofit lolog
gofit.lolog <- function(object, formula, nsim = 100, ...) {
  model <- createCppModel(formula)
  observedNetwork <- object$likelihoodModel$getModel()$getNetwork()
  model$setNetwork(observedNetwork)
  model$calculate()
  ostats <- model$statistics()
  ns <- length(ostats)
  stats <- matrix(0, nrow = nsim, ncol = ns)
  model$calculate
  for (i in 1:nsim) {
    net <- object$likelihoodModel$generateNetwork()$network
    model$setNetwork(net)
    model$calculate()
    stats[i, ] <- model$statistics()
  }
  
  mn <- apply(stats, 2, min)
  mx <- apply(stats, 2, max)
  pv <- apply(sweep(stats, 2, ostats), 2, function(a) {
    if (sd(a) < .Machine$double.eps)
      return(NA)
    2 * pnorm(-abs(mean(a) / sd(a)))
  })
  d <- data.frame(
    obs = ostats,
    min = mn,
    mean = colMeans(stats),
    max = mx,
    pvalue = pv
  )
  
  result <- list(ostats = ostats,
                 stats = stats,
                 summary = d)
  class(result) <- "gofit"
  result
}

#' prints a gofit object
#'
#'
#' @param x The object
#' @param ... passed to print.data.frame
#'
#'
#' @method print gofit
print.gofit <- function(x, ...) {
  print(x$summary, ...)
}

#' Plots a gofit object
#'
#'
#' @param x the gofit object
#' @param y unused
#' @param type type of plot, boxplot or lineplot
#' @param scaling type of scaling of the network statistics. If "std", network statistics are scaling by subtracting off the observed statistics and scaling by the standard deviation. If "sqrt", network statistics are plotted on the square root scale (The square root is the variance stabilizing transformation for a Poisson random variable). The default is "none", where by the network statistics are not scaled. 
#' @param lineAlpha The transparency of the simulated statistics lines
#' @param lineSize The width of the lines
#' @param ... passed to either boxplot or geom_line
#'
#'
#' @examples
#' library(network)
#' data(ukFaculty)
#' 
#' # Delete vertices missing group
#' delete.vertices(ukFaculty, which(is.na(ukFaculty %v% "Group")))
#' 
#' # A dyad independent model
#' fitind <- lolog(ukFaculty ~ edges() + nodeMatch("GroupC") + nodeCov("GroupC"))
#' summary(fitind)
#'
#' # Check gof on degree distribution (bad!)
#' gind <- gofit(fitind, ukFaculty ~ degree(0:50))
#' plot(gind)
#' plot(gind, type="box")
#' @method plot gofit
plot.gofit <- function(x,
                       y,
                       type = c("line", "box"),
                       scaling = c("none","std","sqrt"),
                       lineAlpha = .06,
                       lineSize = 1,
                       ...) {
  type <- match.arg(type)
  stats <- x$stats
  ostats <- x$ostats
  nms <- names(ostats)
  colnames(stats) <- nms
  ylab <- switch(scaling[1],
   sqrt="Statistic (on a square root scale)",
   std="Normalized Statistic",
   "Statistic"
  )
  
  if(ncol(stats) == 1 && type == "line"){
    message("Note: Cannot create line plot with only one statistic. Falling back to boxplot")
    type <- "box"
  }
  
  if (type == "box") {
    stats <- switch(scaling[1],
     sqrt={
      ostats <- sqrt(ostats)
      sqrt(stats)
     },
     std={
      stats <- sweep(stats, 2, ostats)
      ostats <- rep(0, length(ostats))
      apply(stats, 2, function(a) {
        if (sd(a) < .Machine$double.eps)
          return(rep(NA, length(a)))
        a / sd(a)
      })
     },
     {
       stats
     }
    )
    boxplot(stats, ylab = ylab, yaxt='n', ...)
    ylabs <- switch(scaling[1],
     sqrt=graphics::axTicks(2)^2,
     std=graphics::axTicks(2),
     graphics::axTicks(2)
    )
    graphics::axis(side=2, at=graphics::axTicks(2), labels=ylabs)
    points(ostats, col = "red", pch = 16)
    return(invisible(NULL))
  } else{
    Var2 <- value <- Var1 <- xx <- yy <- gg <- NULL #For R CMD check
    stats <- switch(scaling[1],
     sqrt={
      transf <- "sqrt"
      stats
     },
     std={
      transf <- "identity"
      stats <- sweep(stats, 2, ostats)
      ostats <- rep(0, length(ostats))
      apply(stats, 2, function(a) {
        if (sd(a) < .Machine$double.eps)
          return(rep(NA, length(a)))
        a / sd(a)
      })
     },
     {
      transf <- "identity"
      stats
     }
    )
    mstats <- reshape2::melt(stats)
    o <- data.frame(xx = nms, yy = ostats, gg = "observed")
    gg <- ggplot2::ggplot(data = mstats) +
      ggplot2::geom_line(
        ggplot2::aes(x = Var2, y = value, group = Var1),
        alpha = lineAlpha,
        size = lineSize,
        ...
      ) +
      ggplot2::geom_line(
        ggplot2::aes(x = xx, y = yy, group = gg),
        data = o,
        color = "red",
        size = lineSize,
        ...
      ) +
      ggplot2::scale_y_continuous(trans = transf)  +
      ggplot2::theme_bw() + ggplot2::ylab(ylab) + ggplot2::xlab("")  +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor.x = ggplot2::element_blank()
      )
    return(gg)
  }
}
