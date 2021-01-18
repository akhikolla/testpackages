#' Annotate a threshold selection ggplot
#'
#' Annotate a threshold selection ggplot with the number
#' of exceedances of various thresholds.
#'
#' @param p An object produced by ggplot
#' @param x Horizontal axis data containing the full range.
#' @param y Verticle axis data containing the full range.
#' @param data The actual data being considered for GPD modelling.
#' @param textsize The size of the text in the annotations.
#' @export addExcesses
addExcesses <- function(p, x, y, data, textsize){
  x1 <- axisTicks(range(x), log=FALSE)
  yr <- range(y)
  delta <- abs(diff(yr)) * .1
  y1 <- rep(yr[2] + delta, length(x1))
  txt <- sapply(x1, function(u) sum(data > u))
  tx=data.frame(ex="Excesses:", x=min(x), y=y1[1] + delta)
  df <- data.frame(x=x1, y=y1, txt=txt)
  p <- p + geom_text(data=tx, aes(x,y,label=ex), size=textsize, hjust=0)
  p + geom_text(data=df, aes(x, y, label=txt), size=textsize)
}

#' @method ggplot mrl
#' @export
ggplot.mrl <- function(data, mapping, xlab = "Threshold", ylab = "Mean excess", main=NULL,
                       fill="orange", col="blue",
                       rug=TRUE, addNexcesses=TRUE, textsize=4, ..., environment){
  x <- data
  data <- x$data
  x <- x$mrl
  d <- data.frame(th = x[, "threshold"],
                  mrl = x[, "MRL"],
                  xl = x[, "lo"],
                  xu = x[, "hi"])

  k <- !is.na(d$xl)
  poly <- data.frame(x=c(d$th, rev(d$th)), y=c(d$xl, rev(d$xu)))
  poly <- poly[c(k, rev(k)), ]

  d2 <- data.frame(x = data,y=rep(0,length(data)))

  p <- ggplot(poly, aes(x, y)) +
    geom_polygon(fill=fill, alpha=.5) +
    geom_line(data=d, aes(th, mrl), color=col) +
    scale_x_continuous(xlab) +
    scale_y_continuous(ylab) +
    ggtitle(main)

  if (rug){
    p <- p + geom_rug(data=d2,mapping=aes(x,y),sides="b")
  }

  if (addNexcesses){
    p <- addExcesses(p, poly$x, poly$y, data=data, textsize=textsize)
  }

  p
}

#' @method ggplot gpdRangeFit
#' @export
ggplot.gpdRangeFit <- function(data, mapping, xlab = "Threshold", ylab = NULL, main = NULL,
                               fill="orange", col="blue",
                               addNexcesses = TRUE, textsize=4, ..., environment){
  if (missing(ylab)) {
    ylab <- c(expression(hat(sigma)[m]), expression(hat(xi)))
  }  else if (length(ylab) != 2) {
    stop("length of ylab should be 2")
  }
  if (!missing(main) && length(main) != 2) {
    stop("length of main should be 2")
  }

  x <- data
  data <- data$data

  p <- vector("list", 2)

  for (i in 1:2) {
    #        yl <- range(x$hi[, i], x$lo[, i])

    d <- data.frame(th=x$th, par=x$par[, i])
    poly <- data.frame(x=c(x$th, rev(x$th)), y=c(x$lo[, i], rev(x$hi[, i])))

    p[[i]] <- ggplot(poly, aes(x, y)) +
      geom_polygon(fill=fill, alpha=.5) +
      geom_line(data=d, aes(th, par), color=col) +
      scale_x_continuous(xlab) +
      scale_y_continuous(ylab[i]) +
      theme(axis.title.y=element_text(angle=0)) +
      if (!missing(main)) ggtitle(main[i])

    if (addNexcesses)
      p[[i]] <- addExcesses(p[[i]], poly$x, poly$y, data=data, textsize=textsize)
  } # Close for
  p
}
