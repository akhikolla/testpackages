# plot-CirClust.R

#' Plot Method for Circular Data Clustering
#'
#' The \code{plot} method for circular data clustering result
#' object of class \code{CirClust}.
#' It visualizes circular clusters on the input data.
#'
#' @importFrom plotrix draw.circle
#' @importFrom plotrix draw.radial.line
#' @importFrom plotrix arctext
#' @import graphics
#'
#' @param x an object of class as returned by \code{CirClust}
#' @param xlab a character string. The x-axis label for the plot.
#'  Default is no string.
#' @param ylab a character string. The y-axis label for the plot.
#'  Default is no string.
#' @param main a character string. The title for the plot.
#' @param sub a character string. The subtitle for the plot.
#' @param col.clusters a vector of colors, defined either by
#'  integers or by color names. If the length is shorter than
#'  the number of clusters, the colors will be reused. By default
#'  the blue, red3, green3, orange, purple, brown colors are used
#'  in the plot.
#' @param axes the axis will be ploted if set TRUE. Default is FALSE.
#' @param xlim range of the x axis in the plot.
#'    Default is from -1.75 to 1.75.
#' @param ylim range of the y axis in the plot.
#'    Default is from -1.75 to 1.75.
#' @param fill the color to fill inside the ring as the background
#' of data points.
#' @param border the color to draw cluster borders.
#' @param border.lty the line type to draw cluster borders.
#' @param ... other arguments associated with the plot function
#'
#'
#' @return A copy of the input object of class \code{CirClust}.
#'
#'@examples
#' opar <- par(mar=c(1,1,2,1))
#
#' # Example 1. Circular data clustering
#' n <- 100
#' Circumference <- 7
#' O <- runif(n, 0, Circumference)
#' result <- CirClust(O, K=3, Circumference=Circumference)
#' plot(result, fill="mintcream", main="Example 1. Circular clustering")
#'
#'
#' # Example 2. Circular data clustering
#' n <- 40
#' m <- 5
#' O <- c(rnorm(n,mean=5,sd=m), rnorm(n,mean=15,sd=m), rnorm(n,mean=26,sd=m))
#' K <- 3
#' Circumference <- 28
#'
#' result <- CirClust(O, K, Circumference, method = "FOCC")
#'
#' color <- c("royalblue", "green3", "firebrick") # c("#0000CD","#808080", "#DC143C")
#'
#' par(mar=c(1,1,2,1))
#'
#' plot(result, col.clusters = color, fill="floralwhite",
#'      main="Example 2. Circular clustering")
#'
#'
#' # Example 3. Periodic data clustering
#' n <- 100
#' period <- 5.2
#' O <- rnorm(n)
#' result <- CirClust(O, K=5, Circumference=period)
#' plot(result, fill="navy", border="gray", border.lty="dotted",
#'      main="Example 3. Periodic clustering")
#'
#' par(opar)

#' @export
plot.CirClust <- function(
  x,
  xlab = "",
  ylab = "",
  main = NULL,
  sub = "",
  col.clusters = c(
    "blue", "red3", "green3", "orange", "purple", "brown"
  ),
  axes = FALSE,
  xlim = c(-1.75, 1.75),
  ylim = c(-1.75, 1.75),
  fill = "floralwhite",
  border = "gray36",
  border.lty = "dotted",
  ...)
{
  ck <- x

  color <- col.clusters

  # if (is.null(col.clusters))
  # {
  #   color = c("#009270",
  #             "#DC143C",
  #             "#0000CD",
  #             "#000000",
  #             "#c902c6",
  #             "#FA6A03")
  # } else{
  #   color =  col.clusters
  # }

  if (exists(ck$O_name, mode = "numeric")) {
    O <- get(ck$O_name, mode = "numeric")
  } else {
    O <- eval(parse(text = ck$O_name))
  }

  if(is.null(main)) main <- paste0(as.character(ck$O_name))

  plot(
    c(0, 1),
    c(0, 0),
    main = main,
    xlim = xlim,
    ylim = ylim,
    axes = axes,
    xlab = xlab,
    ylab = ylab,
    sub = sub,
    ...
  )
  draw.circle(0, 0, c(1.15, 0.75), col = c(fill, "#FFFFFF"))



  draw.radial.line(0, 0.75, c(0, 0), col = c("#000000"))

  arrows(0.35, 0.15, 0.35, 0.55, length = 0.25, angle = 30)

  text(0.65, -0.075, labels = as.character(0))

  Circumference = ck$Circumference


  U <- unique(ck$cluster)

  count <- 1

  clust <- ck$cluster

  for(i in seq_along(U))
  {
    clust[which(ck$cluster == U[i])] <- count

    count <- count + 1
  }


  for (i in seq_along(ck$cluster))
  {
    angle = O[i] / Circumference

    draw.radial.line(0.78,
                     1.12,
                     c(0, 0),
                     angle = (angle * 2 * pi),
                     col = color[clust[i] %% length(color) + 1])
  }




  for (i in seq_along(ck$Border.mid))
  {
    angle = ck$Border.mid[i] / Circumference

    draw.radial.line(
      0.65,
      1.25,
      c(0, 0),
      angle = (angle * 2 * pi),
      col = border,
      lty = border.lty,
      lwd = 2
    )

  }


  for (i in seq_along(ck$centers))
  {
    angle = ck$centers[i] / Circumference

    arctext(
      paste0("C", i),
      center = c(0, 0),
      radius = 1.25,
      middle = (angle * 2 * pi)
    )

  }

  invisible(ck)
}
