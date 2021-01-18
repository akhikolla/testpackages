
#' Plot Method for Framed Data Clustering
#'
#' The \code{plot} method for framed data clustering result object.
#' It visualizes clusters on the input data that are within a best
#' frame.
#'
#' @import graphics
#' @import stats
#'
#' @param x an object of class \code{FramedClust} as returned by \code{FramedClust}
#' @param xlab a character string. The x-axis label for the plot.
#'  Default is NULL.
#' @param ylab a character string. The y-axis label for the plot.
#'  Default is NULL.
#' @param main a character string. The title for the plot.
#'  Default is NULL.
#' @param sub a character string. The subtitle for the plot.
#'  Default is NULL.
#' @param col.clusters a vector of colors, defined either by integers
#'  or by color names. If the length is shorter than the number
#'   of clusters, the colors will be reused. By default
#'  the blue, red3, green3, orange, purple, brown colors are used
#'  in the plot.
#' @param ... other arguments associated with the plot function
#'
#'
#' @return An object of class \code{"FramedClust"},
#' identical to the input \code{x}
#'
#'@examples
#' N <- 100
#' X <- rchisq(N, 5)
#' K <- 3
#' frame.size <- 40
#'
#' result <- FramedClust(X, K, frame.size)
#'
#' plot(result)


#' @export
plot.FramedClust <- function(x,
                             xlab = NULL,
                             ylab = NULL,
                             main = NULL,
                             sub = NULL,
                             col.clusters =  c(
                               "blue", "red3", "green3", "orange", "purple", "brown"
                             ),
                             ...)
{
  ck <- x

  if (exists(ck$X_name, mode = "numeric")) {
    X <- get(ck$X_name, mode = "numeric")
  } else {
    X <- eval(parse(text = ck$X_name))
  }

  color <- col.clusters

  #if (is.null(col.clusters))
  #{
  #  color = c("#009270",
  #            "#DC143C",
  #            "#0000CD",
  #            "#03faf2",
  #            "#c902c6",
  #            "#FA6A03")
  # } else{
  #  color = col.clusters
  # }

  I <- order(X)

  X <- X[I]

  cluster <- ck$cluster[I]

  unique.clusters <- unique(na.omit(cluster))

  ID <- which(cluster == unique.clusters[1])[1]

  if(ID != 1) ID <- ID - 1


  frame.size <- sum(!is.na(cluster))
  K <- max(cluster, na.rm = TRUE)

  if(is.null(main)) main <- "Optimal framed clustering"
  if(is.null(sub)) sub <- paste0("frame size = ", frame.size, " points, K = ", K)
  if(is.null(xlab)) xlab <- paste0(as.character(ck$X_name))
  if(is.null(ylab)) ylab <- ""
  plot(x = c(min(X), max(X)), y=rep(2,2),
       xlim = c(min(X), max(X)), ylim = c(0,3), type="n",
       xlab = xlab,
       ylab = ylab, yaxt = "n", sub = sub, main=main)

  segments(x0 = X, y0 = 0.1, x1 = X, y1 = 1.9, col = "grey")

  if(ID >= 1)
  {
    xleft = (X[( ID  )] + X[( ID + 1 )])/2
  }
  else
  {
    xleft =  X[( ID + 1 )]
  }


  if((ID + sum(ck$size) ) < length(X))
  {
    xright = (X[(ID + sum(ck$size) )] + X[(ID + sum(ck$size) + 1 )])/2
  }
  else
  {
    xright = X[(ID + sum(ck$size) )]
  }



  rect(xleft = xleft, ybottom = 0, xright = xright,
       ytop = 2, col = NULL, border = "black", lty="dashed")



  for(i in ( ID + 1 ):( ID + sum(ck$size) ) )
  {

    segments(x0 = X[i], y0 = 0.1, x1 = X[i], y1 = 1.9, col = color[cluster[i] %% length(color) + 1])

  }


  o <- order(unique.clusters)

  cluster_name <- paste0("Cluster ", unique.clusters)[o]

  cluster_color <- color[unique.clusters %% length(color) + 1][o]

  legend("topright",
         legend=cluster_name,
         col=cluster_color, lty=1, cex=0.8)



}
