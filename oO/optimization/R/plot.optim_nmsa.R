
##--------------------##
#### Visualization  ####
##--------------------##

plot.optim_nmsa <- function (x, type = "convergence", lower = NA, upper = NA, ...) {

  par_save <- par(mar=c(3.5, 3.5, 1, 1) +0.1) # set graphical parametes
  if (is.null(x$trace)) {
    stop ("No trace matrix found. Check if trace was activated.")
  }

  #--------------------#
  ## Convergence plot ##
  #--------------------#

  if (type == "convergence") {
    plot(x[[3]][, 1], x[[3]][, 2],  type = "n", lwd = 2, axes = FALSE, # opens plot window only with title and axes labeling
         main = "", xlab = "", ylab = "", ...)
    a1 <- axis(1,  cex.axis = 0.9, lwd = 0, lwd.ticks = 1) # adds axes scaling manually according to the trace values
    a2 <- axis(2, cex.axis = 0.9, las = 1, lwd = 0, lwd.ticks = 1)
    mtext("function value", side = 2, line = 2.5, las= 3) ; mtext("iteration", side = 1, line = 2, las = 1)

    abline(h = a2, lty = 2, lwd = 1, col =  grey.colors(1, 0.8))
    abline(v = a1, lty = 2, lwd = 1, col =  grey.colors(1, 0.8)) # adds a grid to the plot
    abline(h = x[[2]], lwd = 2, lty = 6, col = grey.colors(1, 0.5)) # adds the last value (optimum) as a line
    lines(x[[3]][ , 2] ~ x[[3]][ , 1], lwd = 3) # adds the convergence line
    box()
  }

  #----------------#
  ## Contour plot ##
  #----------------#

  if (type == "contour") {
    if (length(x$par)!=2) {stop ("Contour plot only works with two dimensional functions.")}
    data <- x$trace
    fun <- x$fun
    if (any(is.na(c(lower, upper)))) {
      if (any(is.na(c(x$lower, x$upper)))) {
        lower <- apply(data[,c(3, 4)], 2, min)
        upper <- apply(data[,c(3, 4)], 2, max)
      } else {
        upper <- x$upper
        lower <- x$lower
      }
    }# If there are no limits, they are set to the min (max respectively) of the trace matrix.

    # If trace is to long, it becomes very confusing. It is therefore shortened to a length that is nearest to 100 and possible.
    if (dim(data)[1] > 100) {
      cat(length(seq(1, dim(x$trace)[1], round(dim(x$trace)[1]/100))), " iterations are displayed.")
      data <- data[seq(1, dim(data)[1], round(dim(data)[1]/100)),]
    } else {cat(dim(x$trace)[1], "iterations are displayed")}


    # Calculating the step size to determine the pixel density of the contour plot. 250 x 250 pixel is by far enough for high quality color graphics.
    sz_x1 <- (upper[1] - lower[1]) / 250
    sz_x2 <- (upper[2] - lower[2]) / 250
    x1 <- seq(lower[1], upper[1], sz_x1)
    x2 <- seq(lower[2], upper[2], sz_x2)

    # Create an NA matrix y with dimension x1 x x2.
    y <- matrix(rep(NA, length(x1) * length(x2)), ncol = length(x2), byrow = T)

    # Calculate the objective function value of the variables and save them in the matrix y.
    for (i in seq(x1)) {
      for (j in seq(x2)) {
        y[i, j] <- fun(c(x1[i], x2[j]))
      }
    }

    # Image creates the contour plot.
    if (requireNamespace('colorspace', quietly = TRUE)) {
      image(x1, x2 , y, xlim = c(lower[1], upper[1]) ,ylim = c(lower[2], upper[2]), col = rev(colorspace::heat_hcl(500)), axes = F, ann = F, ...)
    } else {
      image(x1, x2 , y, xlim = c(lower[1], upper[1]) ,ylim = c(lower[2], upper[2]), col = rev(heat.colors(500)), axes = F, ann = F, ...)
    }


    # Add axes and frame.
    axis(1, lwd = 0, lwd.ticks = 1, cex.axis = 0.9)
    axis(2, lwd = 0, lwd.ticks = 1, las = 2, cex.axis = 0.9)
    box()
    # Add contour lines and axes labels.
    contour(x1, x2, y, xlim = c(lower[1], upper[1]), ylim = c(lower[2], upper[2]), add= T, nlevels = 20)
    mtext("x1", side = 1, line = 2) ; mtext("x2", side = 2, line = 2.5, las = 3)
    # Add the points of the trace matrix.
    bg_col_points <- gray(seq(0, 1, length.out=dim(data)[1]))
    points(x = data[ ,3], y = data[ ,4], pch = 21, cex = 1.8, col = "black", bg = bg_col_points)
    # Add text to 10 points if there are more than 10 entrys in the trace matrix. Else add text to all points.
    if (dim(data)[1] > 10) {
      for (i in seq(0.1,1,0.1)) {
        if (i > 0.7) { # Choosing the color of the text.
          text(labels = data[trunc(length(data[, 1]) * i), 1], x = data[trunc(length(data[, 1]) * i), 3], y = data[trunc(length(data[, 1]) * i), 4], cex = 0.5, col = "black")
        } else {
          text(labels = data[trunc(length(data[, 1]) * i), 1], x = data[trunc(length(data[, 1]) * i), 3], y = data[trunc(length(data[, 1]) * i), 4], cex = 0.5, col = "white")
        }

      }
    } else {
      text(labels = data[, 1], x = data[, 3], y = data[, 4], cex = 0.5, col = c(rep("white", round(dim(data)[1] * 0.6)), rep("black", round(dim(data)[1] * 0.4))))
    }

    # Add the final objective function variables. The white frame improves the contrast when a trace point is behind the star.
    points(x = x$par[1], y = x$par[2], pch = 8, cex = 1.2, lwd = 3, col = "white")
    points(x = x$par[1], y = x$par[2], pch = 8, cex = 1.2, lwd = 2)
  }


  par(par_save) # Restore the graphical parameters

}
