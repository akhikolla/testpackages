# plot_table.R
#
# Joe Song
# Created: October 17, 2018
#   Adapted from discrete-patterns.R April 3, 2017
#
# Updated: 2018-11-01
#   Added an argument in plot_table() to highlight row maxima
#   with a box.
#
# Updated: 2018-11-05
#   Determine the highlight box value by median value of the table.

plot_table <- function(
  table, xlab="Column", ylab="Row",
  col="green3",
  xaxt="n", yaxt="n", main=NULL,
  show.value=TRUE, value.cex=2,
  highlight=c("row.maxima", "none"),
  highlight.col=col,
  mgp=c(0.5,0,0), mar=c(2,2,3,1.5),
  ...)
{
  if(!is.matrix(table) && !is.data.frame(table)){
    stop("input x must be matrix or data frame\n")
  }

  min.col <- "white"
  box.size <- 0.85
  box.lwd <- 4 # 1.5 #3

  pal <- colorRampPalette(c(min.col, col), space = "Lab")

  if(is.null(main)) {
    main <- deparse(substitute(table))
  }

  table <- as.matrix(table)

  nrow <- nrow(table)
  ncol <- ncol(table)

  m <- t(table)[, rev(seq(nrow)), drop=FALSE]

  opar <- par(mgp=mgp, mar=mar)
  image(m, main=main, xlab=xlab, ylab=ylab,
        xaxt=xaxt, yaxt=yaxt,
        col=pal(ncol*nrow), ...
  )

  highlight <- match.arg(highlight)

  if(highlight == "row.maxima") {

    # med.value <- median(table)
    mid.value <- (min(table) + max(table)) / 2

    for(r in seq(nrow)) {
      p <- which.max(m[, r])

      if(ncol > 1) {
        cx <- (p - 1) / (ncol - 1)
        right <- cx + box.size * 0.5 / (ncol-1)
        left <- cx - box.size * 0.5 / (ncol-1)
      } else {
        cx <- 0
        right <- cx + box.size
        left <- cx - box.size
      }

      if(nrow > 1) {
        cy <- (r - 1) / (nrow - 1)
        top <- cy + box.size * 0.5 / (nrow-1)
        bottom <- cy - box.size * 0.5 / (nrow-1)
      } else {
        cy <- 0
        top <- cy + box.size
        bottom <- cy - box.size
      }

      if(1) {
        # border <- ifelse(m[p,r] >= mid.value, min.col, col)
        # border <- highlight.col
        # rect(left, bottom, right, top, lwd=box.lwd, border=border)
        rect(left, bottom, right, top, lwd=box.lwd * 2, border=min.col)
        rect(left, bottom, right, top, lwd=box.lwd, border=col)

      } else {
        polygon(c(left, cx, right, cx), c(cy, top, cy, bottom),
                border = min.col, lwd=7)
        polygon(c(left, cx, right, cx), c(cy, top, cy, bottom),
                border = highlight.col, lwd=4)
      }

    }
  }

  grid(nx=ncol, ny=nrow, col=col)

  if(show.value) {

    if(ncol == 1) {
      tx <- 0
    } else {
      tx <- rev(rep(seq(ncol), each=nrow)-1) / (ncol-1)
    }

    if(nrow == 1) {
      ty <- 0
    } else {
      ty <- (rep(seq(nrow), times=ncol)-1) / (nrow-1)
    }

    text(tx, ty, rev(as.vector(table)), cex=value.cex)
  }

  par(opar)
}
