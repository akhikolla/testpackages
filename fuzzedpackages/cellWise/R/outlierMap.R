
outlierMap <- function(res, title = "Robust PCA", col = "black", 
                       pch = 16, labelOut = TRUE, id = 3, 
                       xlim = NULL, ylim = NULL, cex = 1,
                       cex.main=1.2, cex.lab=NULL, cex.axis=NULL){
  
  
  if(is.null(res$sd)) { SD = res$SD } else { SD = res$sd }
  if(is.null(res$cutoff.sd)) { cutoffSD = res$cutoffSD } else {
    cutoffSD = res$cutoff.sd }
  if(is.null(res$od)) { OD = res$OD } else { OD = res$od }  
  if(is.null(res$cutoff.od)) { cutoffOD = res$cutoffOD } else {
    cutoffOD = res$cutoff.od }
  if(is.null(SD)) stop(" No score distances are given.")
  if(is.null(cutoffSD)) stop(
    " No cutoff for score distances is given.")
  if(is.null(OD)) stop(" No orthogonal distances are given.")
  if(is.null(cutoffOD)) stop(
    " No cutoff for orthogonal distances is given.")
  if (is.null(xlim)) {
    xlim <- c(0, max(SD) * 1.1)
  }
  if (is.null(ylim)) {
    ylim <- c(0, max(OD) * 1.1)
  }
  if (is.null(pch)) {
    pch <- ifelse(is.null(col), 1, 16)
  }
  if (is.null(col)) {
    col <- "black"
  }
  mycex = cex
  mycex.main = 1
  if(!is.null(cex.main)) { mycex.main = cex.main }
  mycex.lab = 1
  if(!is.null(cex.lab)) { mycex.lab = cex.lab }
  mycex.axis = 1
  if(!is.null(cex.axis)) { mycex.axis = cex.axis }
  if (is.list(col)) {
    for (i in 1:length(col)) {
      if (i == 1) {
        plot(SD[col[[i]]$index], OD[col[[i]]$index], 
             xlab = "", ylab = "", main = "", pch = pch, 
             col = col[[i]]$col, xlim = xlim, ylim = ylim, 
             cex=mycex, cex.axis = mycex.axis)
      }
      points(SD[col[[i]]$index], OD[col[[i]]$index], pch = pch, 
             col = col[[i]]$col, cex=mycex)
    }
  }
  else {
    plot(SD, OD, xlab = "", ylab = "", main = "", pch = pch, 
         col = col, xlim = xlim, ylim = ylim, cex=mycex, 
         cex.axis = mycex.axis)
  }
  title(main = title, line = 1, cex.main = mycex.main)
  title(ylab = "Orthogonal distance", line = 2.3, cex.lab = mycex.lab)
  title(xlab = "Score distance", line = 2.3, cex.lab = mycex.lab)  
  abline(v = cutoffSD)
  abline(h = cutoffOD)
  if (labelOut) { # had to add cellWise::: to make next line work:
    labelDD_cw(SD, OD, id.n.SD = id, id.n.OD = id)
  }
}


labelDD_cw <- function(x,y,id.n.SD=3,id.n.OD=3,off=0.02) 
{ # used in OutlierMap
  xrange <- graphics::par("usr")
  xrange <- xrange[2] - xrange[1]
  if (id.n.SD > 0 && id.n.OD > 0) {
    n <- length(x)
    indexSD <- sort(x, index.return = TRUE)$ix
    indexSD <- indexSD[(n - id.n.SD + 1):n]
    indexOD <- sort(y, index.return = TRUE)$ix
    indexOD <- indexOD[(n - id.n.OD + 1):n]
    lab <- indexOD
    if (is.character(names(y))) {
      lab <- names(y[indexOD])
    }
    graphics::text(x[indexOD] - off * xrange, y[indexOD], lab)
    lab <- indexSD
    if (is.character(names(x))) {
      lab <- names(x[indexSD])
    }
    graphics::text(x[indexSD] - off * xrange, y[indexSD], lab)
  }
}