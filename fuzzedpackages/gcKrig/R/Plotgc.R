#### Partially inherited from geoCount... For plot a geocount data alone. 
plotgc <- function(data = NULL, locs = NULL, bdry=NULL, col=2, pch=1, 
                   textcex = 1, col.regions = gray(90:0/100), 
                   size=c(0.3, 2.7), ...){
  if(!is.null(bdry)){
    if(!is.list(bdry)) bdry <- list(bdry)
    n <- length(bdry)
    tt <- sapply(bdry, function(t) apply(t, 2, range))
    tt1 <- apply(tt, 1, range)
    xrange <- range(tt1[1:4])
    yrange <- range(tt1[5:8])
    plot(xrange, yrange, type="n", ...)
    for(i in 1:n) lines(bdry[[i]])
  } else{
    plot(locs[,1], locs[,2], type="n", ...)
  }
  if(is.null(data)) data <- 0
  points(locs[,1], locs[,2], col=col, pch=pch,
         cex = size[1]+size[2]*(data-min(data))/(max(data)-min(data)))
  
  if(is.null(bdry)){
  print(lattice::levelplot(data ~  locs[,1] + locs[,2], col.regions = col.regions,
                          # xlab = xlab, ylab = ylab, cex = plotcex,
                           panel = function(x,y,z,...) {
                             lattice::panel.levelplot(x,y,z,...)
                             lattice::panel.text(locs[,1], locs[,2], data, 
                                                 col = col, cex = textcex)
                           }))
  }
}
