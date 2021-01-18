draw.ddplot <- function(ddalpha, depth.space, cardinalities,
                        main = "DD plot", xlab = "C1", ylab = "C2", xlim, ylim, 
                        classes = c(1,2), colors = c("red", "blue", "green"), drawsep = T){
  
  if(!missing(ddalpha)){
    
    col = c()
    points = NULL
    for (c in classes){
    points = rbind(points, ddalpha$patterns[[c]]$depths[,classes])
    col = c(col, rep(colors[c], ddalpha$patterns[[c]]$cardinality))
    }
    
    if (class(xlab) != "expression" && class(ylab) != "expression" &&
        xlab == "C1" && ylab == "C2"){
      xlab = ddalpha$patterns[[1]]$name
      ylab = ddalpha$patterns[[2]]$name
    }
    
    if (missing(xlim)) xlim = c(0, max(points[,]))
    if (missing(ylim)) ylim = c(0, max(points[,]))
    
    plot(points, col = col, main = main, xlab = xlab, ylab = ylab, asp = T, xlim = xlim, ylim = ylim)
    
    if(drawsep && ddalpha$methodSeparator %in% c("alpha", "polynomial"))
    {
      gx <- seq(-0.1, 1.2*max(points[,1]), length=100)
      gy <- seq(0, 1.2*max(points[,2]), length=100)
      y <- as.matrix(expand.grid(x = gx, y = gy))   
      
      if (ddalpha$methodSeparator == "alpha") {
      ray = ddalpha$classifiers[[1]]$hyperplane[-1]
      
      funcs = list(function(x) x[1], function(x) x[2], function(x) x[1]^2, function(x) x[1]*x[2], function(x) x[2]^2, function(x) x[1]^3, function(x) x[1]^2*x[2], function(x) x[1]*x[2]^2, function(x) x[2]^3)
      
      depthcontours = apply(y, 1, function(xx) {
        res = 0
        for(i in 1:ddalpha$classifiers[[1]]$dimProperties)(res = res+funcs[[i]](xx)*ray[i])
        res
        })
      } else if (ddalpha$methodSeparator == "polynomial"){
        if (ddalpha$classifiers[[1]]$axis == 0){
          xAxis <- ddalpha$classifiers[[1]]$index1
          yAxis <- ddalpha$classifiers[[1]]$index2
        }else{
          xAxis <- ddalpha$classifiers[[1]]$index2
          yAxis <- ddalpha$classifiers[[1]]$index1
        }
        
       depthcontours = apply(y, 1, function(xx) {
          res = 0
        for(j in 1:ddalpha$classifiers[[1]]$degree){res <- res + ddalpha$classifiers[[1]]$polynomial[j]*xx[xAxis]^j}
          res = res-xx[yAxis]
        })
      }
      contour(gx, gy, matrix(depthcontours, nrow=length(gx), ncol=length(gy)), add=TRUE, levels=0, drawlabels=FALSE, col = "black")
    }
    
  } else if(!missing(depth.space)){
    
    col = c()
    for (c in classes){
      col = c(col, rep(colors[c], cardinalities[c]))
    }
    
    if (missing(xlim)) xlim = c(0, max(depth.space[,classes]))
    if (missing(ylim)) ylim = c(0, max(depth.space[,classes]))
    
    plot(depth.space[,classes], col = col, main = main, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim)
    
  } else stop("Both 'ddalpha' and 'depth.space' area missing")
  
}

plot.ddalpha <- function(x, type = c("ddplot", "depth.contours"), ...){
  type = match.arg(type)
  if(type == "ddplot") 
    draw.ddplot(x, ...)
  if(type == "depth.contours")
    depth.contours.ddalpha(x, ...)
}

plot.ddalphaf <- function(x, type = c("functional.data", "ddplot", "depth.contours"), ...){
  type = match.arg(type)
  if(type == "functional.data")
    plot.functional(list(dataf = x$dataf, labels = lapply(x$data[,ncol(x$data)], function(o){x$labels[[o]]})), ...)
  
  if(class(x$classifier)!="ddalpha")
    stop(type, " is available only for the ddalpha classifier")
  if(type == "ddplot") 
    draw.ddplot(x$classifier, ...)
  if(type == "depth.contours")
    depth.contours.ddalpha(x$classifier, ...)
}
  
  