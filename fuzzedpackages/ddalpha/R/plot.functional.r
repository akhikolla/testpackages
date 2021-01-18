plot.functional <- function(x, main = "Functional data", xlab = "args", ylab = "vals", colors = c("red", "blue", "green", "black", "orange", "pink"), ...) {
  
  if(main == "Functional data" && !is.null(x$name))
    main = x$name
  if(xlab == "args" && !is.null(x$args))
    xlab = x$args
  if(ylab == "vals" && !is.null(x$vals))
    ylab = x$vals
  
  ylims = matrix(unlist(lapply(x$dataf, function(e) (range(e$vals)))), ncol = 2, byrow = TRUE)
  plot(0, type="n", xlim=range(x$dataf[[1]]$args), ylim=c(min(ylims[,1]), max(ylims[,2])), 
       xlab=xlab, ylab=ylab, 
       main = main, ...)
  grid()
  
  if (!is.null(x$labels))
    labs = sort(unlist(unique(unlist(x$labels)))) # second unlist removes factors, other
  else 
    labs = NULL
  
  for (i in 1:length(x$dataf)){
    if (!is.null(labs))
      ind = match(x$labels[[i]],labs)
    else
      ind = 1
    lineColor <- colors[ind]
    
    lines(x$dataf[[i]]$args, x$dataf[[i]]$vals, col=lineColor)
  }
}

lines.functional <- function(x, colors = c("red", "blue", "green", "black", "orange", "pink"), ...) {
  
  if (!is.null(x$labels))
    labs = sort(unlist(unique(x$labels)))
  else 
    labs = NULL
  
  for (i in 1:length(x$dataf)){
    if (!is.null(labs))
      ind = match(x$labels[[i]],labs)
    else
      ind = 1
    lineColor <- colors[ind]
    
    lines(x$dataf[[i]]$args, x$dataf[[i]]$vals, col=lineColor, ...)
  }
}

points.functional <- function(x, colors = c("red", "blue", "green", "black", "orange", "pink"), ...) {
  
  if (!is.null(x$labels))
    labs = sort(unlist(unique(x$labels)))
  else 
    labs = NULL
  
  for (i in 1:length(x$dataf)){
    if (!is.null(labs))
      ind = match(x$labels[[i]],labs)
    else
      ind = 1
    pointColor <- colors[ind]
    
    points(x$dataf[[i]]$args, x$dataf[[i]]$vals, col=pointColor, ...)
  }
}