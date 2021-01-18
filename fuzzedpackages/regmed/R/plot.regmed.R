plot.regmed <- function(x, cex=.9, lwd=1.5, lty=1, ...){
    
  ## input:
  ## x is output from regmed.grid.best (i.e, a single fitted model) or 
  ## output from regmed.fit
  if(class(x) != "regmed") {
      stop("input not regmed class")
  }
       
  x$exposure.name <- ifelse(x$call[[2]] == "x", "exposure", x$call[[2]])
  x$outcome.name <- ifelse(x$call[[4]] == "y",  "outcome", x$call[[4]])

  ## order alpha, beta by decreasing values of abs(alpha*beta)

  alpha <- x$alpha
  beta <- x$beta
  med.labels <- rownames(x$alpha)
  aprd <- abs(alpha*beta)
  ord <- rev(order(aprd))
  alpha <- alpha[ord]
  beta <- beta[ord]
  med.labels <- med.labels[ord]
    
  n.pts <- length(alpha)
    
  ## x dimensions, -3..3, with exposure at -2, mediators at -1..1, outcome at 2.
  ## y dimentions: -1 to n-mediators + 1, with direct at 0, mediators at 1..nmediators (npts)
    
  plot(c(-2.5, 2.5), c(-.5,n.pts + .5),type="n", xaxt = 'n', yaxt = 'n', bty = 'n', ylab = '', xlab = '', ...)
  midy <- (n.pts + 1)/2
  midx <- 0

  x0 <- rep(-2, n.pts)
  y0 <- rep(0, n.pts)
  x1 <- rep(-1, n.pts)
  y1 <- (n.pts:1)
  
#  if(n.pts == 1){
#      x1 <- -1
#      y1 <- 1
#  }
    
  col.alpha <- ifelse(alpha < 0, "blue", ifelse(alpha>0, "red", "gray"))
  arrows(x0, y0, x1, y1, length=.15, lty=lty, col=col.alpha, lwd=lwd, ...)
  
  text(-2, 0, x$exposure.name, cex=cex, pos=1, ...)
  text(x1+1, y1 + .2, med.labels, cex=cex*.9, pos=1, ...)
  
  x0 <- rep(x1+2, n.pts)
  y0 <- (n.pts:1)
  x1 <- rep(x1+3, n.pts)
  y1 <- rep(0, n.pts)
  col.beta <- ifelse(beta < 0, "blue", ifelse(beta>0, "red", "gray"))
  
#  if(n.pts == 1){
#    x0 <- 1
#    y0 <- 1
#    x1 <- 2
#    y1 <- 0
#  }
  
  arrows(x0, y0, x1, y1, length=.15, lty=lty, col=col.beta, lwd=lwd, ...)
  
  delta.col <- ifelse(x$delta < 0, "blue", ifelse(x$delta>0, "red", "gray")) 
  arrows(-2,0, 2, 0, length=.15, lty=lty, lwd=lwd, col=delta.col, ...)
    
  text(2, 0, x$outcome.name, cex=cex, pos=1, ...)
    
  invisible()
}
