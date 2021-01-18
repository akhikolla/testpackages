plot.PCMRS <- function(x, col.PCMRS = "red", col.PCM = "black", columns = NULL, ...){

  delta <- x$delta
  delta.PCM <- x$delta.PCM
  
  p <- nrow(delta)
  q <- ncol(delta)
  
  if(is.null(columns)){
    cols <- floor(sqrt(p))
  }else{
    cols <- columns
  }
  rows <- ceiling((p)/cols)
  
limits <- range(c(delta.PCM, delta))

  layout(matrix(1:(rows*cols),nrow=rows,byrow=TRUE))
  
  for(i in 1:nrow(delta)){
    plot(1:q, delta.PCM[i,],pch=19, ylab="",ylim=limits,xlab="",cex=2,xaxt="n",type="b", 
         main =rownames(delta)[i], col = col.PCM, ...)
    axis(1,1:q)
    lines(1:q,delta[i,],col=col.PCMRS,lwd=2,type="b")
  }
  
layout(1)

  invisible(x)
}