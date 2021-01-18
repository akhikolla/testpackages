

plot.qbnmix <- function(x, ...)
  {
    for (i in 1:length(x$mu))
      {
        plot.zoo(x$mu[[i]],type="l",xlab="",ylab="",screens="single",
                 main="",col=1:ncol(x$mu[[i]]),
                 lty=1:ncol(x$mu[[i]]))
        mtext(bquote(paste(mu,"_", .(i))),line=0.25)
      }

    for (i in 1:ncol(x$w))
      {
        plot(x$w[,i],type="l",xlab="",ylab="",main=paste("w_",i,sep=""))
      }

    for (i in 1:ncol(x$w))
      {
        plot(x$alpha[,i],type="l",xlab="",ylab="",main="")
        mtext(bquote(paste(alpha,"_", .(i))),line=0.25)
      }
  }

