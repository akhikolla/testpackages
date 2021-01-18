

plot.tvpreg <- function(x, ...)
  {
    opar <- par(no.readonly=TRUE)
    on.exit(par(opar))                   

    inc <- vector()
    inc[1] <- 1
    for (j in 1:7)
      {
        inc[j+1] <- floor(j * nrow(x$coef)/7)
      }
    labs <- rownames(x$coef)[inc]
    
    for (i in 1:ncol(x$coef))
      {
        crint <- 1.645 * as.vector(x$R[,i])^0.5
        m1 <- min(x$coef[,i]-crint,na.rm=TRUE)
        m2 <- max(x$coef[,i]+crint,na.rm=TRUE)

        plot(x$coef[,i],col="blue",ylim=c(m1,m2),axes=FALSE,
             xaxt='n',xlab='',ylab='',type="l",main="",lty=1)
        par(new=TRUE)
        plot(rep(0,length(x$coef[,i])),col="red",ylim=c(m1,m2),axes=FALSE,
             xaxt='n',xlab='',ylab='',type="l",main="",lty=1)
        par(new=TRUE)
        plot(x$coef[,i]+crint,col="black",ylim=c(m1,m2),axes=FALSE,
             xaxt='n',xlab='',ylab='',type="l",main="",lty=2)
        par(new=TRUE)
        plot(x$coef[,i]-crint,col="black",ylim=c(m1,m2),axes=TRUE,
            xaxt='n',xlab='',ylab='',type="l",main=colnames(x$coef)[i],lty=2)
            axis(1,at=inc,labels=labs)
      }
  }

