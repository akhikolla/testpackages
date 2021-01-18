

plot.mixest <- function(x, ...)
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

    col <- 1:ncol(x$rvi)
    
    par(xpd=TRUE)
    
    if (x$parameters[2]==0 | x$parameters[2]==1) 
      {
        for (i in 1:((ncol(x$rvi)-1)))
          {
            plot(x$rvi[,i],type="l",col=col[i],ylim=c(0,1),axes=FALSE,
                 xaxt='n',xlab="",ylab="",main="")
            par(new=TRUE)
          }
        plot(x$rvi[,i+1],type="l",col=col[i+1],ylim=c(0,1),axes=TRUE,
             xaxt='n',xlab="",ylab="",main="RVI")
        axis(1,at=inc,labels=labs)
        legend('bottom',inset=c(0,-0.35),colnames(x$coef),lty=rep(1,(i+1)),
               col=col,ncol=6,cex=0.6) 
      }  
    else
      {
        temp <- replace(x$rvi,which(x$rvi==0),NA)
        for (i in 1:((ncol(x$rvi)-1)))
          {
            plot(i/ncol(temp)*temp[,i],type="p",pch=15,col=col[i],ylim=c(0,1),
                 axes=FALSE,xaxt='n',xlab="",ylab="",main="")
            par(new=TRUE)
          }
        plot(temp[,i+1],type="p",pch=15,col=col[i+1],ylim=c(0,1),
             axes=TRUE,xaxt='n',yaxt='n',xlab="",ylab="",main="inclusion of variables")
        axis(1,at=inc,labels=labs)
        legend('bottom',inset=c(0,-0.35),colnames(x$coef),pch=rep(c(15),(i+1)), 
               col=col[1:(i+1)],ncol=6,cex=0.6) 
      }
 
    for (i in 1:(ncol(x$weights)-1))
      {
        plot(x$weights[,i],col=i,ylim=c(0,1),axes=FALSE,
             xaxt='n',xlab='',ylab='',type="l",main="",lty=1)
        par(new=TRUE)
      }
    plot(x$weights[,i+1],col=i+1,ylim=c(0,1),axes=TRUE,
         xaxt='n',xlab='',ylab='',type="l",main="weights",lty=1)
    axis(1,at=inc,labels=labs)
  }

