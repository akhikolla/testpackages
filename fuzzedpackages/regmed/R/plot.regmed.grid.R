plot.regmed.grid <- function(x, as.log=FALSE, ...){
  ## input : 
  ## x  is output of regmed.grid
  ## log indicates if log Lambda should be plotted for mediator coefficient plot

  if(class(x) != "regmed.grid") {
      stop("input not regmed.grid class")
  }
    
  ### create log lambda ###

  lambda<-x$grid.dat$lambda
  ### fix up log lambda when lambda is 0, this is just for plotting purpose ###
  log.lambda<-ifelse(lambda>0,log(lambda),log(min(lambda[lambda>0])/2))

  ### plot of BIC vs. Lambda ###

  if(as.log){

    plot(log.lambda, x$grid.data$bic, xlab="Log Lambda", ylab="BIC", ...)

  }else{

    plot(lambda, x$grid.data$bic, xlab="Lambda", ylab="BIC", ...)

  }

  ### plot of coefficients of mediators vs. Lambda ###

  ### all alphas and betas over grid in single matrices ###
  p.beta<-do.call(cbind,lapply(x$fit.list, function(x) x$beta))
  p.alpha<-do.call(cbind,lapply(x$fit.list, function(x) x$alpha))

  plot.lambda <- lambda
  ### change log lambda so lines are plotted in the same direction as lambda ###
  if(as.log) plot.lambda<-min(log.lambda)-log.lambda
		
  matplot(c(-plot.lambda,plot.lambda),
          rbind(t(p.beta),t(p.alpha)),col=1:nrow(p.beta),type="n",
          xlab=paste0(ifelse(as.log,"Log ",""),"Lambda"),ylab="Coefficients",axes=F,
          xlim=c(-max(abs(plot.lambda)),max(abs(plot.lambda))), ...)
  abline(v=0,lwd=2)
  abline(h=0,lwd=2)

  matplot(c(-plot.lambda,plot.lambda),
          rbind(t(p.beta),t(p.alpha)),col=1:nrow(p.beta),type="l",
          xlab="Lambda",ylab="Coefficients",axes=FALSE,add=TRUE, ...)

  ### labels for x-axis ticks ###
  x.ticks<-abs(axTicks(1))
  if(as.log) x.ticks<-round(min(log.lambda) + abs(axTicks(1)))

  axis(side=1,at=axTicks(1),labels=x.ticks)
  axis(side=2)
  mtext(side=3,at=par("usr")[1:2]/2,c("Beta","Alpha"), ...)

  box(bty="o", ...)

  invisible()
}
