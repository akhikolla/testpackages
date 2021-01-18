##############################################################
### Plot transition hazards 

plot.cumhaz.l1mstate = function(x,type=c("single","separate"),cols,
                         xlab="Years since transplant",ylab="Cumulative hazard",ylim,lwd=3,lty,
                         legend,legend.pos,bty="o",...)
{
  #process x object obtained from l1ms_haz function
  x1 = x$Haz
  x1$time = x1$time/365
  unt = sort(unique(x1$time))
  K = unique(x1$trans)
  k=min(K)
  xk = x1[x1$trans==k,]
  xk_new = matrix(rep(0,3*length(unt)), nrow = length(unt))
  xk_new[,1] = unt
  xk_new[,2] = ifelse(xk_new[,1] %in% xk$time, xk$Haz, 0)
  xk_new[,2] = cumsum(xk_new[,2])
  xk_new[,3] = rep(k, length(unt))
  x1_new = xk_new
  
  for(k in K){
    if(k!=min(K)){
      xk = x1[x1$trans==k,]
      xk_new = matrix(rep(0,3*length(unt)), nrow = length(unt))
      xk_new[,1] = unt
      xk_new[,2] = ifelse(xk_new[,1] %in% xk$time, xk$Haz, 0)
      xk_new[,2] = cumsum(xk_new[,2])
      xk_new[,3] = rep(k, length(unt))
      x1_new = rbind(x1_new, xk_new)
    }
  }
  x1_new = data.frame(x1_new)
  names(x1_new)=c("time", "cumHaz", "trans")
  
  #now plot processed results of cumulative hazard
  msf1 = x1_new
  Q = max(msf1$trans)
  msft = unique(msf1$time) # the time points
  nt = length(msft)
  msfp = matrix(msf1[,2],nrow=nt) # the cumulative hazards in matrix 
  
  if (missing(legend)) legend = as.character(1:Q) 
  if (type=="single") {
    if (missing(cols)) cols = 1:Q
    if (missing(ylim)) ylim = c(0,max(msfp)+0.02)
    if (missing(lwd)) lwd = 1
    if (missing(lty)) lty = rep(1,Q)
    plot(msft,msfp[,1],type="s",ylim=ylim,xlab=xlab,ylab=ylab,col=cols[1],lwd=lwd,lty=lty[1],...)
    if(Q>1){
      for (q in 2:Q) lines(msft,msfp[,q],type="s",col=cols[q],lwd=lwd,lty=lty[q],...)
    }
    
    if (missing(legend.pos)) legend("topleft",legend=legend,col=cols,lwd=lwd,lty=lty,bty=bty)
    else legend(legend.pos[1],legend.pos[2],legend=legend,col=cols,lwd=lwd,lty=lty,bty=bty)
  }
  else if (type=="separate") {
    if (missing(cols)) cols = rep(1,Q)
    if (missing(lwd)) lwd = 1
    if (missing(lty)) lty = 1
    for (q in 1:Q) {
      if (missing(ylim)) plot(msft,msfp[,q],type="s",xlab=xlab,ylab=ylab,col=cols[q],lwd=lwd,...)
      else plot(msft,msfp[,q],type="s",ylim=ylim,xlab=xlab,ylab=ylab,col=cols[q],lwd=lwd,...)
      title(main=paste("Transition",q))
    }
  }
  return(invisible())
}
