plot.bess=function(x,type=c("loss","coefficients","both"),breaks=TRUE,K=NULL, ...)
{
  object=x
  type <- match.arg(type)
  s.list=object$s.list

  if(is.null(K)) K=s.list[length(s.list)]
  if(object$family=="bess_gaussian") dev=object$mse else dev=object$deviance

  beta=object$beta
  s_order=order(s.list)
  s.list=s.list[s_order]
  dev=dev[s_order]
  beta=beta[,s_order]
  beta=cbind(rep(0,nrow(object$beta)),beta)

  if(type=="loss")
  {
    plot_loss(dev,s.list,K,breaks, mar = c(3,4,3,4))
  }
  if(type=="coefficients")
  {
    plot_solution(beta,c(0, s.list),K,breaks, mar = c(3,4,3,4))
  }
  if(type=="both")
  {
    layout(matrix(c(1,2),2,1,byrow=TRUE),heights=c(0.45,0.55), widths=1)
    oldpar <- par(las=1, mar=c(2,4,2,4), oma=c(2.5,0.5,1.5,0.5))
    plot_loss(dev,s.list,K,breaks,show_x = FALSE)
    plot_solution(beta, c(0, s.list), K,breaks)
    par(oldpar)
    par(mfrow=c(1,1))
  }
}


plot_loss <- function(loss,df,K,breaks=TRUE,show_x=TRUE, mar = c(0,4,2,4)){

  plot.new()                            # empty plot
  plot.window(range(df), range(loss), xaxs="i")
  oldpar <- par(mar = mar,              # no bottom spacing
                lend="square")          # square line ends
  par(new=TRUE)                         # add to the plot
  if(show_x)
  {
    plot(df, loss, type = "b", ylab=expression(L(beta)),
       xlim=c(0,max(df)))
  }else
  {
    plot(df, loss, type = "b", ylab=expression(L(beta)),
         xlim=c(0,max(df)), xaxt='n')
  }
  title(xlab='Model size', line = 2)
  if(breaks)abline(v=K, col="orange", lwd=1.5, lty=2) ## add a vertical line
  grid()
  axis(2)
  #axis(4, pos=par("usr")[1], line=0.5)  # this would plot them 'inside'
  # box()                                 # outer box

  par(oldpar)
}

plot_solution <- function(beta, df, K, breaks = TRUE, mar = c(3,4,0,4)){
  p <- nrow(beta)
  plot.new()                            # empty plot
  plot.window(range(df), range(beta), xaxs="i")

  oldpar <- par(mar=mar,         # no top spacing
                lend="square")          # square line ends
  par(new=TRUE)                         # add to the plot

  plot(df, beta[1,], type="l",col=1, xlim=c(0,max(df)),xlab="",
       ylim=range(beta),ylab=expression(beta))
  title(xlab='Model size', line = 2)
  for(i in 2:p){
    lines(df, beta[i,], col=i,xlim=c(0,p+1))
  }
  if(breaks) abline(v=K, col="orange", lwd=1.5, lty=2) ## add a vertical line
  #matplot(df, t(beta), lty = 1, ylab="",  xaxs="i",type = "l",xlim=c(0,p+1))

  nnz = p
  xpos = max(df)-0.8
  pos = 4
  xpos = rep(xpos, nnz)
  ypos = beta[, ncol(beta)]
  text(xpos, ypos, 1:p, cex = 0.8, pos = pos)

  grid()
  axis(2)
  box()                                 # outer box

  par(oldpar)
}


