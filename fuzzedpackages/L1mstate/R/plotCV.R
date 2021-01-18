#############################################
#####  Plot likelihood path             #####
#############################################
plot.l1mstateCV = function(x,...){
  error.bars =function(object, upper, lower, width = 0.02, ...)
  {
    xlim = range(object)
    barw = diff(xlim) * width
    segments(object, upper, object, lower, ...)
    segments(object - barw, upper, object + barw, upper, ...)
    segments(object - barw, lower, object + barw, lower, ...)
    range(upper, lower)
  }
  
  res = x$fit
  plot.args=list(x=log(res[,1]),y=res[,2],ylim=range((res[,2]+res[,3]),(res[,2]-res[,3])),
                 xlab="log(Lambda)",ylab="CV Partial log-likelihood",type="n",
                 cex.axis = 1.25, font = 2, cex.lab = 1.5, col.lab = '#993333', font.lab=2)
  do.call("plot",plot.args)
  error.bars(log(res[,1]),res[,2]+res[,3],res[,2]-res[,3],width=0.01,col="darkgrey")
  points(log(res[,1]),res[,2],pch=20,col="red")
  lambda_pcvl = which(res[,4]=="pcvl")
  if(length(lambda_pcvl)==0){
    lambda_pcvl = which(res[,4]=="pcvl=cvmax")
    lambda_max = which(res[,4]=="pcvl=cvmax")
  }else{
    lambda_max = which(res[,4]=="cvmax")
  }
  abline(v=log(res[,1][lambda_max]),lty=3)
  abline(v=log(res[,1][lambda_pcvl]),lty=3)
  invisible()
}