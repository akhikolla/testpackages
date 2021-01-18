#' @importFrom minqa bobyqa
#' @importFrom stats dnorm
#' @importFrom stats optimize
#' @importFrom stats cor
#' @importFrom stats weighted.mean
polysFast <- function(x, M, w, ML=FALSE) {
  M <- as.numeric(as.factor(M))
  uM <- sort(unique(M))
  mapTheta <- function(v) {
    vv <- cumsum(c(v[1],exp(v[-1])))
    c(NA,-Inf,vv,Inf)
  }
  theta0 <- sapply(uM[-length(uM)],function(z) qnorm(weighted.mean(M<=z, w)) )
  if(ML) {
    #values = mainF(x, M, w, theta0)
    #theta0 <- imapThetaFast2(theta(M));
    #bob <- bobyqa(par=c(imapCor(cor(x,M)),imapTheta(theta0)), fn=optF(x,M,w))
    temp = w/sum(w)
    temp2 = fixxFast(x, temp)
    temp3 = sum(temp*dnorm(temp2,log=TRUE))
    bob <- suppressWarnings(bobyqa(par=c(atanh(cor(x,M)),imapThetaFast2(theta0)),
                  fn=optFFast, x=temp2,  w=temp, M=M, temp3 = temp3))
    return(  tanh(bob$par[1]))
  } else {
    temp = w/sum(w)
    temp2 = fixxFast(x, temp)
    temp3 = sum(temp*dnorm(temp2,log=TRUE))

    values = mainF(x, M, w, theta0)
     #opt <- optimize(mainF, interval = imapCorFast(cor(x,M)) + c(-3,3), x,w,temp1, temp2, temp3)
    opt <- suppressWarnings(optimize(optFcFast, interval=unlist(values[1]),
                    x=temp2, w=temp, theta0=(imapThetaFast2(theta0)),
                    M=M, temp3= temp3))
    return( tanh(opt$minimum) )
  }
}


