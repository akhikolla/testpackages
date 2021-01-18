setClass(
  Class='statp',
  slots=c(
    stat='numeric',
    pvalue='numeric'
  )
)

setClass(
  Class='adjp',
  slots=c(
    stat='vector',
    pvalue='vector',
    adjp='vector'
  )
)



##' @title SV1F model with one factor simulation
##' @description Simulate stochastic volatility with one factor model (no jump) with given length and other parameters
##' @param M number of interverals to be simulated
##' @param m number of time points within each interval
##' @param p0 start price
##' @param mu drift
##' @param v0 volatility parameter
##' @param beta0 underlying Brownian motion intercept paramter
##' @param beta1 underlying Brownian motion slope parameter
##' @param alphav volatility parameter
##' @param cov Brownian motion correlation
##' @return simulated time series
##' @useDynLib JumpTest
##' @importFrom Rcpp evalCpp
##' @importFrom MASS mvrnorm
##' @examples SV1F(1200,390)
##' @references Chernov, M., et al. (2003). "Alternative models for stock price dynamics." Journal of Econometrics 116(1): 225-257.
##' @export
SV1F <- function(M,m,p0=3,mu=.03,v0=5,beta0=0,beta1=.125,alphav=-.1,cov=-.62)
{
  n <- M*m
  dt <- 1/M
  z <- mvrnorm(n+1,rep(0,2),matrix(c(1,cov,cov,1),2,2))
  pp <- pvc0(n, p0, mu*dt, beta0, beta1, v0, sqrt(dt), 1+alphav*dt, z)
  return(pp[-1])
}



##' @title SV1FJ model simulation
##' @description Simulate Stochastic Volitility model with one factor model (including jump) with given length and other parameters
##' @param M number of interverals to be simulated
##' @param m number of time points within each interval
##' @param p0 start price
##' @param lam frequency of jump
##' @param mu drift
##' @param v0 volatility parameter
##' @param beta0 underlying Brownian motion intercept paramter
##' @param beta1 underlying Brownian motion slope parameter
##' @param alphav volatility parameter
##' @param cov Brownian motion correlation
##' @return simulated time series
##' @useDynLib JumpTest
##' @importFrom Rcpp evalCpp
##' @importFrom MASS mvrnorm
##' @importFrom stats rpois rnorm
##' @examples SV1FJ(1000,390)
##' @references Chernov, M., et al. (2003). "Alternative models for stock price dynamics." Journal of Econometrics 116(1): 225-257.
##' @export
SV1FJ <- function(M,m,p0=3,lam=.2,mu=.03,v0=.5,beta0=0,beta1=.125,alphav=-.1,cov=-.62){
  n <- M*m
  dt <- 1/M
  z <- mvrnorm(n,rep(0,2),matrix(c(1,cov,cov,1),2,2))
  n1 <- rpois(n,lam*dt)
  jl <- which(n1>0)
  M0 <- sapply(n1[jl],rnorm,n=1,mean=0)
  M <- rep(0,n)
  M[jl] <- M0
  pp <- pvc(n, p0, mu*dt, beta0, beta1, v0, sqrt(dt), 1+alphav*dt, z, M)
  return(cbind(pp[-1], n1))
}

##' @title SV2F model simulation
##' @description Simulate Stochastic Volitility model with two factors model (no jump) with given length and other parameters
##' @param M number of interverals to be simulated
##' @param m number of time points within each interval
##' @param p.0 start price
##' @param mu drift
##' @param v.1 volatility parameter
##' @param v.2 volatility parameter
##' @param beta.0 underlying Brownian motion intercept paramter
##' @param beta.1 underlying Brownian motion slope parameter
##' @param beta.2 underlying Brownian motion slope parameter
##' @param alpha.1 volatility parameter
##' @param alpha.2 volatility parameter
##' @param beta.v2 second factor Brownian motion slope parameter
##' @param r1 correlation to first factor
##' @param r2 correlation to second factor
##' @return simulated time series
##' @useDynLib JumpTest
##' @importFrom Rcpp evalCpp
##' @importFrom MASS mvrnorm
##' @examples SV2F(1000,390)
##' @references Chernov, M., et al. (2003). "Alternative models for stock price dynamics." Journal of Econometrics 116(1): 225-257.
##' @export
SV2F<-function(M,m,p.0=3,mu=.03,v.1=.5,v.2=.5,beta.0=-1.2,beta.1=.04,beta.2=1.5,alpha.1=-.137*exp(-2),alpha.2=-1.386,beta.v2=.25,r1=-.3,r2=-.3){
  n <- M*m
  dt <- 1/M
  z <- mvrnorm(n,rep(0,3),matrix(c(1,r1,r2,r1,1,0,r2,0,1),3,3))
  mt <- mu*dt
  st <- sqrt(dt)
  v1xs <- 1+alpha.1*dt
  v2xs <- 1+alpha.2*dt
  pp2 <- pv2(n,mt,beta.0,beta.1,beta.2,p.0,v.1,v.2,st,z,v1xs,v2xs,beta.v2)
  return(pp2[-1])
}

##' @title SVJ model with one factor simulation
##' @description Simulate stochastic volatility model (with jump) with given length and other parameters
##' @param M number of interverals to be simulated
##' @param m number of time points within each interval
##' @param p0 start price
##' @param lambda frequency of jump
##' @param mu drift
##' @param v0 starting volatility
##' @param b volatility parameter
##' @param alpha volatility parameter
##' @param sigma volatility parameter
##' @param sigma1 jump size parameter
##' @return simulated time series
##' @useDynLib JumpTest
##' @importFrom Rcpp evalCpp
##' @importFrom stats rnorm rchisq
##' @examples SVJ(390,1200)
##' @references Yen, Y.-M. (2013). "Testing Jumps via False Discovery Rate Control." PloS one 8(4): e58365.
##' @export
SVJ <- function(M,m,p0=3,lambda=.2,mu=.05,v0=0,b=.2,alpha=.015,sigma=.05,sigma1=1){
  n <- M*m
  dt <- 1/M
  z0 <- rnorm(n)
  z1 <- rnorm(n)
  n1 <- sqrt(rpois(n,lambda*dt)*sigma1)
  M0 <- sapply(n1,rnorm,n=1,mean=0)
  d <- 4*b*alpha/sigma^2-1
  c0 <- .25*sigma^2*(1-exp(-alpha*dt))/alpha
  x0 <- rchisq(n,d)
  ppp <- lp(n,p0,mu,v0,dt,alpha,c0,z0,z1,M0,x0)
  return(cbind(ppp[-1],n1))
}

##' @title SV model with one factor simulation
##' @description Simulate stochastic volatility model (np jump) with given length and other parameters
##' @param M number of interverals to be simulated
##' @param m number of time points within each interval
##' @param p0 start price
##' @param mu drift
##' @param v0 starting volatility
##' @param b volatility parameter
##' @param alpha volatility parameter
##' @param sigma volatility parameter
##' @return simulated time series
##' @useDynLib JumpTest
##' @importFrom Rcpp evalCpp
##' @importFrom stats rnorm rchisq
##' @examples SV(390,1200)
##' @references Yen, Y.-M. (2013). "Testing Jumps via False Discovery Rate Control." PloS one 8(4): e58365.
##' @export
SV <- function(M,m,p0=3,mu=.05,v0=0,b=.2,alpha=.015,sigma=.05){
  n <- M*m
  dt <- 1/M
  z0 <- rnorm(n)
  z1 <- rnorm(n)
  d <- 4*b*alpha/sigma^2-1
  c0 <- .25*sigma^2*(1-exp(-alpha*dt))/alpha
  x0 <- rchisq(n,d)
  ppp <- lp2(n,p0,mu,v0,dt,alpha,c0,z0,z1,x0)
  return(ppp)
}






bns <- function(dp){
  mu43 <- 1/sqrt(pi)*2^(2/3)*gamma((4/3+1)/2)
  a <- (pi/2)^2+pi-5
  l <- length(dp)
  rv <- sum(dp^2)
  bv <- pi/2*l/(l-1)*sum(abs(dp[-length(dp)])*abs(dp[-1]))
  tp <- mu43^(-3)*l^2/(l-2)*sum((abs(dp[1:(length(dp)-2)])*abs(dp[2:(length(dp)-1)])*abs(dp[3:length(dp)]))^(4/3))
  b <- tp/bv^2
  z.rjt <- sqrt(l)*(rv-bv)/rv/sqrt(a*max(1,b))
  return(z.rjt)
}

amin <- function(r){
  RV <- sum(r^2)
  n <- length(r)
  r1 <- r[-1]
  r2 <- r[-n]
  minr <- pmin(abs(r1),abs(r2))
  mRV <- pi/(pi-2)*n/(n-1)*sum(minr^2)
  mRQ <- pi/(3*pi-8)*n^2/(n-1)*sum(minr^4)
  return((1-mRV/RV)/sqrt(1.81/n*max(1,mRQ/mRV^2)))
}

amed <- function(r){
  RV <- sum(r^2)
  n <- length(r)
  r1 <- r[3:n]
  r2 <- r[2:(n-1)]
  r3 <- r[1:(n-2)]
  medr <- apply(cbind(abs(r1),abs(r2),abs(r3)),1,median)
  meRV <- pi/(6-4*sqrt(3)+pi)*n/(n-2)*sum(medr^2)
  meRQ <- 3*pi/(9*pi+72-52*sqrt(3))*n^2/(n-2)*sum(medr^4)
  return((1-meRV/RV)/sqrt(.96/n*max(1,meRQ/meRV^2)))
}


##' @title Nonparametric jump test for a long period
##' @description perform nonparametric jump test for many intervals, and saved in vectors
##' @param retmat log return matrix, with intervals saved in columns
##' @param method jump test methods, chosen from "BNS", "Amed", and "Amin"
##' @return \item{stat}{test statistics}
##' @return \item{pvalue}{p-value}
##' @return \item{adjp}{adjusted p-values via 'BH' method}
##' @examples orip <- matrix(runif(3000),1000,3)
##' testres <- jumptestperiod(orip)
##' ts <- testres@stat
##' pv <- testres@pvalue
##' adjpv <- testres@adjp
##' @importFrom methods new
##' @importFrom stats pnorm p.adjust
##' @references Barndorff-Nielsen, O. E. and N. Shephard (2006). "Econometrics of testing for jumps in financial economics using bipower variation." Journal of financial Econometrics 4(1): 1-30.
##' @references Andersen, T. G., et al. (2012). "Jump-robust volatility estimation using nearest neighbor truncation." Journal of Econometrics 169(1): 75-93.
##' @references Dumitru, A.-M. and G. Urga (2012). "Identifying jumps in financial assets: a comparison between nonparametric jump tests." Journal of Business & Economic Statistics 30(2): 242-255.
##' @export
jumptestperiod <- function(retmat,method='BNS'){
  if(method=='BNS'){
    stat <- apply(retmat,2,bns)
  }else if(method=='Amed'){
    stat <- apply(retmat,2,amed)
  }else{
    stat <- apply(retmat,2,amin)
  }
  pvalue <- 1-pnorm(stat)
  adjp <- p.adjust(pvalue,'BH')
  return(new('adjp',stat=stat,pvalue=pvalue,adjp=adjp))
}


##' @title Nonparametric jump test for each interval
##' @description perform nonparametric jump test for each given interval (day)
##' @param ret log return vector
##' @param method jump test methods, chosen from "BNS", "Amed", and "Amin"
##' @return \item{stat}{test statistics}
##' @return \item{pvalue}{p-value}
##' @examples orip <- runif(100)
##' testres <- jumptestday(orip)
##' ts <- testres@stat
##' pv <- testres@pvalue
##' @importFrom methods new
##' @importFrom stats pnorm
##' @export
jumptestday <- function(ret,method='BNS'){
  if(method=='BNS'){
    stat <- bns(ret)
  }else if(method=='Amed'){
    stat <- amed(ret)
  }else{
    stat <- amin(ret)
  }
  pvalue <- 1-pnorm(stat)
  return(new('statp',stat=stat,pvalue=pvalue))
}





##' @title p-values matrix to be pooled
##' @description generate p-value matrix with given methods (at least 2)
##' @param retmat log return matrix by columns
##' @param method jump test methods, chosen from "BNS", "Amed", and "Amin"
##' @return a p-values matrix
##' @examples orip <- matrix(runif(3000),1000,3)
##' pmatrix <- pcombine(orip,c('BNS','Amed','Amin'))
##' @export
pcombine <- function(retmat,method){
  if(length(method)<2){
    stop('At least 2 method')
  }
  p0 <- list()
  for(i in 1:length(method)){
    p0[[i]] <- jumptestperiod(retmat,method[i])@pvalue
  }
  res <- do.call(cbind,p0)
  colnames(res) <- method
  return(res)
}



##' @title p-values pooling and adjustment
##' @description Pooling input p-values and perfrom FDR adjustments
##' @param pmat p-values matrix stored by columns
##' @param method pooling methods, see details
##' @details for p-values poolings, we provided six methods. "FI" for Fisher's method, "FD" for Fisher's with correlation adjustments, "SI" for Stouffer's method, "SD" for Stouffer's method with correlation adjustments, "MI" for minimum p-value methods, and "MA" for maximum p-value method
##' @return \item{stat}{pooled test statistcs}
##' @return \item{pvalue}{pooled p-values}
##' @return \item{adjp}{pooled p-values via "BH" adjustments}
##' @examples orip <- matrix(runif(3000),1000,3)
##' pvobj <- ppool(orip)
##' pvalue <- pvobj@pvalue
##' padjust <- pvobj@adjp
##' @importFrom methods new
##' @importFrom stats cor median p.adjust pbeta pchisq pnorm qnorm quantile
##' @references Benjamini, Y. and Y. Hochberg (1995). "Controlling the false discovery rate: a practical and powerful approach to multiple testing." Journal of the Royal Statistical Society. Series B (Methodological): 289-300.
##' @references Chang, L.-C., et al. (2013). "Meta-analysis methods for combining multiple expression profiles: comparisons, statistical characterization and an application guideline." BMC bioinformatics 14(1): 368.
##' @references Won, S., et al. (2009). "Choosing an optimal method to combine P-values." Statistics in medicine 28(11): 1537-1553.
##' @references Alves, G., & Yu, Y. K. (2014). Accuracy evaluation of the unified P-value from combining correlated P-values. PloS one, 9(3), e91225.
##' @export
ppool <- function(pmat,method='SD'){
  nmet <- ncol(pmat)
  pmat[pmat<1e-5] <- 1e-5
  pmat[pmat>1-1e-5] <- 1-1e-5
  if(method=='SD'){
    tf <- apply(pmat,1,max)
    cs <- qnorm(1-pmat)
    ccs <- cs[-which(tf<quantile(tf,.2)),]
    crs <- cor(ccs)
    er <- sum(crs[upper.tri(crs)])
    stat <- rowSums(cs)/sqrt(nmet+2*er)
    pv <- 1-pnorm(stat)
  }else if(method=='FD'){
    tf <- apply(pmat,1,max)
    cf <- -2*log(pmat)
    ccf <- cf[-which(tf<quantile(tf,.2)),]
    cr <- cor(ccf)
    cv <- 3.263*cr+.71*cr^2+.027*cr^3
    nmet2 <- 2*nmet
    cc <- (nmet2+sum(cv[upper.tri(cv)]))/nmet2
    ff <- nmet2^2/(nmet2+sum(cv[upper.tri(cv)]))
    stat <- rowSums(cf)/cc
    pv <- 1-pchisq(rowSums(cf)/cc,ff)
  }else if(method=='SI'){
    cs <- qnorm(1-pmat)
    stat <- rowSums(cs)/sqrt(nmet)
    pv <- 1-pnorm(stat)
  }else if(method=='FI'){
    stat <- -2*rowSums(log(pmat))
    pv <- 1-pchisq(stat,2*nmet)
  }else if(method=='MI'){
    stat <- apply(pmat,1,min)
    pv <- pbeta(stat,1,nmet)
  }else{
    stat <- apply(pmat,1,max)
    pv <- pbeta(stat,nmet,1)
  }
  pvadj <- p.adjust(pv,'BH')
  return(new('adjp',stat=stat,pvalue=pv,adjp=pvadj))
}