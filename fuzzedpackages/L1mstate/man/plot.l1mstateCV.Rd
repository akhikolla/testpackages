\name{plot.l1mstateCV}
\alias{plot.l1mstateCV}
\title{
plot the cross-validation curve produced by cv.l1mstateR
}
\description{
Plots the cross-validation curve, and upper and lower standard deviation curves, as a function of the lambda values used.
}
\usage{
\method{plot}{l1mstateCV}(x,...)
}
\arguments{
  \item{x}{fitted cv.l1mstateR object}
  \item{\dots}{Further arguments to plot}
}
\value{No return value.}
\details{
A plot is produced, and nothing is returned.
}

\author{
Xuan Dang \email{xuandang11289@gmail.com}
}

\examples{
set.seed(1001)
p <- 9
times <- 1:p
rho <- 0.5
H <- abs(outer(times, times, "-"))
C <- 1 * rho^H
C[cbind(1:p, 1:p)] <- C[cbind(1:p, 1:p)] 
sigma <- matrix(C,p,p)
mu <- rep(0,p)

beta12 <- c(-.65,-.65,-.65,0,0,-.65,-.65,0,0)
beta13 <- c(-.65,-.65,0,0,0,0,-.65,0,0) 
beta23 <- c(0,-.65,-.65,0,0,-.65,-.65,0,-.65)

N <- 200
x <- mvrnorm(n=N, mu, sigma)
col_names  <-  c(sprintf("X\%d", seq(1,dim(x)[2])))
colnames(x)  <-  col_names
N12  <-  N-2*N\%/\%5
N13 <- N-N12
hx12 <- exp(x[1:N12,] \%*\% beta12) 
ty12 <- rexp(N12,hx12) 
hx13 <- exp(x[(N12+1):N,] \%*\% beta13) 
ty13 <- rexp(N13,hx13)
hx23 <- exp(x[1:N12,] \%*\% beta23) 
ty23 <- rexp(N12,hx23) 
t2 <- rep(0,N)
t2[1:N12]  <-  ty12
t2[(N12+1):N]  <-  ty13
s2  <-  rep(0,N)
s2[1:N12] <-  rep(1,N12)
t3 <- rep(0,N)
t3[1:N12]  <-  ty12+ty23
t3[(N12+1):N]  <-  ty13
#censor case for state 3
tcens3 <- rbinom(n=N, prob = 0.3, size = 1)
s3 <- 1-tcens3
    
# pre-process data into long-format using msprep
dt <- data.frame(illt=t2, ills=s2, dt=t3, ds=s3, x)
tmat <- matrix(c(NA,NA,NA,1,NA,NA,2,3,NA), nrow = 3)
longdt <- msprep(time=c(NA,"illt","dt"),status=c(NA,"ills","ds"),
                 keep = col_names, data=dt,trans=tmat)

cv.l1fit <- cv.l1mstateR(longdt, nlambda = 100, nfolds = 10)
plot.l1mstateCV(cv.l1fit)
}

\keyword{ multi-state models}% use one of  RShowDoc("KEYWORDS")
\keyword{ regularization }% __ONLY ONE__ keyword per line
