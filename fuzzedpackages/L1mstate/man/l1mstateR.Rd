\name{l1mstateR}
\alias{l1mstateR}
\title{
Fit multi-state models with lasso regularization
}
\description{
Fit a multi-state models via penalized partial likelihood. The regularization path is computed for the lasso at a path of values for the regularization paramter lambda. Can deal with right-censoring and left-truncated data.
}
\usage{
l1mstateR(longdt, lambda=NULL, nlambda=100, rlambda=NULL, thresh=1e-7, maxit=1e+5)
}
\arguments{
  \item{longdt}{input in long-format structure}
  \item{lambda}{A user supplied lambda sequence. Typical usage is to have the program compute its own lambda sequence based on nlambda and rlambda. Supplying a value of lambda overrides this. WARNING: use with care. Avoid supplying a single value for lambda.  Supply instead a decreasing sequence of lambda values. l1mstateR relies on its warms starts for speed, and its often faster to fit a whole path than compute a single fit.}
  \item{nlambda}{The number of lambda values- default is 100.}
  \item{rlambda}{Smallest value for lambda, as a fraction of the maximum lambda, the (data derived) entry value (i.e. the smallest value for which all coefficients are zero). The default depends on the sample size relative to the number of covariates. If sample size>#covariates, the default is {0.001}, close to zero.  If sample size>#covariates, the default is 0.01.}
  \item{thresh}{Convergence threshold for one-step coordinate descent. Defaults value is 1E-7.}
  \item{maxit}{Maximum number of passes over the data for all lambda values; default is 1E+5.}
}
\details{
  Please note that input has to be long-format structure.
}
\value{
  \item{aBetaSTD}{A list of coefficients in standardized form, each one corresponds to each lambda value.}
  \item{aBetaO}{A list of coefficients in original form, each one corresponds to each lambda value.}
  \item{fit}{A matrix of lambda values and log-likelihood values}
  \item{numcovs}{Number of covariates}
  \item{numtrans}{Number of transitions}
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

l1fit <- l1mstateR(longdt, nlambda = 100)
}
\keyword{ multi-state models}% use one of  RShowDoc("KEYWORDS")
\keyword{ regularization }% __ONLY ONE__ keyword per line
