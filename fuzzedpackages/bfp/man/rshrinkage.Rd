\name{rshrinkage}
\alias{rshrinkage}

\title{Sample from the model-specific posterior of the shrinkage factor t = g / (1 + g)} 
\description{
  Sample from the model-specific posterior of the shrinkage factor t = g
  / (1 + g), using inverse sampling. This function is internal only
  because the arguments are not very user-friendly.
}
\usage{
rshrinkage(n=1, R2, nObs, p, alpha)
}

\arguments{
  \item{n}{number of samples (default: 1)}
  \item{R2}{coefficient of determination in the model}
  \item{nObs}{number of observations used to fit the model}
  \item{p}{number of design matrix columns without counting the
    intercept}
  \item{alpha}{used hyperparameter for hyper-g prior} 
}
\value{
  \code{n} posterior shrinkage factor samples
}
\author{Daniel Saban\'es Bov\'e}
\examples{
## construct a BayesMfp object
set.seed(29)
nObs <- 15

x1 <- rnorm (n=nObs)
x2 <- rbinom (n=nObs, size=20, prob=0.5) 
x3 <- rexp (n=nObs)

y <- rt (n=nObs, df=2)

test <- BayesMfp (y ~ bfp (x2, max = 3) + bfp(x3, max=3), nModels = 200, method="exhaustive")

## get the best found model
map <- test[[1]]

## produce many samples from the shrinkage factor
samples <- bfp:::rshrinkage(n=1e+5,
			    R2=map$R2,
		            nObs=nObs,
			    p=length(unlist(map$powers)),
			    alpha=4)
hist(samples,
     nclass=50, prob=TRUE)

## compare the mean with the analytical mean
abline(v=c(sampled <- mean(samples),
	   true <- map$postExpectedShrinkage),
       col=c(1, 2))
stopifnot(abs(sampled - true) < 0.0001)
}
\keyword{internal}
\keyword{distribution}

