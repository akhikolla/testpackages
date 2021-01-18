\name{sumlog}
\alias{sumlog}
\title{Log-sum from log data}
\description{
Obtain log(sum(x)) from log(x), without passing to exponentials. 
It is based on the fact that 
log(a + b) = log(a) + log (1 + exp(log(b) - log(a))). 
}
\usage{
sumlog(x,lower=-745,upper=709)
}
\arguments{
  \item{x}{Vector of log-values}
  \item{lower}{Value such that exp(lower-epsilon)=0}
  \item{upper}{Value such that exp(upper+epsilon)=Inf}
}
\details{
   This function computes the logarithm of the sum of exp(x), without
   passing through exponentials. It shall be used to avoid under/over
   flow. It has proven useful in computing the likelihood of finite
   mixture models, normalization constants, importance sampling, etc. 
   It is described in the appendix of Farcomeni (2012).}
\value{A scalar equal to \code{log(sum(exp(x)))}.}
\references{
Farcomeni, A. (2012)  Quantile Regression for longitudinal data based
on latent Markov subject-specific parameters. \emph{Statistics and Computing}, \bold{22}, 141-152
}
\author{Alessio Farcomeni \email{alessio.farcomeni@uniroma1.it}, Andy
Leung \email{andy.leung@stat.ubc.ca}}
\examples{
# complete underflow without sumlog 
x <- c(-750,-752)
log(sum(exp(x)))
sumlog(x)

# imprecise sum 
x <- c(-745,-752)
log(sum(exp(x)))
sumlog(x)

# no issues 
x <- c(log(3),log(2))
log(5)
log(sum(exp(x)))
sumlog(x)
}


