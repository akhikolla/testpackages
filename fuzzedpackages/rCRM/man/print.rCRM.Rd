\name{print.rCRM}
\alias{print.rCRM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print a rCRM Object
}
\description{
Print a summary of results
}
\usage{
\method{print}{rCRM}(x, digits = 3, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{fitted rCRM object
}
  \item{digits}{number of digits in printout
}
  \item{\dots}{additional print arguments
}
}
\details{
The performed model is printed, followed by the estimated probability of DLT from a fitted \code{rCRM} object.
}
\value{
The data frame above is silently returned
}
\author{
Maintainer: Xiang Li <xli256@its.jnj.com>
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{rCRM}}
}
\examples{
set.seed(1213)

dose0=c(1:6)
prob0=c(0.007, 0.086, 0.294, 0.545, 0.731, 0.841)

m=3; Y=NULL; X=NULL
for (i in 1:length(dose0)) {
  Y=c(Y, rbinom(m, size=1, prob=prob0[i]))
  X=c(X, rep(i, m))
}

fiti=rCRM(X, Y, dose0, tp=0.3, mlambda=10)
# attributes(fiti)
fiti
}

\keyword{Print}


