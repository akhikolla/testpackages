\name{cond.fun.chisq.test}
\alias{cond.fun.chisq.test}

\title{
Conditional Functional Chi-Squared Test
}
\description{
Asymptotic chi-squared test to determine the model-free functional dependency of effect variable \eqn{Y} on a cause variable \eqn{X}, conditioned on a third variable \eqn{Z}.
}

\usage{
cond.fun.chisq.test(x, y, z=NULL, data=NULL, log.p = FALSE,
                    method = c("fchisq", "nfchisq"))
}

\arguments{
  \item{x}{
vector or character; either a discrete random variable (cause) represented as vector or a character column name in \code{data}.
}
  \item{y}{
vector or character; either a discrete random variable (effect) represented as vector or a character column name in \code{data}.
}
  \item{z}{
vector or character; either a discrete random variable (condition) represented as vector or a character column name in \code{data}. In case of \code{NULL}, a \code{fun.chisq.test} on a contingency table, with \code{x} as row variable and \code{y} as column variable, is returned. See \code{\link{fun.chisq.test}} for details. The default is \code{NULL}.
}

  \item{data}{
  a data frame containing three or more columns whose names can be used as values for \code{x}, \code{y} and \code{z}. In case of \code{NULL}, \code{x}, \code{y} and \code{z} must be vectors. The default is \code{NULL}.
}

\item{log.p}{
  logical; if \code{TRUE}, the p-value is given as \code{log(p)}. Taking the log improves numerical precision when the p-value is close to zero. The default is \code{FALSE}.
  }

  \item{method}{
  a character string to specify the method to compute the conditional functional chi-squared test statistic and its p-value. The options are \code{"fchisq"} (default) and \code{"nfchisq"}. See Details.
}

}
\details{
The conditional functional chi-squared test introduces the concept of conditional functional depedency, where the functional association between two variables (\code{x} and \code{y}) is tested conditioned on a third variable (\code{z}) \insertCite{zhang2014nonparametric}{FunChisq}. Two methods are provided to compute the chi-squared statistic and its p-value. When \code{method = "fchisq"}, the p-value is computed using the chi-squared distribution; when \code{method = "nfchisq"}, a normalized statistic is obtained by shifting and scaling the original chi-squared statistic and a p-value is computed using the standard normal distribution \insertCite{Box2005}{FunChisq}. The normalized test is more conservative on the degrees of freedom.
}

\value{
A list with class "\code{htest}" containing the following components:

\item{statistic}{the conditional functional chi-squared statistic if \code{method = "fchisq"}; or the normalized conditional functional chi-squared statistic if \code{method = "nfchisq"}.}

\item{parameter}{degrees of freedom for the conditional functional chi-squared statistic.}

\item{p.value}{p-value of the conditional functional test. If \code{method = "fchisq"}, the p-value is computed by an asymptotic chi-squared distribution; if \code{method = "nfchisq"}, the p-value is computed by the standard normal distribution.}

\item{estimate}{an estimate of the conditional function index between 0 and 1. The value of 1 indicates strong functional dependency between \code{x} and \code{y}, given \code{z}. It is asymmetrical with respect to whether \code{x} was chosen as the cause of effect \code{y} or vice versa.}
}

\references{

\insertAllCited{}

}

\author{
Sajal Kumar and Mingzhou Song
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
See (unconditional) functional chi-squared test \code{\link{fun.chisq.test}}.
}

\examples{
# Generate a relationship between variables X and Z
xz = matrix(c(30,2,2, 2,2,40, 2,30,2),ncol=3,nrow=3,
            byrow = TRUE)
# Re-construct X
x = rep(c(1:nrow(xz)),rowSums(xz))
# Re-construct Z
z = c()
for(i in 1:nrow(xz))
  z = c(z,rep(c(1:ncol(xz)),xz[i,]))

# Generate a relationship between variables Z and Y
# Make sure Z retains its distribution
zy = matrix(c(4,30, 30,4, 4,40),ncol=2,nrow=3,
            byrow = TRUE)
# Re-construct Y
y = rep(0,length(z))
for(i in unique(z))
  y[z==i] = rep(c(1:ncol(zy)),zy[i,])

# Tables
table(x,z)
table(z,y)
table(x,y)

# Conditional functional dependency
# Y = f(X) | Z should be false
cond.fun.chisq.test(x=x,y=y,z=z)
# Z = f(X) | Y should be true
cond.fun.chisq.test(x=x,y=z,z=y)
# Y = f(Z) | X should be true
cond.fun.chisq.test(x=z,y=y,z=x)
}

\keyword{conditional functional dependency}
