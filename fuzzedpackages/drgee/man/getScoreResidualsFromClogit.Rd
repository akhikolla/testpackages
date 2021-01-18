\encoding{latin1}
\name{getScoreResidualsFromClogit}
\alias{getScoreResidualsFromClogit}

\title{
Get the (conditional) score residuals from conditional logistic regression
}
\description{
  A function to obtain the conditional score residuals from conditional
  logistic regression as well as the sum of derivatives of these
  residuals with respect to the parameters evaluated at the estimated
  parameters.
}
\usage{
getScoreResidualsFromClogit(fit, coefs, y, x, id)
}
\arguments{
  \item{fit}{
    A fitted object. Could be any object that has a
    class.function \code{coefficients}, but is intended to be used with
    \code{clogit} in the \code{survival} package.
  }
  \item{coefs}{
    An estimate of the parameters. If the argument \code{fit} is supplied,
    this will be replaced by \code{coefficients(fit)}.
  }
  \item{y}{
    The observed outcomes.
  }
  \item{x}{
    The design matrix.
  }
  \item{id}{
    The cluster identification variable.
  }
}
\details{

  \code{getScoreResidualsFromClogit} calculates the residuals from the
  conditional score equations used in conditional logistic regression as
  we as the sums of  their derivatives with respect to the paramaters.

  This is useful if one wants to use a calculate the sandwich estimate
  of a variance where the uncertainty in the estimation of a
  conditional odds ratio needs to be taken into account.
}
\value{

  \code{getScoreResidualsFromClogit} return a list containing:

  \item{U }{
    An \eqn{n} x \eqn{p} matrix, where \eqn{n} is the number of observations
    and \eqn{p} is the length of the estimated parameter.}

  \item{dU.sum }{
    An \eqn{p} x \eqn{p} matrix containing the sums of the
    derivatives. The rows correspond to the estimating equations and
    the columns correspond to the partial derivatives.}
  
  \item{dU.mean}{The mean cluster sum of the derivatives of the score
    residuals from the conditional likelihood. 
  }

}

\author{
Johan Zetterqvist, Arvid \enc{Sjölander}{Sjolander}
}

\seealso{
  \code{\link{robVcov}}.
}

\keyword{cluster}
