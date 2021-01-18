\name{vcov.mlegc}
\alias{vcov.mlegc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Covariance Matrix of the Maximum Likelihood Estimates
}


\usage{
\method{vcov}{mlegc}(object, digits = max(3, getOption("digits") - 3), ...)
          }

\description{
	Calculate covariance and correlation matrix.
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class \code{mlegc} inherited from function \code{\link{mlegc}}.}

  \item{digits}{integer indicating the number of decimal places (round) or
  significant digits (signif) to be used. }
  
  \item{...}{other arguments passed to \code{vcov}.}

}


\value{
The estimated variance-covariance matrix and estimated correlation matrix.
}

\seealso{
\code{\link{mlegc}}
}

\author{
Zifei Han \email{hanzifei1@gmail.com}
}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Estimation}
