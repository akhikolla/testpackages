\name{optimMethods}
\alias{optimMethods}
\title{Optimization Methods (or Algorithms) for the \code{mle}
Method}

\usage{

optimMethods(optimFun = c("both", "nloptr::nloptr", "stats::optim"))

}
\arguments{

  \item{optimFun}{Value of the corresponding formal argument
  of the \code{mle} method, or \code{"both"}. In the later
  case the full list of algorithms will be obtained.}

}
\value{

  A data frame with two columns: \code{optimFun} and \code{optimMethod}.

}
\section{Caution}{

  The optimisation method given in the argument \code{optimMethod} of
  the \code{mle} method sould be compliant with the \code{compGrad}
  argument. Only a small number of possibilies have been tested,
  including the default values.

}
\description{
  
  Optimization methods (or algorithms) for the \code{mle} method.

}
\seealso{

  \code{\link{mle-methods}}.
}
\examples{
optimMethods()
}
