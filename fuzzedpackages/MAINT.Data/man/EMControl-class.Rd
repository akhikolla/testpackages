\name{EMControl-class}
\docType{class}
\alias{EMControl-class}

\title{EM algorithm control parameters for fitting Gaussian mixtures to interval data.}

\description{This class contains the control parameters  for the EM algorithm  used in estimation of Gaussian mixtures by function \code{\link{Idtmclust}}.  .}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("EMControl", ...)}
or by calling the constructor-function \code{EMControl}.
}

\section{Slots}{
  \describe{
    \item{\code{nrep}}{Number of replications (different randomly generated starting points) of the EM algorithm.}
    \item{\code{maxiter}}{Maximum number of iterations in each replication of the EM algorithm.}
    \item{\code{convtol}}{Numeric tolerance for testing the convergence of the EM algorithm. Convergence is assumed when the log-likelihood changes less than convtol.}
    \item{\code{protol}}{Numeric tolerance for the mixture proportions. Proportions below protol, considered to be zero, are not allowed.}
    \item{\code{seed}}{Starting value for random generator.}
  }
}


\seealso{
  \code{\link{EMControl}}
}

\keyword{classes}
\keyword{EM algorithm}
\keyword{model based clustering}
\keyword{interval data}

