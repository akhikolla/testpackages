
% File  man/lm.br_internal.Rd
% Part of the lm.br package

\name{lm.br_internal}

\alias{Cpp_Clmbr}
\alias{Cpp_Clmbr-class}
\alias{Rcpp_Cpp_Clmbr-class}

\docType{class}


\title{Internal lm.br class \code{"Cpp_Clmbr"} }

\description{
Internal S4 class.

Class of data + model objects used for broken line regression.  An object instance stores pre-calculated quantities for a particular data set and model, so that methods can quickly compute exact-confidence intervals and regions for the changepoint.

}

\details{
  This class is not meant for the user.
}


\section{Extends}{
Class \code{"\linkS4class{C++Object}"}, directly.

All reference classes extend and inherit methods from \code{"\linkS4class{envRefClass}"}.

}


\section{Methods}{
  \describe{
    \item{\code{ci}:}{ printout confidence interval }
    \item{\code{cr}:}{ printout confidence region }
    \item{\code{mle}:}{ maximum likelihood estimates of parameters }
    \item{\code{sety}:}{ set new y-values in a Cpp_Clmbr object }
    \item{\code{sl}:}{ significance level for a postulate changepoint value }
  }
}


\seealso{ lm.br }

\examples{
showClass( Cpp_Clmbr )
}


\keyword{ classes }
\keyword{ internal }

