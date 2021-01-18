\name{summary.NTA}
\alias{summary.NTA}
\alias{print.summary.NTA}
\title{
Summarizing the results of novel testing approach
}
\description{
\code{summary }method for class \code{"NTA"}.
}
\usage{
\method{summary}{NTA}(object, pless=0.05,...)
\method{print}{summary.NTA}(x, ...)
}

\arguments{
  \item{object}{an object of class \code{NTA}, a result of a call to \code{\link[vita]{NTA}}.}
 \item{pless}{ print only p-values less than pless.  Default is \code{pless=0.05}.}
  \item{x}{an object of class \code{summary.NTA}, a result of a call to \code{summary.NTA}.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
\code{print.summary.NTA} tries to be smart about formatting the permutation variable importance values,
  pvalue and  gives "significance stars".
}
\value{
\item{cmat}{a p x 2 matrix with columns for the \emph{original permutation variable importance} values
            and corresponding p-values.}
\item{pless}{ p-values less than pless}
\item{call}{the matched call to \code{NTA}.}
}


\seealso{
 \code{\link[vita]{NTA}}
 }


