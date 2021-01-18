\name{summary.PimpTest}
\alias{summary.PimpTest}
\alias{print.summary.PimpTest}
\title{
Summarizing PIMP-algorithm outcomes
}
\description{
\code{summary }method for class \code{"PimpTest"}.
}
\usage{
\method{summary}{PimpTest}(object, pless=0.05,...)
\method{print}{summary.PimpTest}(x, ...)
}
\arguments{
 \item{object}{an object of class \code{PimpTest}, a result of a call to \code{\link[vita]{PimpTest}}.}
 \item{pless}{ print only p-values less than pless.  Default is \code{pless=0.05}.}
  \item{x}{an object of class \code{summary.PimpTest}, a result of a call to \code{summary.PimpTest}.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
\code{print.summary.PimpTest} tries to be smart about formatting the VarImp, pvalue etc. and  gives "significance stars".
}
\value{
\item{cmat}{a p x 3 matrix with columns for the \code{mean(PerVarImp)},\code{sd(PerVarImp)} and
            the the p-values of the Kolmogorov-Smirnov Tests. }
\item{cmat2}{a p x 2 matrix with columns for the \emph{original permutation variable importance} values
            and corresponding p-value.}
\item{para }{ Shall the null distribution be modelled by a Gaussian distribution? }
\item{pless}{ p-values less than pless}
\item{call}{the matched call to \code{PimpTest}.}
\item{call.PIMP}{the matched call to \code{PIMP}.}
\item{type}{one of regression, classification}
}


\seealso{
\code{\link[vita]{PimpTest}}, \code{\link[vita]{PIMP}}
}

