\name{handover}
\alias{handover}
\alias{handoff}
\alias{malpractice}
\docType{data}
\title{Dataset on communication breakdown in handover between physicians}
\description{

  Object \code{handover} is a likelihood function corresponding to a
  dataset arising from 69 medical malpractice claims and concerns
  handover (or hand-off) between physicians.  This dataset was analysed
  by Lin et al. (2009), and further analysed by Altham and Hankin
  (2010).  The computational methods are presented in the
  \sQuote{hyperdirichlet} and \sQuote{aylmer} packages and a further
  discussion is given in the \dQuote{integration} vignette of the
  \sQuote{hyper2} package.

}
\usage{data(handover)}
\details{
  Details are provided in the integration vignette, and the cited papers.
}
\references{
  \itemize{
    \item Y. Lin and S. Lipsitz and D. Sinha and A. A. Gawande and
    S. E. Regenbogen and C. C. Greenberg, 2009. \dQuote{Using {B}ayesian
      \eqn{p}-values in a \eqn{2\times 2}{2x2} table of matched pairs with
      incompletely classified data}. \emph{Journal of the Royal Statistical
      Society, Series C}, 58:2
    \item P. M. E. Altham and R. K. S. Hankin, 2010. \dQuote{Using
      recently developed software on a \eqn{2\times 2}{2x2} table of
      matched pairs with incompletely classified data}.  \emph{Journal of
      the Royal Statistical Society, series C}, 59(2): 377-379
    \item R. K. S. Hankin 2010. \dQuote{A generalization of the Dirichlet
      distribution}.   \emph{Journal of Statistical software}, 33:11
    \item L. J. West and R. K. S. Hankin 2008. \dQuote{Exact tests for
      two-way contingency tables with structural zeros}.   \emph{Journal
      of Statistical software}, 28:11
  }
}

\examples{
data(handover)
maxp(handover)
}
\keyword{datasets}
