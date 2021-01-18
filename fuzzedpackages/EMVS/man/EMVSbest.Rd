\name{EMVSbest}
\alias{EMVSbest}
\title{Select the Best Model with EMVS}
\description{
\code{EMVSbest} outputs indices of the variables included in the model  with the highest posterior probability found.
  }
\usage{
EMVSbest(result)
}
\arguments{
  \item{result}{List object outputed by the EMVS procedure}
 }

\value{
\item{log_g_function}{The highest log-g-function found along the regularization path}
\item{indices}{The indices of the variables included in the best model found}
}

\author{Veronika Rockova\cr 
Maintainer: Veronika Rockova \email{vrockova@wharton.upenn.edu}}
\references{Rockova, V. and George, E. I.  (2014)
  \emph{EMVS: The EM Approach to Bayesian Variable Selection},   
\url{http://amstat.tandfonline.com/doi/abs/10.1080/01621459.2013.869223?journalCode=uasa20#preview}\cr
  \emph{Journal of the American Statistical Association}
}


\seealso{\code{EMVS}, \code{EMVSsummary}, \code{EMVSplot}}

\keyword{Bayesian variable selection}
\keyword{Spike and slab}
