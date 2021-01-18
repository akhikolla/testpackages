\name{biochem}
\alias{biochem}
\docType{data}
\title{Biochemical Data}
\description{
The data consists of levels of biochemical components in brains of mice.
The treatment group received a drug and the control group a placebo.
}
\usage{data(biochem)}

\format{
  A data frame with 22 observations on the following 3 variables.
  \describe{
    \item{\code{comp.1}}{First biochemical component}
    \item{\code{comp.2}}{Second biochemical component}
    \item{\code{group}}{Factor with the two levels 'Control' and 'Treat'}
  }
}

\source{Table 1 in Brown and Hettmansperger (1987).}

\references{
\cite{Brown, B. M. and Hettmansperger, T. P. (1987). Affine invariant rank methods in the
bivariate location model. \emph{Journal of the Royal Statistical Society, Series B}, \bold{49}, 301--310.}
}



\examples{
data(biochem)
ojaMedian(biochem[,1:2])
}
\keyword{datasets}

