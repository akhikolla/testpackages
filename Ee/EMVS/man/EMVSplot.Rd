\name{EMVSplot}
\alias{EMVSplot}
\title{Spike-and-slab Dynamic Posterior Exploration}
\description{
\code{EMVSplot} procedure plots the solution path of the estimated regression coefficients (posterior modes) for different \code{v0} values. 
  }
\usage{
EMVSplot(result, plot_type=c("both","reg","gf"), omit.zeroes = FALSE)
}
\arguments{
  \item{result}{List object outputed by the EMVS procedure}
  \item{plot_type}{Plot type: "both" for plotting both the regularization path together with the associated log g function, "reg" only for the regulaziation plot, "gf" only for the log g function.
}
  \item{omit.zeroes}{Logical: \code{TRUE} or \code{FALSE}. If \code{TRUE}, only the selected coefficients are plotted, the remaining coefficients set to zero}
 }

\details{
Coefficients that are not thresholded out are depicted in blue, the rest in red. Log g function computed only for models with at most 1 000 predictors.
} 

\author{Veronika Rockova\cr 
Maintainer: Veronika Rockova \email{vrockova@wharton.upenn.edu}}
\references{Rockova, V. and George, E. I.  (2014)
  \emph{EMVS: The EM Approach to Bayesian Variable Selection},   
\url{http://amstat.tandfonline.com/doi/abs/10.1080/01621459.2013.869223?journalCode=uasa20#preview}\cr
  \emph{Journal of the American Statistical Association}
}


\seealso{\code{EMVS}, \code{EMVSsummary}, \code{EMVSbest}}

\keyword{Bayesian variable selection}
\keyword{Spike and slab}
