\name{logLik.bess}
\alias{logLik.bess}
\title{Extract the loglikelihood from a "bess" object.}
\description{
  Similar to other logLik methods, which returns loglikelihood from a fitted "\code{bess}" object.
}
\usage{
  \method{logLik}{bess}(object,...)
}
\arguments{
  \item{object}{Output from the \code{bess} function.}
  \item{\dots}{Additional arguments affecting the predictions produced.}
}
\value{The value of the loglikelihood extracted from the "bess" object.}

\author{
  Canhong Wen, Aijun Zhang, Shijie Quan, and Xueqin Wang.
}
\references{
 Wen, C., Zhang, A., Quan, S. and Wang, X. (2020). BeSS: An R Package for Best Subset Selection in Linear, Logistic and Cox Proportional Hazards Models, \emph{Journal of Statistical Software}, Vol. 94(4).  doi:10.18637/jss.v094.i04.
}


\seealso{
\code{\link{bess}}, \code{\link{bess.one}}
}

\examples{

data <- gen.data(500, 20, family = "gaussian", 10, 0.2, 1)
fit <- bess(data$x, data$y, family = "gaussian")
logLik(fit)

}
