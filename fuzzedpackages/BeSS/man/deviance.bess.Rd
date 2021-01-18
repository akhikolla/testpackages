\name{deviance.bess}
\alias{deviance.bess}
\title{Extract the deviance from a "bess" object.}
\description{
  Similar to other deviance methods, which returns deviance from a fitted "\code{bess}" object.
}
\usage{
  \method{deviance}{bess}(object,...)
}
\arguments{
  \item{object}{Output from the \code{bess} function or the \code{bess.one} function.}
  \item{\dots}{Additional arguments affecting the predictions produced.}
}
\value{The value of the deviance extracted from the "bess" object.}

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
deviance(fit)

}
