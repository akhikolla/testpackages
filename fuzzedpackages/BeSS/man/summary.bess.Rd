\name{summary.bess}
\alias{summary.bess}
\title{summary method for a "bess" object}
\description{
  Print a summary of the "\code{bess}" object.
}
\usage{
  \method{summary}{bess}(object, ...)
}
\arguments{
  \item{object}{a "\code{bess}" object}
  \item{\dots}{additional print arguments}
}

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
summary(fit)

}
