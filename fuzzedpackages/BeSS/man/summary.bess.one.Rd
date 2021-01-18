\name{summary.bess.one}
\alias{summary.bess.one}
\title{summary method for a "bess.one" object}
\description{
  Print a summary of the "\code{bess.one}" object.
}
\usage{
  \method{summary}{bess.one}(object, ...)
}
\arguments{
  \item{object}{a "\code{bess.one}" object}
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
fit <- bess.one(data$x, data$y, s = 10, family = "gaussian")
summary(fit)

}
