\name{plot.bess}
\alias{plot.bess}
\title{Produces a coefficient profile plot of the coefficient or loss function paths}
\description{
  Produces a coefficient profile plot of the coefficient or loss paths for a fitted "bess" object.
}
\usage{
  \method{plot}{bess}(x, type=c("loss","coefficients","both"), breaks=TRUE, K=NULL, ...)
}
\arguments{
  \item{x}{a "bess" project}
  \item{type}{Either "both", "solutionPath" or "loss"}
  \item{breaks}{If TRUE, then vertical lines are drawn at each break point in the coefficient paths}
  \item{K}{which break point should the vertical lines drawn at}
  \item{\dots}{Other graphical parameters to plot}
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
#--------------linear model--------------#

data <- gen.data(500, 20, family = "gaussian", 10, 0.2, 1)
fit <- bess(data$x, data$y, family = "gaussian")
plot(fit, type = "both")

}
