\name{bet.plot}
\alias{bet.plot}
\title{Plotting Binary Expansion Testing (2-dimensions)}
\usage{
bet.plot(X, d, unif.margin = FALSE, cex=0.5, ...)
}
\arguments{
\item{X}{a matrix with two columns.}

\item{d}{depth of BET.}

\item{unif.margin}{
  logicals. If \code{TRUE} the data has been uniformed based on empirical cumulative distribution function. Default to be \code{FALSE} and the function uniforms the data.
  }

\item{cex}{number indicating the amount by which plotting text and symbols should be scaled relative to the default.}

\item{...}{graphical parameters to plot}
}
\description{
\code{bet.plot} shows the cross interaction of the strongest asymmetry, which the BET returns with the rejection of independence null.
This function only works for the test on two variables, that is, \code{X} can only have two columns.
There are \eqn{2^{2d} - 1} nontrivial binary variables in the \eqn{\sigma}-field and \eqn{(2^d - 1)^2} of them are cross interactions, whose positive regions are in plotted in white and whose negative regions are plotted in blue.
\code{plot.bet} shows the cross interaction where the difference of number of observations in the positive and negative region is largest.
}
\examples{
v <- runif(128, -pi, pi)
X1 <- cos(v) + 2.5 * rnorm(128, 0, 1/20)
X2 <- sin(v) + 2.5 * rnorm(128, 0, 1/20)
bet.plot(cbind(X1, X2), 3)
}
