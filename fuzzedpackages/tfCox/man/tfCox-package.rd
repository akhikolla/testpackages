\name{tfCox-package}
\alias{tfCox-package}
\docType{package}
\title{
Fit the Additive Trend Filtering Cox Model
}
\description{
This package is called tfCox or trend filtering for Cox model, which is proposed in Jiacheng Wu & Daniela Witten (2019) Flexible and Interpretable Models for Survival Data, Journal of Computational and Graphical Statistics, DOI: 10.1080/10618600.2019.1592758. It provides an approach to fit additive Cox model in which each component function is estimated to be piecewise polynomial with adaptively-chosen knots.

Function \code{\link{tfCox}} fits the trend filtering Cox model for a range of tuning parameters. Function \code{\link{cv_tfCox}} returns the optimal tuning parameter selected by K-fold cross validation.
}

\details{
\tabular{ll}{
Package: \tab tfCox\cr
Type: \tab Package\cr
Version: \tab 0.1.0\cr
Date: \tab 2019-05-20\cr
License: \tab GPL (>= 2)\cr
}
The package includes the following functions:
\code{\link{tfCox}}, \code{\link{cv_tfCox}}, \code{\link{plot.tfCox}}, \code{\link{plot.cv_tfCox}}, \code{\link{predict.tfCox}}, \code{\link{summary.tfCox}}, \code{\link{summary.cv_tfCox}}, \code{\link{sim_dat}}, \code{\link{plot.sim_dat}}.
}
\author{
Jiacheng Wu
Maintainer: Jiacheng Wu <wujiacheng1992@gmail.com>
}
\references{
Jiacheng Wu & Daniela Witten (2019) Flexible and Interpretable Models for Survival Data, Journal of Computational and Graphical Statistics, DOI: 10.1080/10618600.2019.1592758
}
\keyword{ package }
