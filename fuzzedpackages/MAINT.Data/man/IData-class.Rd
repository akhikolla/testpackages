\name{IData-class}
\Rdversion{1.1}
\docType{class}
\alias{IData-class}
\alias{head,IData-method}
\alias{print,IData-method}
\alias{summary,IData-method}
\alias{tail,IData-method}
\alias{MidPoints,IData-method}
\alias{LogRanges,IData-method}
\alias{Ranges,IData-method}
\alias{NbMicroUnits,IData-method}

\title{Class IData}
\description{
  A data-array of interval-valued data is an array where each of the NObs rows, corresponding to each entity under analysis, contains the observed intervals of the NIVar descriptive variables.
}

\section{Slots}{
  \describe{
    \item{\code{MidP}:}{A data-frame of the midpoints of the observed intervals}
    \item{\code{LogR}:}{A data-frame of the logarithms of the ranges of the observed intervals}
    \item{\code{ObsNames}:}{An optional vector of names assigned to the individual observations.}
    \item{\code{VarNames}:}{An optional vector of names to be assigned to the Interval-valued Variables.}
    \item{\code{NObs}:}{Number of entities under analysis (cases)}
    \item{\code{NIVar}:}{Number of interval variables}
    \item{\code{NbMicroUnits}:}{An integer vector with the number of micro data units by interval-valued observation (or an empty vector, if not applicable)}
  }
}

\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "IData")}: show S4 method for the IData-class.}
    \item{nrow}{\code{signature(x = "IData")}: returns the number of statistical units (observations).}
    \item{ncol}{\code{signature(x = "IData")}: returns the number of of Interval-valued variables.}
    \item{dim}{\code{signature(x = "IData")}: returns a vector with the of number statistical units as first element, and the number of Interval-valued variables as second element.}
    \item{rownames}{\code{signature(x = "IData")}: returns the row (entity) names for an object of class IData.}
    \item{colnames}{\code{signature(x = "IData")}: returns column (variable) names for an object of class IData.}
    \item{names}{\code{signature(x = "IData")}: returns column (variable) names for an object of class IData.}
    \item{MidPoints}{\code{signature(Idt = "IData")}: returns a data frame with MidPoints for an object of class IData.}
    \item{LogRanges}{\code{signature(Idt = "IData")}: returns a data frame with LogRanges for an object of class IData.}
    \item{Ranges}{\code{signature(Idt = "IData")}: returns an data frame with Ranges for an object of class IData.}
    \item{NbMicroUnits}{\code{signature(Idt = "IData")}: returns an integer vector with the number of micro data units by interval-valued observation for an object of class IData.}
    \item{head}{\code{signature(x = "IData")}: head S4 method for the IData-class.}
    \item{tail}{\code{signature(x = "IData")}: tail S4 method for the IData-class.}
    \item{plot}{\code{signature(x = "IData")}: plot S4 methods for the IData-class.}
    \item{mle}{\code{signature(x = "IData")}: Maximum likelihood estimation.}
    \item{fasttle}{\code{signature(x = "IData")}: Fast trimmed maximum likelihood estimation.}
    \item{fulltle}{\code{signature(x = "IData")}: Exact trimmed maximum likelihood estimation.}
    \item{RobMxtDEst}{\code{signature(x = "IData")}: Robust estimation of distribution mixtures for interval-valued data.}
    \item{MANOVA}{\code{signature(x = "IData")}: MANOVA tests on the interval-valued data.}
    \item{lda}{\code{signature(x = "IData")}: Linear Discriminant Analysis using maximum likelihood parameter estimates of Gaussian mixtures.}
    \item{qda}{\code{signature(x = "IData")}: Quadratic Discriminant Analysis using maximum likelihood parameter estimates of Gaussian mixtures.}
    \item{Roblda}{\code{signature(x = "IData")}: Linear Discriminant Analysis using robust estimates of location and scatter.}
    \item{Robqda}{\code{signature(x = "IData")}: Quadratic Discriminant Analysis using robust estimates of location and scatter.}
    \item{snda}{\code{signature(x = "IData")}: Discriminant Analysis using maximum likelihood parameter estimates of SkewNormal mixtures.}
  }
}

\author{Pedro Duarte Silva <psilva@porto.ucp.pt>\cr 
        Paula Brito <mpbrito.fep.up.pt>
}

\references{
Azzalini, A. and Dalla Valle, A. (1996), The multivariate skew-normal distribution. \emph{Biometrika} \bold{83}(4), 715--726.\cr

Brito, P., Duarte Silva, A. P. (2012), Modelling Interval Data with Normal and Skew-Normal Distributions. \emph{Journal of Applied Statistics} \bold{39}(1), 3--20.\cr

Duarte Silva, A.P., Filzmoser, P. and Brito, P. (2017), Outlier detection in interval data. \emph{Advances in Data Analysis and Classification}, 1--38.\cr

Noirhomme-Fraiture, M., Brito, P. (2011), Far Beyond the Classical Data Models: Symbolic Data Analysis. \emph{Statistical Analysis and Data Mining} \bold{4}(2), 157--170.
}

\seealso{
\code{\link{IData}}, \code{\link{AgrMcDt}}, \code{\link{mle}}, \code{\link{fasttle}}, \code{\link{fulltle}}, \code{\link{RobMxtDEst}}, 
\code{\link{MANOVA}}, \code{\link{lda}}, \code{\link{qda}}, \code{\link{Roblda}}, \code{\link{Robqda}}
}

\keyword{classes}
\keyword{interval data}
