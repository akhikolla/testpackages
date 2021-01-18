\name{ImpS}
\alias{ImpS}
\title{Imputed S-estimator}
\description{
  Computes the simple three-step estimator as described in the rejoinder of Agostinelli et al. (2015).
}
\usage{
ImpS(x, alpha=0.95, method=c("bisquare","rocke"), init=c("emve","emve_c"), ...)
}
\arguments{
  \item{x}{a matrix or data frame.}
  \item{alpha}{quantile of the reference distribution in the univariate filter step (see \code{\link{gy.filt}}). Default is 0.95.}
  \item{method}{which loss function to use: 'bisquare', 'rocke'.}
  \item{init}{
      type of initial estimator. Currently this can either be "emve" (EMVE with uniform sampling, see Danilov et al., 2012) or
      "emve_c" (EMVE_C with cluster sampling, see Leung and Zamar, 2016). Default is "emve".
  }
  \item{...}{optional, additional arguments to be passed to \code{\link{GSE}}.}
}
\details{
   This function computes the simple three-step estimator as described in the rejoinder in Agostinelli et al. (2015). The procedure has three steps:

   In Step I, the method flags and removes cell-wise outliers using the Gervini-Yohai univariate only filter (see \code{\link{gy.filt}}).

   In Step II, the method imputes the filtered cells using coordinate-wise medians. 
   
   In Step III, the method applies MVE-S to the filtered and imputed data from Step II (see \code{\link{GSE}}).

}
\value{
    The following gives the major slots in the output S4 object: 
    \tabular{ll}{
        \code{mu} \tab Estimated location. Can be accessed via \code{\link{getLocation}}. \cr
        \code{S} \tab Estimated scatter matrix. Can be accessed via \code{\link{getScatter}}. \cr 
        \code{xf} \tab Filtered data matrix from the first step of 2SGS. Can be accessed via \code{\link{getFiltDat}}. \cr
    }
}
\references{
    Agostinelli, C., Leung, A. , Yohai, V.J., and Zamar, R.H. (2015) Robust estimation of multivariate location and scatter in the presence of cellwise and casewise contamination. TEST.
}
\author{Andy Leung \email{andy.leung@stat.ubc.ca}, Claudio Agostinelli, Ruben H. Zamar, Victor J. Yohai}
\seealso{
  \code{\link{GSE}}, \code{\link{gy.filt}}
}
