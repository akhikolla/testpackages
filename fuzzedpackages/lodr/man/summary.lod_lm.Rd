\name{summary.lod_lm}
\alias{summary.lod_lm}
\alias{print.summary.lod_lm}
\title{
  Summarizing Linear Model Fits with Covariates Subject to a Limit of Detection
}
\description{
  \code{summary} method for class "\code{lod_lm}"
}
\usage{
\method{summary}{lod_lm}(object, ...)

\method{print}{summary.lod_lm}(x, ...)
}
\arguments{
  \item{object}{An object of class "\code{lod_lm}", usually, a result of a call to \link[lodr]{lod_lm}}

  \item{x}{An object of class "\code{summary.lod_lm}", usually, a result of a call to \link[lodr]{summary.lod_lm}}

  \item{...}{further arguments passed to or from other methods.}

  }
\details{
    \code{print.summary.lod_lm} prints a table containing the coefficient estimates, standard errors, etc. from the \code{lod_lm} fit.
  }
\value{
  The function \code{summary.lod_lm} returns a list of summary statistics of the fitted linear model given in \code{object}, using the components (list elements) "\code{call}" and "\code{terms}" from its argument, plus

  \item{residuals}{residuals computed by \code{lod_lm}}
  \item{coefficients}{a \code{p x 4} matrix  for the estimated coefficient, its standard error, t-statistic and corresponding (two-sided) p-value.}
  \item{sigma}{the square root of the estimated variance of the random error.}
  \item{df}{degrees of freedom, a vector \code{(p, n-p)}, where \code{p} is the number of regression coefficients and n is the sample size of the data used in the model fitting}
}
\author{
Kevin Donovan, kmdono02@ad.unc.edu.

Maintainer: Kevin Donovan <kmdono02@ad.unc.edu>
}
\references{
  \insertRef{may}{lodr}
}

\seealso{
  The model fitting function \code{\link[lodr]{lod_lm}}, \code{\link[base]{summary}}.
}
\examples{
library(lodr)
## Using example dataset provided in lodr package: lod_data_ex
## 3 covariates: x1, x2, x3 with x2 and x3 subject to a lower limit of
## detection of 0

## nSamples set to 100 for computational speed/illustration purposes only.  
## At least 250 is recommended.  Same for boots=0; results in NAs returned for standard errors

fit <- lod_lm(data=lod_data_ex, frmla=y~x1+x2+x3, lod=c(0,0),
                  var_LOD=c("x2", "x3"), nSamples=100, boots=0)
 summary(fit)
}
