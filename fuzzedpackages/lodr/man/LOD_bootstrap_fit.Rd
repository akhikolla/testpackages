\name{LOD_bootstrap_fit}
\alias{LOD_bootstrap_fit}
\title{
  Rcpp Code for Computing Standard Errors When Fitting Linear Models with Covariates Subject to a Limit of Detection (LOD)
}
\description{
  \code{LOD_bootstrap_fit} calls Rcpp code to compute linear model regression parameter standard errors in C++, taking into account covariates with limits of detection per the method detailed in \insertCite{may;textual}{lodr}.
}
\usage{
LOD_bootstrap_fit(num_of_boots, y_data, x_data, no_of_samples, threshold, 
max_iterations, LOD_u_l)

}
\arguments{
  \item{num_of_boots}{number denoting the number of bootstrap resamples to use to compute the regression parameter standard errors.}
  \item{y_data}{numeric vector consisting of data of the model's outcome variable.}
  \item{x_data}{column-named matrix consisting of data of the model's covariates with each column representing one covariate, with values outside of the limit(s) of detection marked as \code{NA}.  A columns of ones must be included if the model has an intercept term.  Note that for valid inference, order of the covariates/columns in the matrix must be as follows from left to right: those with no LOD followed by those with an LOD.}
  \item{no_of_samples}{an integer specifying the number of samples to generate for each subject with covariate values outside of their limits of detection.  For more details, see \insertCite{may;textual}{lodr}.}
  \item{threshold}{number denoting the minimum difference in the regression parameter estimates needed for convergence of the model fitting procedure.}
  \item{max_iterations}{number denoting the maximum number of iterations allowed in the model fitting procedure.}
  \item{LOD_u_l}{numeric matrix consisting of the lower and upper limits of detection for all covariates in the model as the columns, with each covariate containing its own row, in the same order as the covariates in \code{x_data}.  If no limit of detection exists, the corresponding matrix entry is marked with an \code{NA}.  An entry for the intercept (\code{NA} in each column) must be included if applicable.}
  }
\details{
  This function is used to complete the standard error computations done when fitting a linear model by calling \link[lodr]{lod_lm}; the standard error computations are done in C++ to minimize computation time.
}
\value{
  \code{LOD_bootstrap_fit} returns a list which each component being a numeric vector consisting of the last iteration's regression parameter estimates when fitting the model on a bootstrap resample of the input data.
}
\author{
Kevin Donovan, kmdono02@ad.unc.edu.

Maintainer: Kevin Donovan <kmdono02@ad.unc.edu>
}
\references{
  \insertRef{may}{lodr}
}

\seealso{
  \code{\link[lodr]{lod_lm}} is the recommended function for fitting a linear model with covariates subject to limits of detection, which uses \code{LOD_fit}.  \code{\link[lodr]{LOD_fit}} is used to compute the regression parameter estimates.
}
\examples{
library(lodr)
## Using example dataset provided in lodr package: lod_data_ex
## 3 covariates: x1, x2, x3 with x2 and x3 subject to a lower limit of
## detection of 0

# Replace values marked as under limit of detection using 0 with NA,
# add column of ones for intercept
lod_data_with_int <-
  as.matrix(cbind("Intercept"=rep(1, dim(lod_data_ex)[1]), lod_data_ex))

lod_data_ex_edit <-
  apply(lod_data_with_int, MARGIN = 2, FUN=function(x){ifelse(x==0, NA, x)})

# Fit model with bootstrap procedure, report regression parameter estimate standard errors
LOD_matrix <- cbind(c(NA, NA, -100, -100), c(NA, NA, 0, 0))

## no_of_samples set to 50 for computational speed/illustration purposes only.  
## At least 250 is recommended.  
## Same for num_of_boots=5; at least 25 is recommended
 
bootstrap_fit_object <-
  LOD_bootstrap_fit(num_of_boots=5, y_data=lod_data_ex_edit[,2],
                    x_data=lod_data_ex_edit[,-2],
                    no_of_samples=50,
                    threshold=0.001, max_iterations=100, LOD_u_l=LOD_matrix)

boot_SEs <- apply(do.call("rbind", bootstrap_fit_object), 2, sd)
names(boot_SEs) <- names(lod_data_with_int[,-2])

boot_SEs
}
