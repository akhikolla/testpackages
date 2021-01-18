\name{LOD_fit}
\alias{LOD_fit}
\title{
  Rcpp Code for Fitting Linear Models with Covariates Subject to a Limit of Detection (LOD)
}
\description{
  \code{LOD_fit} calls Rcpp code to compute linear model regression parameter estimates in C++, taking into account covariates with limits of detection per the method detailed in \insertCite{may;textual}{lodr}.
}
\usage{
LOD_fit(y_data, x_data, mean_x_preds, beta, sigma_2_y, sigma_x_preds, no_of_samples, 
threshold, max_iterations, LOD_u_l)

}
\arguments{
  \item{y_data}{numeric vector consisting of data of the model's outcome variable}
  \item{x_data}{column-named matrix consisting of data of the model's covariates with each column representing one covariate, with values outside of the limit(s) of detection marked as \code{NA}.  A columns of ones must be included if the model has an intercept term.}
  \item{mean_x_preds}{numeric vector consisting of initial estimates of the means for each covariate, in the same order as the covariates in \code{x_data}}
  \item{beta}{numeric vector consisting of initial estimates of the regression parameters for each covariate, in the same order as the covariates in \code{x_data}}
  \item{sigma_2_y}{an initial estimate of the variance of the outcome variable}
  \item{sigma_x_preds}{numeric matrix consisting of an initial estimate of the covariance matrix for the model's covariates, in the same order as the covariates in \code{x_data}}
  \item{no_of_samples}{an integer specifying the number of samples to generate for each subject with covariate values outside of their limits of detection.  For more details, see \insertCite{may;textual}{lodr}.}
  \item{threshold}{number denoting the minimum difference in the regression parameter estimates needed for convergence of the model fitting procedure.}
  \item{max_iterations}{number denoting the maximum number of iterations allowed in the model fitting procedure.}
  \item{LOD_u_l}{numeric matrix consisting of the lower and upper limits of detection for all covariates in the model as the columns, with each covariate containing its own row, in the same order as the covariates in \code{x_data}.  If no limit of detection exists, the corresponding matrix entry is marked with an \code{NA}.  An entry for the intercept (\code{NA} in each column) must be included if applicable.}
  }
\details{
  This function is used to complete the model fitting computations done when calling \link[lodr]{lod_lm}; the fitting computations are done in C++ to minimize computation time.
}
\value{
  \code{LOD_fit} returns a list containing the following components:

  \item{y_expand_last_int}{a numeric vector consisting of the outcome data with duplicate entries for subjects with covariates outside of their limits of detection per the corresponding resampling procedure, from the last iteration of the model fitting procedure.}
  \item{x_data_return_last_int}{a numeric matrix consisting of the covariate data with sampled values for covariates of subjects with covariates outside of their limits of detection, from the last iteration of the model fitting procedure.}
  \item{beta_estimates}{a numeric matrix consisting of the regression parameter estimates from each iteration of the model fitting procedure.}
  \item{beta_estimate_last_iteration}{a numeric vector consisting of the regression parameter estimates from the last iteration of the model fitting procedure.}
}
\author{
Kevin Donovan, kmdono02@ad.unc.edu.

Maintainer: Kevin Donovan <kmdono02@ad.unc.edu>
}
\references{
  \insertRef{may}{lodr}
}

\seealso{
  \code{\link[lodr]{lod_lm}} is the recommended function for fitting a linear model with covariates subject to limits of detection, which uses \code{LOD_fit}.  \code{\link[lodr]{LOD_bootstrap_fit}} is used to compute regression parameter estimate standard errors using bootstrap resampling.
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
  data.frame(apply(lod_data_with_int, MARGIN = 2, FUN=function(x){ifelse(x==0, NA, x)}))

# Fit linear model to dataset with only subjects without covariates under
# limit of detection to get initial estimate for the regression parameters.
beta_inital_est <- coef(lm(y~x1+x2+x3, data=lod_data_ex_edit))

# Get initial estimates of mean vector and covariance matrix for covariates and variance of outcome,
# again using data from subjects without covariates under limit of detection

mean_x_inital <- colMeans(lod_data_ex_edit[,c(-1,-2)], na.rm = TRUE)
sigma_x_inital <- cov(lod_data_ex_edit[,c(-1,-2)], use="pairwise.complete.obs")
sigma_2_y_inital <- sigma(lm(y~x1+x2+x3, data=lod_data_ex_edit))^2

# Fit model, report regression parameter estimates from last iteration
LOD_matrix <- cbind(c(NA, NA, -100, -100), c(NA, NA, 0, 0))

## no_of_samples set to 100 for computational speed/illustration purposes only.  
## At least 250 is recommended.
 
fit_object <-
LOD_fit(y_data=lod_data_ex_edit[,2],
        x_data=as.matrix(lod_data_ex_edit[,-2]),
        mean_x_preds=mean_x_inital, beta=beta_inital_est, sigma_2_y=sigma_2_y_inital,
        sigma_x_preds=sigma_x_inital, no_of_samples=100,
        threshold=0.001, max_iterations=100, LOD_u_l=LOD_matrix)

fit_object$beta_estimate_last_iteration
}