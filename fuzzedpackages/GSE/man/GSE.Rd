\name{GSE}
\alias{GSE}
\title{ Generalized S-Estimator in the presence of missing data}
\description{
  Computes the Generalized S-Estimate (GSE) -- a robust estimate of location and scatter 
  for data with contamination and missingness. 
}
\usage{
GSE(x, tol=1e-4, maxiter=150, method=c("bisquare","rocke"), 
    init=c("emve","qc","huber","imputed","emve_c"), mu0, S0, ...)
}
\arguments{
  \item{x}{a matrix or data frame. May contain missing values, but cannot contain columns with completely missing entries. }
  \item{tol}{tolerance for the convergence criterion. Default is 1e-4.}
  \item{maxiter}{maximum number of iterations for the GSE algorithm. Default is 150.}
  \item{method}{which loss function to use: 'bisquare', 'rocke'.}
  \item{init}{
      type of initial estimator. Currently this can either be "emve" (EMVE with uniform sampling, see Danilov et al., 2012),
      "qc" (QC, see Danilov et al., 2012), "huber" (Huber Pairwise, see Danilov et al., 2012), 
      "imputed" (Imputed S-estimator, see the rejoinder in Agostinelli et al., 2015), or
      "emve_c" (EMVE_C with cluster sampling, see Leung and Zamar, 2016).
      Default is "emve". If \code{mu0} and \code{S0} are provided, this argument is ignored.
  }
  \item{mu0}{optional vector of initial location estimate}
  \item{S0}{optional matrix of initial scatter estimate}
  \item{...}{optional arguments for computing the initial estimates (see \code{\link{emve}}, \code{\link{HuberPairwise}}).}
}
\details{
   This function computes GSE (Danilov et al., 2012) and GRE (Leung and Zamar, 2016). The estimator requires a robust positive definite 
   initial estimator. This initial estimator is required to ``re-scale" the partial square mahalanobis distance for 
   the different missing pattern, in which a single scale parameter is not enough. This function currently allows two 
   main initial estimators: EMVE (the default; see \code{\link{emve}} and Huberized Pairwise 
   (see \code{\link{HuberPairwise}}). GSE using Huberized Pairwise with sign psi function is referred to as QGSE in Danilov et al. (2012). 
   Numerical results have shown that GSE with EMVE as 
   initial has better performance (in both efficiency and robustness), but computing time can be longer. 
}
\value{
  An S4 object of class \code{\link{GSE-class}} which is a subclass of the virtual class \code{\link{CovRobMissSc-class}}. The
  output S4 object contains the following slots:
  
    \tabular{ll}{
        \code{mu} \tab Estimated location. Can be accessed via \code{\link{getLocation}}. \cr
        \code{S} \tab Estimated scatter matrix. Can be accessed via \code{\link{getScatter}}. \cr 
        \code{sc} \tab Generalized S-scale (GS-scale). Can be accessed via \code{\link{getScale}}. \cr
        \code{pmd} \tab Squared partial Mahalanobis distances. Can be accessed via \code{\link{getDist}}. \cr
        \code{pmd.adj} \tab Adjusted squared partial Mahalanobis distances. Can be accessed via \code{\link{getDistAdj}}. \cr
        \code{pu} \tab Dimension of the observed entries for each case. Can be accessed via \code{\link{getDim}}. \cr
        \code{mu0} \tab Estimated initial location.\cr
        \code{S0} \tab Estimated initial scatter matrix. \cr 
        \code{ximp} \tab Input data matrix with missing values imputed using best linear predictor. Not meant to be accessed. \cr
        \code{weights} \tab Weights used in the estimation of the location. Not meant to be accessed. \cr
        \code{weightsp} \tab First derivative of the weights used in the estimation of the location. Not meant to be accessed. \cr
        \code{iter} \tab Number of iterations till convergence. Not meant to be accessed. \cr
        \code{eps} \tab relative change of the GS-scale at convergence. Not meant to be accessed. \cr
        \code{call} \tab Object of class \code{"language"}. Not meant to be accessed. \cr
        \code{x} \tab Input data matrix. Not meant to be accessed. \cr
        \code{p} \tab Column dimension of input data matrix. Not meant to be accessed. \cr
        \code{estimator}\tab Character string of the name of the estimator used. Not meant to be accessed. \cr
    }
}
\references{
    Agostinelli, C., Leung, A. , Yohai, V.J., and Zamar, R.H. (2015) Robust estimation of multivariate location and scatter in the presence of cellwise and casewise contamination. \emph{TEST}.

    Danilov, M., Yohai, V.J., Zamar, R.H. (2012). Robust Esimation of Multivariate
    Location and Scatter in the Presence of Missing Data. \emph{Journal of the American Statistical Association}
    \bold{107}, 1178--1186.
    
    Leung, A. and Zamar, R.H. (2016). Multivariate Location and Scatter Matrix Estimation Under Cellwise and Casewise Contamination.
    Submitted.    
}
\author{Andy Leung \email{andy.leung@stat.ubc.ca}, Ruben H. Zamar, Mike Danilov, Victor J. Yohai}
\seealso{
  \code{\link{emve}}, \code{\link{HuberPairwise}}, \code{\link{GSE-class}}, \code{\link{generate.casecontam}}
}
\examples{
set.seed(12)

## generate 10-dimensional data with 10% casewise contamination
n <- 100
p <- 10
A <- matrix(0.9, p, p)
diag(A) <- 1
x <- generate.casecontam(n, p, cond=100, contam.size=10, contam.prop=0.1, A=A)$x

## introduce 5% missingness
pmiss <- 0.05
nmiss <- matrix(rbinom(n*p,1,pmiss), n,p)
x[ which( nmiss == 1 ) ] <- NA

## Using EMVE as initial
res.emve <- GSE(x)
slrt( getScatter(res.emve), A) ## LRT distances to the true covariance

## Using QC as initial
res.qc <- GSE(x, init="qc")
slrt( getScatter(res.qc), A) ## in general performs worse than if EMVE used as initials

}
