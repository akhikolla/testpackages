\name{TSGS}
\alias{TSGS}
\title{Two-Step Generalized S-Estimator for cell- and case-wise outliers}
\description{
  Computes the Two-Step Generalized S-Estimate (2SGS) -- a robust estimate of location and scatter
  for data with cell-wise and case-wise contamination.
}
\usage{
TSGS(x, filter=c("UBF-DDC","UBF","DDC","UF"),
    partial.impute=FALSE, tol=1e-4, maxiter=150, method=c("bisquare","rocke"),
    init=c("emve","qc","huber","imputed","emve_c"), mu0, S0)
}
\arguments{
  \item{x}{a matrix or data frame.}
  \item{filter}{the filter to be used in the first step (see Leung et al. (2016)). Default is 'UBF-DDC'. For all filters, the default parameters are used.}
  \item{partial.impute}{whether partial imputation is used prior to estimation (see details). The default is \code{FALSE}.}
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
}
\details{
   This function computes 2SGS as described in Agostinelli et al. (2015) and Leung and Zamar (2016). The procedure has two major steps:

   In Step I, the method filters (i.e., flags and removes) cell-wise outliers using Gervini-Yohai univariate filter (Agostinelli et al., 2015)
   or univariate-bivariate filter (Leung et al., 2016) or univariate-bivariate-plus-DDC filter (Leung et al., 2016; Rousseeuw and Van den Bossche, 2016).
   The filtering step can be called on its own by using the function \code{\link{gy.filt}} or \code{\link[cellWise]{DDC}}.

   In Step II, the method applies GSE or GRE (GSE with a Rocke-type loss function), which has been specifically designed to deal with
   incomplete multivariate data with case-wise outliers, to the filted data coming from Step I. The second step can be called on its own
   by using the function \code{\link{GSE}}.

   The 2SGS method is intended for continuous variables, and requires that the number of observations
   n be relatively larger than 5 times the number of variables p for desirable performance (see the rejoinder in Agostinelli et al., 2015).
   In our numerical studies, for n too small relative to p, 2SGS may experience a lack of convergence, especially for filtered data
   sets with a proportion of complete observations less than
   1/2 + (p+1)/n. To overcome this problem, partial imputation prior to estimation is proposed
   (see the rejoinder in Agostinelli et al., 2015). The procedure is rather ad hoc, but initial numerical experiements
   show that partial imputation may work. Further research on this topic is still needed.
   By default, partial imputation is not used, unless specified.

   In general, we warn users to use 2SGS with caution for data set with n relatively smaller than 5 times p.

   The application to the chemical data set analyzed in Agostinelli et al. (2015) can be found in \code{\link{geochem}}.

   The tools that were used to generate contaminated data in the simulation study in Agostinelli et al. (2015) can be found in \code{\link{generate.cellcontam}} and \code{\link{generate.casecontam}}.
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

    Leung, A., Yohai, V.J., Zamar, R.H. (2016). Multivariate Location and Scatter Matrix Estimation Under Cellwise and Casewise Contamination.
    arXiv:1609.00402.

    Rousseeuw P.J., Van den Bossche W. (2016). Detecting deviating data cells. arXiv:1601.07251
}
\author{Andy Leung \email{andy.leung@stat.ubc.ca}, Claudio Agostinelli, Ruben H. Zamar, Victor J. Yohai}
\seealso{
  \code{\link{GSE}}, \code{\link{gy.filt}}, \code{\link[cellWise]{DDC}}, \code{\link{generate.cellcontam}}, \code{\link{generate.casecontam}}
}
\examples{
set.seed(12345)

# Generate 5% cell-wise contaminated normal data
# using a random correlation matrix with condition number 100
x <- generate.cellcontam(n=100, p=10, cond=100, contam.size=5, contam.prop=0.05)

## Using MLE
slrt( cov(x$x), x$A)

## Using Fast-S
slrt( rrcov:::CovSest(x$x)@cov, x$A)

## Using 2SGS
slrt( TSGS(x$x)@S, x$A)


# Generate 5% case-wise contaminated normal data
# using a random correlation matrix with condition number 100
x <- generate.casecontam(n=100, p=10, cond=100, contam.size=15, contam.prop=0.05)

## Using MLE
slrt( cov(x$x), x$A)

## Using Fast-S
slrt( rrcov:::CovSest(x$x)@cov, x$A)

## Using 2SGS
slrt( TSGS(x$x)@S, x$A)
}
