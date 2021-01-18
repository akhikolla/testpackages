\name{EditImputeCont-package}
\alias{EditImputeCont-package}
\alias{EditImputeCont}
\docType{package}
\title{
Simultaneous Edit-Imputation for Continuous Microdata Package
}
\description{
The package implements an integrated editing and imputation for continuous microdata under linear constraints. It relies on a Bayesian nonparametric hierarchical modeling approach in which the joint distribution of the data is estimated by a flexible joint probability model. The generated edit-imputed data are guaranteed to satisfy all imposed edit rules, whose types include ratio edits, balance edits and range restrictions.
}
\details{
\tabular{ll}{
Package: \tab EditImputeCont\cr
Type: \tab Package\cr
License: \tab GPL (>= 3)\cr
}
}
\author{
Quanli Wang, Hang J. Kim, Jerome P. Reiter, Lawrence H. Cox and Alan F. Karr

Maintainer: Hang J. Kim <hangkim0@gmail.com>
}
\references{
Hang J. Kim, Lawrence H. Cox, Alan F. Karr, Jerome P. Reiter and Quanli Wang (2015). "Simultaneous Edit-Imputation for Continuous Microdata", Journal of the American Statistical Association, DOI: 10.1080/01621459.2015.1040881.
}
\seealso{\code{\link{readData}}, \code{\link{createModel}}, \code{\link{multipleEI}}} 

\keyword{package}

\examples{
library(EditImputeCont)

## read the toy example data, which has two ratio edits and a balance edit
data(SimpleEx)
data1 = readData(Y.original=SimpleEx$D.obs, ratio=SimpleEx$Ratio.edit, 
	range=NULL, balance=SimpleEx$Balance.edit)

## create and initialize the model with 15 DP mixture components
# model1 = createModel(data.obj=data1, K=15)

## Run an iteration of MCMC

# model1$Iterate()

# dim(model1$Y.edited)
## [1]   1000   4  # Edit-imputed datasets of n=1000 records with p=4 variables

## Please see the example in the demo folder for more detailed explanation
}
