\name{fccaXX}
\alias{fccaXX}
\title{
Canonical correlation analysis between two groups of mixed functional and scalar variables
}
\description{
This function carries out the canonical correlation analysis between two groups of mixed functional and scalar variables. Three different representing methods can be used for the functional coefficients. The tuning parameters should be specified in the arguments \code{control1} and \code{control2} for the two groups \code{xL1} and \code{xL2}, respectively.
}
\usage{
fccaXX(xL1,xL2,centre=TRUE,method=c('basis','gq','raw'),control1=list(),
       control2=list(),tol=1e-7)
}
\arguments{
  \item{xL1}{
The mixed scalar and functional variables. For any number and any type of variables, \code{xL1} should be a list. Each item of the list should correspond to one variable. 
}
  \item{xL2}{
Same as \code{xL1}.
}
  \item{centre}{
Logic argument. Default is \code{TRUE}, which means the variables do need to be centred.
}
  \item{method}{
The representative methods for the functional coefficients. The method could be one of the 'basis', 'gq' and 'raw' for basis function expression, Gaussian quadrature and representative data points, respectively. 
}
  \item{control1}{
List of elements that controls the details of the functional coefficients for \code{xL1}. See details for more information. See the argument \code{control} in function \code{\link{fccaGen}} for details.
}
  \item{control2}{
Similar to \code{control1}. 
}
  \item{tol}{
The threshold to decide whether the correlation is to small to be non-zero.   
}
}

\details{
This function uses Moore-Penrose generalized inverse in the calculation to avoid sigular problem.
}

\value{
  \item{corr}{All the non-zero canonical correlation.
}
  \item{coef1}{The corresponding coefficients (weights) for the \code{xL1}.
}
  \item{coef2}{The corresponding coefficients (weights) for the \code{xL2}.
}
}

\examples{
# library(flars)
# library(fda)
## Generate some data sets.
# dataL1=data_generation(seed = 1,uncorr = FALSE,nVar = 8,nsamples = 120,
#       var_type = 'm',cor_type = 1)
# dataL1=dataL1$x

# dataL2=data_generation(seed = 2,uncorr = FALSE,nVar = 8,nsamples = 120,
#       var_type = 'm',cor_type = 1)
# dataL2=dataL2$x

## cross validation
# outCV=fccaXXcv(xL1 = dataL1[1:2], xL2 = dataL2[1:2], method = 'basis'
#      ,alpha = 10^seq(-6,0,len=5))

# cvCor=outCV$cor  
# calculate the correlation
# out=fccaXX(dataL1, dataL2, method = 'basis',control1 = list(pen1=
#      outCV$alpha[which.max(cvCor)]),control2 = list(pen1=
#      outCV$alpha[which.max(cvCor)]))

}
