\name{gy.filt}
\alias{gy.filt}
\title{Gervini-Yohai filter for detecting cellwise outliers}
\description{
  Flags cellwise outliers detected using Gervini-Yohai filter as described in Agostinelli et al. (2015) and Leung and Zamar (2016). 
}
\usage{
gy.filt(x, alpha=c(0.95,0.85), bivarQt=0.99, bivarCellPr=0.1, miter=5)
}
\arguments{
  \item{x}{a matrix or data frame.}
  \item{alpha}{a vector of the quantiles of the univariate and bivariate reference distributions, respectively. 
  Filtering turns off when alpha is 0. For univariate filtering only, \code{alpha=c(0.95,0)}.
  Default value is \code{c(0.95,0.85)}.}
  \item{bivarQt}{quantile of the binomial model for the number of flagged pairs in the bivariate filter. Default is 0.99.}
  \item{bivarCellPr}{probability of the binomial model for the number of flagged pairs in the bivariate filter. Default is 0.1.}
  \item{miter}{maximum number of iteration of filtering. Default value is 5.}
}
\details{

   This function implements the univariate filter and the univariate-plus-bivariate filter as described in Agostinelli et al. (2015)
   and Leung and Zamar (2016), respectively. 

   In the univariate filter, outliers are flagged by comparing the empirical tail distribution of each marginal 
   with a reference (normal) distribution using Gervini-Yohai approach. 
   
   In the univiarate-plus-bivariate filter, outliers are first flagged by applying the univariate filter. Then, the bivariate filter
   is applied to flag any additional outliers. In the bivariate filter, outliers are flagged by comparing the empirical tail distribution 
   of each bivariate marginal with a reference (chi-square with 2 d.f.) distribution using Gervini-Yohai approach. The number of 
   flagged pairs associated with each cell approximately follows a binomial model under Independent Cellwise Contamination Model. 
   A cell is additionally flagged if the number of flagged pairs exceeds a large quantile of the binomial model. 

}
\value{
    a matrix or data frame of the filtered data.
}
\references{
    Agostinelli, C., Leung, A. , Yohai, V.J., and Zamar, R.H. (2015) Robust estimation of multivariate location and scatter in the presence of cellwise and casewise contamination. TEST.
    
    Leung, A. and Zamar, R.H. (2016). Multivariate Location and Scatter Matrix Estimation Under Cellwise and Casewise Contamination.
    Submitted.        
}
\author{Andy Leung \email{andy.leung@stat.ubc.ca}, Claudio Agostinelli, Ruben H. Zamar, Victor J. Yohai}
\seealso{
  \code{\link{TSGS}}, \code{\link{generate.cellcontam}}
}
\examples{
set.seed(12345)

# Generate 5% cell-wise contaminated normal data 
x <- generate.cellcontam(n=100, p=10, cond=100, contam.size=5, contam.prop=0.05)$x

## Using univariate filter only
xna <- gy.filt(x, alpha=c(0.95,0))
mean(is.na(xna))

## Using univariate-and-bivariate filter
xna <- gy.filt(x, alpha=c(0.95,0.95))
mean(is.na(xna))
}


