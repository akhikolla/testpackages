\name{dCov}
\title{Distance Covariance Statistic}
\alias{dCov}
\description{
  Computes distance covariance statistic,
  in which Xs are quantitative and Y are categorical and return the measures of dependence.
}
\usage{
  dCov(x, y, alpha)
}
\arguments{
  \item{x}{ data }
  \item{y}{ label of data or response variable}
  \item{alpha}{ exponent on Euclidean distance, in (0,2]}
}
\details{
  \code{dCov} calls \code{dcov} function from energy package to compute distance covariance statistic.
  The sample size (number of rows) of the data must agree with the length of the label vector, and samples must not contain missing values. Arguments
  \code{x}, \code{y} are treated as data and labels.
  
  The distance covariance (\cite{Sezekley07}) is extended from Euclidean space to general metric spaces by Lyons (2013). Based on that idea, we define the discrete metric   
  \deqn{d(y, y^\prime) =|y-y^\prime|:=  I(y\neq y^\prime),}
  where \eqn{I (\cdot)} is the indicator function. Equipped with this set difference metric on the support of \eqn{Y} and Euclidean 
  distance on the support of \eqn{\mathbf{X}}, the corresponding distance covariance and distance correlation for numerical \eqn{\mathbf{X}} and categorical \eqn{Y} variables are as follows. 
  
  Let \eqn{A=(a_{ij})} be a symmetric, \eqn{n \times n}, centered distance matrix of sample \eqn{\bf x_1,\cdots, \bf x_n}. The \eqn{(i,j)}-th entry of \eqn{A} is \eqn{a_{ij}-\frac{1}{n-2}a_{i\cdot}-\frac{1}{n-2}a_{\cdot j} + \frac{1}{(n-1)(n-2)}a_{\cdot \cdot}} if \eqn{i \neq j} and 0 if \eqn{i=j},
    where \eqn{a_{ij} = \|\bf x_i-\bf x_j\|^{\alpha}}, \eqn{a_{i\cdot} = \sum_{j=1}^n a_{ij}}, \eqn{a_{\cdot j} = \sum_{i=1}^n a_{ij}}, and \eqn{a_{\cdot \cdot}=\sum_{i,j=1}^n a_{ij}}. Similarly, using the set difference metric, a symmetric, \eqn{n \times n}, centered distance matrix is calculated for samples \eqn{y_1,\cdots, y_n} and denoted by \eqn{B = (b_{ij})}. Unbiased estimators of \eqn{\mbox{dCov}(\bf X,Y;\alpha)} is   
      
      \eqn{\frac{1}{n(n-3)}\sum_{i\ne j}A_{ij}B_{ij}}.
}    
\value{
  \code{dCov} returns the sample distance covariance between data \code{x} and label \code{y}. 
}
\seealso{
  \code{\link{dCor}} \code{\link{KdCov}}  \code{\link{KdCor}}
}
\references{
  
  Lyons, R. (2013). Distance covariance in metric spaces. The Annals of Probability, 41 (5), 3284-3305. 
  
  
  Rizzo, M.L. and Szekely, G.J., (2017). Energy: E-Statistics: Multivariate Inference via the Energy of Data (R Package), Version 1.7-0.
  
  Szekely, G. J., Rizzo, M. L. and Bakirov, N. (2007). Measuring and testing dependence by correlation of distances. Annals of Statistics, 35 (6), 2769-2794.
}
\examples{
  x <- iris[,1:4]
  y <- unclass(iris[,5])
  dCov(x, y, alpha = 1)
}
\keyword{ multivariate }
\concept{ independence }
\concept{ distance correlation }
\concept{ distance covariance }
\concept{ energy statistics }
