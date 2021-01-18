\name{dCor}
\title{Distance Covariance and Correlation Statistics}
\alias{dCor}
\description{
  Computes distance covariance and correlation statistics, in which Xs are quantitative and Ys are categorical and return the measures of dependence.  
}
\usage{
  dCor(x, y, alpha)
}
\arguments{
  \item{x}{ data }
  \item{y}{ label of data or univariate response variable}
  \item{alpha}{ exponent on Euclidean distance, in (0,2]}
}
\details{
  
  The sample size (number of rows) of the data must agree with the length of the label vector, and samples must not contain missing values. Arguments
  \code{x}, \code{y} are treated as data and labels.
  
  
  \code{dCor} calls \code{dcor} function from energy package which computes the distance correlation between X and Y where both are numerical variables. If Y is categorical,  the set difference metric on the support of \eqn{Y} is used. That is, \eqn{d(y, y^\prime) =|y-y^\prime|:=  I(y\neq y^\prime),}
  where \eqn{I (\cdot)} is the indicator function. Then the sample distance correlation between data and labels is computed as follows. 
  
  Let \eqn{A=(a_{ij})} be a symmetric, \eqn{n \times n}, centered distance matrix of sample \eqn{\mathbf x_1,\cdots, \mathbf x_n}. The \eqn{(i,j)}-th entry of \eqn{A} is \eqn{a_{ij}-\frac{1}{n-2}a_{i\cdot}-\frac{1}{n-2}a_{\cdot j} + \frac{1}{(n-1)(n-2)}a_{\cdot \cdot}} if \eqn{i \neq j} and 0 if \eqn{i=j},
    where \eqn{a_{ij} = \|\mathbf x_i-\mathbf x_j\|^{\alpha}}, \eqn{a_{i\cdot} = \sum_{j=1}^n a_{ij}}, \eqn{a_{\cdot j} = \sum_{i=1}^n a_{ij}}, and \eqn{a_{\cdot \cdot}=\sum_{i,j=1}^n a_{ij}}. Similarly, using the set difference metric, a symmetric, \eqn{n \times n}, centered distance matrix is calculated for samples \eqn{y_1,\cdots, y_n} and denoted by \eqn{B = (b_{ij})}. Unbiased estimators of \eqn{\mbox{dCov}(\mathbf X,Y;\alpha)}, \eqn{\mbox{dCov}(\mathbf X, \mathbf X;\alpha)} and \eqn{\mbox{dCov}(\mathbf Y, \mathbf Y;\alpha)}  are given  respectively as,  \eqn{\frac{1}{n(n-3)}\sum_{i\ne j}A_{ij}B_{ij}}, \eqn{\frac{1}{n(n-3)}\sum_{i\ne j}A_{ij}^2} and \eqn{\frac{1}{n(n-3)}\sum_{i\ne j}B_{ij}^2}. Then the distance correlation is 
      
      \deqn{{dCor}(\mathbf{X}, Y; \alpha) = \frac{\mbox{ dCov}(\mathbf{X}, Y, \alpha)}{ \sqrt{\mbox{ dCov}(\mathbf{X},\mathbf{X};\alpha)} \sqrt{\mbox{ dCov}(Y,Y)}}.}
    
}
\value{
  \code{dCor} returns the sample distance variance of \code{x}, distance variance of \code{y}, distance covariance of \code{x} and \code{y} and distance correlation of \code{x}, \code{y}. 
}
\seealso{
  \code{\link{dCov}}  \code{\link{KdCov}}  \code{\link{KdCor}} 
}
\references{
  
  Lyons, R. (2013). Distance covariance in metric spaces. The Annals of Probability, 41 (5), 3284-3305. 
  
  Szekely, G. J., Rizzo, M. L. and Bakirov, N. (2007). Measuring and testing dependence by correlation of distances. Annals of Statistics, 35 (6), 2769-2794. 
  
  Rizzo, M.L. and Szekely, G.J., (2017). Energy: E-Statistics: Multivariate Inference via the Energy of Data (R Package), Version 1.7-0.
  
  
}
\examples{
  x <- iris[,1:4]
  y <- unclass(iris[,5])
  dCor(x, y, alpha = 1)
}
\keyword{ multivariate }
\concept{ independence }
\concept{ Gini distance correlation }
\concept{ Gini distance covariance }
\concept{ Gini statistics }
