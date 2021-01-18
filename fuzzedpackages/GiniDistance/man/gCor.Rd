\name{gCor}
\title{Gini Distance Covariance and Correlation Statistics}
\alias{gCor}
\description{
  Computes Gini distance covariance and correlation statistics,
  in which Xs are quantitative, Y are categorical, alpha is exponent on the Euclidean distance and returns the measures of dependence.
}
\usage{
  gCor(x, y, alpha)
}
\arguments{
  \item{x}{ data }
  \item{y}{ label of data or univariate response variable}
  \item{alpha}{ exponent on Euclidean distance, in (0,2)}
}
\details{
  \code{gCor} compute Gini distance correlation statistics.
  It is a self-contained R function returning a measure of dependence statistics. 
  
  The sample size (number of rows) of the data must agree with the length of the label vector, and samples must not contain missing values. Arguments
  \code{x}, \code{y} are treated as data and labels. \code{alpha} if missing by default is 1, otherwise it is exponent on the Euclidean distance.
  
  Suppose a sample data \eqn{ {\mathcal{D}} =\{(\mathbf{x}_i,y_i)\} } for \eqn{i = 1,...,n} available. The sample counterparts can be easily computed. Let \eqn{{\mathcal{I}}_k } be the index set of sample points with \eqn{y_i =L_k}, then  \eqn{p_k} is estimated by the sample proportion of that category, that is, \eqn{\hat{p}_k= \frac{n_k}{n}} where \eqn{n_k} is the number of elements in \eqn{{\mathcal{I}}_k}. With a given \eqn{\alpha \in (0,2)}, a point estimator of \eqn{\rho_g(\alpha)} is given as follows.  
  \deqn{\hat{\Delta}_k(\alpha)= {n_k \choose 2}^{-1} \sum_{i<j \in {\mathcal{I}}_k} \|\mathbf{x}_i -\mathbf{x}_j\| ^{\alpha},} 
  \deqn{\hat{\Delta}(\alpha)={n \choose 2}^{-1} \sum_{1=i<j=n} \|\mathbf{x}_i -\mathbf{x}_j\| ^{\alpha},} 
  \deqn{gCor=\hat{\rho}_g (\alpha)= 1-\frac{\sum_{k=1}^K \hat p_k \hat{\Delta}_k(\alpha)}{\hat{\Delta}(\alpha)}.} 
}
\value{
  \code{gCor} returns the sample Gini distance covariacne and correlation between \code{x} and \code{y}.
}

\references{
  Dang, X., Nguyen, D., Chen, Y. and Zhang, J. (2019). A new Gini correlation between quantitative and qualitative variables. Submitted to Journal of American Statistics Association.   

}

\seealso{
  \code{\link{gmd}} \code{\link{gCov}} \code{\link{KgCov}}  \code{\link{KgCor}}
} 

\examples{
  x <- iris[,1:4]
  y <- unclass(iris[,5])
  gCor(x, y, alpha = 1)
}

\keyword{ multivariate }
\concept{ independence }
\concept{ Gini distance correlation }
\concept{ Gini distance covariance }
\concept{ Gini statistics }
