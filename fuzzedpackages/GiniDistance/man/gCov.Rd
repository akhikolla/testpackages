\name{gCov}
\title{ Gini Distance Covariance Statistics}
\alias{gCov}
\description{
  Computes Gini distance covariance statistics,
  in which Xs are quantitative, Y are categorical, alpha is an exponent on Euclidean distance and returns the measures of dependence.
}
\usage{
  gCov(x, y, alpha)
}
\arguments{
  \item{x}{ data }
  \item{y}{ label of data or univariate response variable}
  \item{alpha}{ exponent on Euclidean distance, in (0,2]}
}
\details{
  \code{gCov} compute Gini distance covariance statistics.
  It is a self-contained R function returning a measure of dependence statistics. 
  
  The sample size (number of rows) of the data must agree with the length of the label vector, and samples must not contain missing values. Arguments
  \code{x}, \code{y} are treated as data and labels. \code{alpha} if missing by default is 1, otherwise it is exponent on the Euclidean distance.
  
  Gini distance covariance is a new measure of dependence between random vectors and its labels. 
  For all distributions with finite first moments, Gini distance
  correlation gCov has the following fundamental properties:
    
  (1) gCov(X,Y) is defined for \eqn{X} in arbitrary dimension quantitive variable and \eqn{Y} a univariate categorical variable.
  
  (2) gCov(X,Y)=0 characterizes independence of \eqn{X} and
  \eqn{Y}.
    
  Gini distance covariance satisfies \eqn{0 \le gCov(X,Y)}, and
  \eqn{gCov = 0} only if \eqn{X} and \eqn{Y} are independent. Gini distance
  covariance gCov provides a new approach to the problem of
  testing the joint independence of random vectors. The formal
  definitions of the population coefficients gCov is given in (DNCZ 2018). The empirical Gini distance covariance \eqn{gCov_n(X,Y; alpha)}  is the nonnegative number computed as follows.
      
  Suppose a sample data \eqn{ {\mathcal{D}} =\{(\mathbf{x}_i,y_i)\} } for \eqn{i = 1,...,n} available. The sample counterparts can be easily computed. Let \eqn{{\mathcal{I}}_k } be the index set of sample points with \eqn{y_i =L_k}, then \eqn{p_k} is estimated by the sample proportion of that category, that is, \eqn{\hat{p}_k= \frac{n_k}{n}} where \eqn{n_k} is the number of elements in \eqn{{\mathcal{I}}_k}. With a given \eqn{\alpha \in (0,2)}, a point estimator of \eqn{\rho_g(\alpha)} is given as follows.  
  \deqn{\hat{\Delta}_k(\alpha)= {n_k \choose 2}^{-1} \sum_{i<j \in {\mathcal{I}}_k} \|\mathbf{x}_i -\mathbf{x}_j\| ^{\alpha},} 
  \deqn{\hat{\Delta}(\alpha)={n \choose 2}^{-1} \sum_{1=i<j=n} \|\mathbf{x}_i -\mathbf{x}_j\| ^{\alpha},} 
  \deqn{{gCov}= \hat{\Delta}(\alpha)-\sum_{k=1}^K \hat p_k \hat{\Delta}_k(\alpha).} 
      
}
\value{
  \code{gCov} returns the sample Gini distance covariance
}
\seealso{
  \code{\link{gCor}} \code{\link{gmd}}  \code{\link{KgCov}}  \code{\link{KgCor}}  
}

\references{
  Dang, X., Nguyen, D., Chen, Y. and Zhang, J., (2019). A new Gini correlation between quantitative and qualitative variables, 
  \emph{Journal of the American Statistical Association (submitted)},
  \url{https://arxiv.org/pdf/1809.09793.pdf}
}
  
\examples{
  x <- iris[,1:4]
  y <- unclass(iris[,5])
  gCov(x, y, alpha = 1) 
}
\keyword{ multivariate }
\concept{ independence }
\concept{ Gini distance correlation }
\concept{ Gini distance covariance }
\concept{ Gini statistics }
