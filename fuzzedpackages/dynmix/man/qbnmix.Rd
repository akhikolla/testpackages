\name{qbnmix}
\alias{qbnmix}

\title{Estimates Normal Mixtures.}

\description{
This function performs a recursive clustering for normal mixtures. Quasi-Bayesian approximation is performed.
}

\usage{
qbnmix(y,m=2,mu0=NULL,R0=NULL)
}

\arguments{
\item{y}{{\code{\link[base]{matrix}}} of observations, rows correspond to observations, columns correspond to tuples}
\item{m}{{\code{\link[base]{numeric}}} specifying the number of components (clusters), by default \code{m=2} is taken}
\item{mu0}{optional, initial means, should be a \code{\link[base]{list}} of \code{m} matrices, each of them having one row and \code{ncol(y)} columns, if not specified random values are taken}
\item{R0}{optional, initial covariance matrices, should be a \code{\link[base]{list}} of \code{m} matrices, each of them having \code{ncol(y)} rows and \code{ncol(y)} columns, if not specified identity matrices are taken}
}

\value{object of class {\code{qbnmix}}, i.e., {\code{\link[base]{list}}} of 
\item{$mu}{{\code{\link[base]{list}}} of estimated means}
\item{$R}{{\code{\link[base]{list}}} of estimated covariance matrices (from last step only)}
\item{$alpha}{{\code{\link[base]{matrix}}} of estimates of mixing weights (components columnwise)}
\item{$w}{{\code{\link[base]{matrix}}} of posterior probabilities (components columnwise)}
\item{$mu0}{{\code{\link[base]{list}}} of initial means matrices}
\item{$R0}{{\code{\link[base]{list}}} of initial covaraince matrices}
}

\examples{
R <- list(matrix(c(1,0.3,0,
                   0.3,0.3,0,
                   0,0,0.15),3,3),
          matrix(c(1,0,0,
                   0,0.5,0,
                   0,0,0.2),3,3))
data <- rbind(MASS::mvrnorm(n=180,c(5,2,3),R[[1]]),
              MASS::mvrnorm(n=20,c(1,2,3),R[[2]]))
data <- data[sample(nrow(data)),]
                  
mu0 <- list(matrix(c(4.8689,1.9417,3.0175),nrow=1,ncol=3),
            matrix(c(1.0182,1.9903,2.8847),nrow=1,ncol=3))     
est <- qbnmix(y=data,mu0=mu0)
}

\source{
Karny, M., Kadlec, J., Sutanto, E.L., 1998, Quasi-Bayes estimation applied to normal mixture,
\emph{Preprints of The 3rd European IEEE Workshop on Computer-Intensive Methods in Control and Data Processing}, Rojicek, J., Valeckova, M., Karny, M., Warwick K. (eds.), UTIA AV CR, 77--82.   
}
