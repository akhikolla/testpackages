\name{apca}
\alias{apca}
\title{Asymptotic Principal Component Analysis
}
\description{Perform asymptotic PCA for a data set. Typically for cases in which 
the number of variables is greater than the number of data points.
}
\usage{
apca(da, m)
}
\arguments{
  \item{da}{A T-by-k data set matrix, where T is the sample size and k is the 
dimension
}
  \item{m}{The number of common factors
}
}
\details{Perform the PCA analysis of interchanging the roles of variables and observations.
}
\value{
\item{sdev}{Square root of the eigenvalues}
\item{factors}{The common factors}
\item{loadings}{The loading matrix}
}
\references{Tsay (2014, Chapter 6). Multivariate Time Series Analysis with R and Financial Applications. John Wiley. Hoboken, NJ.
}
\author{Ruey S. Tsay
}
\examples{
rtn=matrix(rnorm(1200),12,100)
sp100=apca(rtn,3)
}




