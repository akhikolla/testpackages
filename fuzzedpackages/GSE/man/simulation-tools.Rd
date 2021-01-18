\name{simulation-tools}
\alias{generate.randcorr}
\alias{generate.cellcontam}
\alias{generate.casecontam}
\title{Data generator for simulation study on cell- and case-wise contamination}
\description{
	Includes the data generator for the simulation study on cell- and case-wise contamination that appears on Agostinelli et al. (2014). 
}
\usage{
generate.randcorr(cond, p, tol=1e-5, maxits=100) 

generate.cellcontam(n, p, cond, contam.size, contam.prop, A=NULL)

generate.casecontam(n, p, cond, contam.size, contam.prop, A=NULL)
}
\arguments{
  \item{cond}{desired condition number of the random correlation matrix. The correlation matrix will be used to generate multivariate normal samples in \code{generate.cellcontam} and \code{generate.cellcontam}.}
  \item{tol}{tolerance level for the condition number of the random correlation matrix. Default is \code{1e-5}.}
  \item{maxits}{integer indicating the maximum number of iterations until the condition number of the random correlation matrix is within a tolerance level. Default is 100.}
  \item{n}{integer indicating the number of observations to be generated.}
  \item{p}{integer indicating the number of variables to be generated.}
  \item{contam.size}{size of cell- or case-wise contamination. For cell-wise outliers, random cells in a data matrix are replaced by \code{contam.dist}. 
  For case-wise outliers, random cases in a data matrix are replaced by \code{contam.dist} times \eqn{v} where \eqn{v}}
  \item{contam.prop}{proportion of cell- or case-wise contamination.}
  \item{A}{correlation matrix used for generating data. If \code{A} is \code{NULL}, 
  a random correlation matrix is generated. Default is \code{NULL}.}
}
\details{
	Details about how the correlation matrix is randomly generated and how the contaminated data is generated can be found in Agostinelli et al. (2014). 
}
\value{
	\code{generate.randcorr} gives the random correlation matrix in dimension \code{p} and with condition number \code{cond}. 
	
	\code{generate.cellcontam} and \code{generate.casecontam} give the multivariate normal sample that is either cell-wise 
	or case-wise contaminated as described in Agostinelli et al. (2014). The contaminated sample is returned as components of a list with components
	\tabular{ll}{
	\code{x} \tab multivariate normal sample with cell- or case-wise contamination. \cr
	\code{u} \tab \code{n} by \code{p} matrix of 0's and 1's with 1's correspond to an outlier. A row of 1's correspond to a case-wise outlier.\cr 
	\code{A} \tab random correlation matrix with a specified condition number. \cr
	}
}
\references{
 	Agostinelli, C., Leung, A. , Yohai, V.J., and Zamar, R.H. (2014) Robust estimation of multivariate location and scatter in the presence of cellwise and casewise contamination. arXiv:1406.6031[math.ST]
}
\author{Andy Leung \email{andy.leung@stat.ubc.ca}, Claudio Agostinelli, Ruben H. Zamar, Victor J. Yohai}
\seealso{
  \code{\link{TSGS}}
}



