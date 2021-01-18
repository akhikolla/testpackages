\name{Robda-methods}
\docType{methods}
\alias{Roblda}
\alias{Roblda-methods}
\alias{Roblda,IData-method}
\alias{Robqda}
\alias{Robqda-methods}
\alias{Robqda,IData-method}

\title{Robust Discriminant Analysis of Interval Data}

\description{Roblda and Robqda perform linear and quadratic discriminant analysis of Interval Data based on robust estimates of location and scatter.}

\usage{

\S4method{Roblda}{IData}( x, grouping, prior="proportions", CVtol=1.0e-5, egvtol=1.0e-10,
  subset=1:nrow(x), CovCase=1:4, SelCrit=c("BIC","AIC"), silent=FALSE, 
  CovEstMet=c("Pooled","Globdev"), SngDMet=c("fasttle","fulltle"), k2max=1e6,
  Robcontrol=RobEstControl(), \dots )

\S4method{Robqda}{IData}( x, grouping, prior="proportions", CVtol=1.0e-5, 
  subset=1:nrow(x), CovCase=1:4, SelCrit=c("BIC","AIC"), silent=FALSE,
  SngDMet=c("fasttle","fulltle"), k2max=1e6, Robcontrol=RobEstControl(), \dots )

}

\arguments{

  \item{x}{An object of class \code{\linkS4class{IData}} with the original Interval Data.}

  \item{grouping}{Factor specifying the class for each observation.}

  \item{prior}{The prior probabilities of class membership.  If unspecified, the class proportions for the training set are used.  If
present, the probabilities should be specified in the order of the factor levels.}

  \item{CVtol}{Tolerance level for absolute value of the coefficient of variation of non-constant variables. When a MidPoint or LogRange has an absolute value within-groups coefficient of variation below CVtol, it is considered to be a constant.}

  \item{egvtol}{Tolerance level for the eigenvalues of the product of the inverse within by the between covariance matrices. When a  eigenvalue has an absolute value below egvtol, it is considered to be zero.}

  \item{subset}{An index vector specifying the cases to be used in the analysis.}

  \item{CovCase}{Configuration of the variance-covariance matrix: a set of integers between 1 and 4.}

  \item{SelCrit}{The model selection criterion.}

  \item{silent}{A boolean flag indicating wether a warning message should be printed if the method fails.}
 
  \item{CovEstMet}{Method used to estimate the common covariance matrix in \code{Roblda} (Robust linear discriminant analysis). Alternatives are \dQuote{Pooled} (default) for a pooled average of the the robust within-groups covariance estimates, and \dQuote{Globdev} for a global estimate based on all deviations from the groups multivariate l_1 medians. See Todorov and Filzmoser (2009) for details.}
 
  \item{SngDMet}{Algorithm used to find the robust estimates of location and scatter. Alternatives are \dQuote{fasttle} (default) and \dQuote{fulltle}.}
 
  \item{k2max}{Maximal allowed l2-norm condition number for correlation matrices. Correlation matrices with condition number above k2max are considered to be numerically singular, leading to degenerate results.}

  \item{Robcontrol}{A control object (S4) of class \code{\link{RobEstControl-class}}
    containing estimation options - same as these provided in the function
    specification. If the control object is supplied, the parameters from it
    will be used. If parameters are passed also in the invocation statement, they will
    override the corresponding elements of the control object.}
 
  \item{\dots}{Other named arguments.}
}

\references{
Duarte Silva, A.P. and Brito, P. (2015), Discriminant analysis of interval data: An assessment of parametric and distance-based approaches. \emph{Journal of Classification} \bold{39}(3), 516--541.\cr

Duarte Silva, A.P., Filzmoser, P. and Brito, P. (2017), Outlier detection in interval data. \emph{Advances in Data Analysis and Classification}, 1--38.\cr
}

\seealso{\code{\link{lda}}, \code{\link{qda}}, \code{\link{snda}}, \code{\linkS4class{IData}}, \code{\linkS4class{RobEstControl}}}

\keyword{methods}
\keyword{Interval Data}
\keyword{Discriminant Analysis Interval Data}

\examples{
# Create an Interval-Data object containing the intervals for 899 observations 
# on the temperatures by quarter in 60 Chinese meteorological stations.

ChinaT <- IData(ChinaTemp[1:8],VarNames=c("T1","T2","T3","T4"))

#Robust Linear Discriminant Analysis

\dontrun{

ChinaT.rlda <- Roblda(ChinaT,ChinaTemp$GeoReg)
cat("Temperatures of China -- robust lda discriminant analysis results:\n")
print(ChinaT.rlda)

#Estimate error rates by ten-fold cross-validation with 5 replications 

CVrlda <- DACrossVal(ChinaT,ChinaTemp$GeoReg,TrainAlg=Roblda,CovCase=CovCase(ChinaT.rlda),
   CVrep=5)
summary(CVrlda[,,"Clerr"])

#Robust Quadratic Discriminant Analysis

ChinaT.rqda <- Robqda(ChinaT,ChinaTemp$GeoReg)
cat("Temperatures of China -- robust qda discriminant analysis results:\n")
print(ChinaT.rqda)

#Estimate error rates by ten-fold cross-validation with 5 replications 

CVrqda <- DACrossVal(ChinaT,ChinaTemp$GeoReg,TrainAlg=Robqda,CovCase=CovCase(ChinaT.rqda),
   CVrep=5)
summary(CVrqda[,,"Clerr"])

}

}
