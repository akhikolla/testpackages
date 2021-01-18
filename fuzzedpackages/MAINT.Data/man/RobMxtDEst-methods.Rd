\name{RobMxtDEst-methods}
\docType{methods}
\alias{RobMxtDEst}
\alias{RobMxtDEst-methods}
\alias{RobMxtDEst,IData-method}

\title{Methods for Function RobMxtDEst in Package `MAINT.Data'}
\description{RobMxtDEst estimates mixtures of distribution for interval-valued data using robust methods.}

\usage{

\S4method{RobMxtDEst}{IData}(Idt, grouping, Mxt=c("Hom","Het"), CovEstMet=c("Pooled","Globdev"),
    CovCase=1:4, SelCrit=c("BIC","AIC"), Robcontrol=RobEstControl(),
    l1medpar=NULL, \dots)
}

\arguments{
  \item{Idt}{An IData object representing interval-valued entities.}
  \item{grouping}{Factor indicating the group to which each observation belongs to.}
  \item{Mxt}{Indicates the type of mixing distributions to be considered. Current alternatives are \dQuote{Hom} (homocedastic) and \dQuote{Het} (hetereocedasic).}
  \item{CovEstMet}{Method used to estimate the common covariance matrix. Alternatives are \dQuote{Pooled} (default) for a pooled average of the the robust within-groups covariance estimates, and \dQuote{Globdev} for a global estimate based on all deviations from the groups multivariate l_1 medians. \cr
See Todorov and Filzmoser (2009) for details..}
  \item{CovCase}{Configuration of the variance-covariance matrix: a set of integers between 1 and 4.}
  \item{SelCrit}{The model selection criterion.}
  \item{Robcontrol}{A control object (S4) of class \code{\link{RobEstControl-class}} containing estimation options - same as these provided in the function specification. If the control object is supplied, the parameters from it will be used. If parameters are passed also in the invocation statement, they will override the corresponding elements of the control object.}
  \item{l1medpar}{List of named arguments to be passed to the function \code{pcaPP::l1median} (in package pcaPP) used to find the multivariate l_1 medians. Possible components are \sQuote{MaxStep}, \sQuote{ItTol} and \sQuote{trace} (see the documentation of \code{pcaPP::l1median} for details). \cr
If kept at NULL (default) the defaults of \code{pcaPP::l1median} will be used.}
  \item{\dots}{Other named arguments.}
}

\value{An object of class IdtMxNDRE, containing the estimation results.}

\keyword{methods}
\keyword{RobMxtDEst}
\keyword{interval data}

\references{
Brito, P., Duarte Silva, A. P. (2012), Modelling Interval Data with Normal and Skew-Normal Distributions. \emph{Journal of Applied Statistics} \bold{39}(1), 3--20.\cr

Duarte Silva, A.P., Filzmoser, P. and Brito, P. (2017), Outlier detection in interval data. \emph{Advances in Data Analysis and Classification}, 1--38.\cr

Todorov V. and Filzmoser P. (2009), An Object Oriented Framework for Robust Multivariate Analysis. \emph{Journal of Statistical Software} \bold{32}(3), 1--47.
}

\seealso{ \code{\linkS4class{IdtMxNDRE}}, \code{\link{RobEstControl}}.}

\examples{
# Create an Interval-Data object containing the intervals for 899 observations 
# on the temperatures by quarter in 60 Chinese meteorological stations.

ChinaT <- IData(ChinaTemp[1:8],VarNames=c("T1","T2","T3","T4"))

\dontrun{

# Estimate robustly an homoscedastic mixture, with mixture components defined by regions

ChinaHomMxtRobE <- RobMxtDEst(ChinaT,ChinaTemp$GeoReg)
print(ChinaHomMxtRobE)

# Estimate robustly an heteroscedastic mixture, with mixture components defined by regions

ChinaHetMxtRobE <- RobMxtDEst(ChinaT,ChinaTemp$GeoReg,Mxt="Het")
print(ChinaHetMxtRobE)


}

}

