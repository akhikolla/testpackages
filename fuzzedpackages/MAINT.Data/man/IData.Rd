\name{IData}
\alias{IData}
\alias{is.IData}
\title{Interval Data objects}
\description{
IData creates IData objects from data frames of interval bounds or MidPoint/LogRange values of the interval-valued observations.
}

\usage{
IData(Data, 
Seq = c("LbUb_VarbyVar", "MidPLogR_VarbyVar", "AllLb_AllUb", "AllMidP_AllLogR"), 
VarNames=NULL, ObsNames=row.names(Data))
}

\arguments{
  \item{Data}{a data frame or matrix of interval bounds or MidPoint/LogRange values.}
  \item{Seq}{the format of \sQuote{Data} data frame. Available options are:\cr
\dQuote{LbUb_VarbyVar}: lower bounds followed by upper bounds, variable by variable.\cr
\dQuote{MidPLogR_VarbyVar}: MidPoints followed by LogRanges, variable by variable.\cr
\dQuote{AllLb_AllUb}: all lower bounds followed by all upper bounds, in the same variable order.\cr
\dQuote{AllMidP_AllLogR}: all MidPoints followed all LogRanges, in the same variable order.\cr
}
  \item{VarNames}{An optional vector of names to be assigned to the Interval-Valued Variables.}
  \item{ObsNames}{An optional vector of names assigned to the individual observations.}
}

\details{Objects of class \code{\linkS4class{IData}} describe a data set of \sQuote{NObs} observations on \sQuote{NIVar} Interval-valued variables. This function creates an interval-data object from a data-frame with either the lower and upper bounds of the observed intervals or by their midpoints and log-ranges.}

\seealso{
  \code{\linkS4class{IData}}, \code{\link{AgrMcDt}}
}

\keyword{interval data}

\examples{
ChinaT <- IData(ChinaTemp[1:8],VarNames=c("T1","T2","T3","T4"))
cat("Summary of the ChinaT IData object:\n")  ; print(summary(ChinaT))
cat("ChinaT first ant last three observations:\n")  
print(head(ChinaT,n=3))
cat("\n...\n")
print(tail(ChinaT,n=3))
}



