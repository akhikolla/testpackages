\name{EventElement}
\alias{EventElement}

\title{ Element entry as OpLine line item. }

\description{
Provided for completeness of package only.  Oprating line entries are more easily
entered as part of an OpLine system using a spreadsheet, then using RExcel to 
"Put R Value -> dataframe" into an R session.
}

\usage{
EventElement(element_name,OpLine,EventID,FD,FP1,FP2,FP3,RD,RP1,RP2,RP3,Seed) 
}

\arguments{
\item{element_name}{ a string that will used to identify the element as a row.name}
\item{OpLine}{ an identifying integer for the OpLine membership.}
\item{EventID}{ an identifying integer for this element }
\item{FD}{ a single capital letter in quotation marks indicating the failure 
distribution (i.e. "E"-Exponential, "N"-Normal, "W"-Weibull) }
\item{FP1}{ The first parameter for the failure distribution (MTTF for "E",
mean for "N", characteristic life for "W").}
\item{FP2}{ The second parameter for the failure distribution, else 0 }
\item{FP3}{ The third parameter for the failure distribution (only used as a
translation parameter for "W", else 0 )}
\item{RD}{ a single capital letter in quotation marks indicating the repair 
distribution (i.e. "N"-Normal, "W"-Weibull, "L"-Lognormal) 
note that "E" is not accepted, use "W" with shape of 1.0 instead.}
\item{RP1}{ The first parameter for the failure distribution (mean for "N",
characteristic life for "W", logmean for "L")}
\item{RP2}{ The second parameter for the repair distribution }
\item{RP3}{ The third parameter for the failure distribution (only used as a
translation parameter for "W", or "L" else 0 )}
\item{Seed}{ An integer value to be used as the random seed for this element.}
}


\value{
Returns a single row dataframe suitable for combination using \code{rbind} with other elements into a single Operating Line.
}

\references{
Carazas et. al.,"Availability Analysis of Gas Turbines Used in Power Plants",International Journal of Thermodynamics, Vol. 12 (No.1), March 2009
}

\examples{
LRU1 <- EventElement("GasTurbine2",1,101,"W",2562.5,0.95,0, "L",1.4,0.86,0,87)

}

\keyword{ input }

