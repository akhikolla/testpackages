\name{oil}
\alias{oil}

\title{Crude Oil Data.}

\description{Selected data from oil market.}

\usage{data(oil)}

\format{\code{oil} is \code{\link[base]{matrix}} object such that columnwise are
\itemize{
\item WTI -- WTI spot price in USD per barrel 
\item MSCI -- MSCI World Index 
\item TB3MS -- U.S. 3-month treasury bill secondary market rate 
\item TWEXM -- Trade weighted U.S. dollar index (Mar, 1973 = 100) 
\item PROD -- U.S. product supplied for crude oil and petroleum products in thousands of barrels 
}
}

\details{The data are in monthly frequency. They cover the period between Feb, 1990 and Dec, 2016. MSCI, TB3MS, TWEXM and PROD are lagged one period back.
}

\source{The data are provided by Federal Reserve Bank of St. Louis, MSCI and U.S. Energy Information Administration. 

\url{https://www.eia.gov}

\url{https://fred.stlouisfed.org}

\url{https://www.msci.com/end-of-day-data-search}
}

\examples{
data(oil)
}
