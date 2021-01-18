\name{dccPre}
\alias{dccPre}
\title{Preliminary Fitting of DCC Models
}
\description{This program fits marginal GARCH models to each component 
of a vector return series and returns the standardized return series for 
further analysis. The garchFit command of fGarch package is used.
}
\usage{
dccPre(rtn, include.mean = T, p = 0, cond.dist = "norm")
}
\arguments{
  \item{rtn}{A T-by-k data matrix of k-dimensional asset returns
}
  \item{include.mean}{A logical switch to include a mean vector. 
Deafult is to include the mean.
}
  \item{p}{VAR order for the mean equation
}
  \item{cond.dist}{The conditional distribution of the innovations.
Default is Gaussian. 
}
}
\details{The program uses fGarch package to estimate univariate GARCH 
model for each residual series after a VAR(p) fitting, if any.
}
\value{
\item{marVol}{A matrix of the volatility series for each return series}
\item{sresi}{Standardized residual series}
\item{est}{Parameter estimates for each marginal volatility model}
\item{se.est}{Standard errors for parameter estimates of marginal volatility models}
}
\references{Tsay (2014, Chapter 7). Multivariate Time Series Analysis with 
R and Financial Applications. John Wiley. Hoboken, NJ.
}
\author{Ruey S. Tsay
}
\note{fGarch package is used
}
\seealso{dccFit
}




