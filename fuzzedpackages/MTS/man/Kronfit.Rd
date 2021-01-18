\name{Kronfit}
\alias{Kronfit}
\title{Fitting a VARMA Model via Kronecker Index
}
\description{Perform estimation of a VARMA model specified by the Kronecker 
indices
}
\usage{
Kronfit(da, kidx, include.mean = T, fixed = NULL, Kpar=NULL, 
  seKpar=NULL, prelim = F, details = F, thres = 1)
}
\arguments{
  \item{da}{Data matrix (T-by-k) of a k-dimensional time series
}
  \item{kidx}{The vector consisting of Kronecker indices
}
  \item{include.mean}{A logical switch for including the mean vector in 
estimation. Default is to include the mean vector.
}
  \item{fixed}{A logical matrix used to set zero parameter constraints. 
This is used mainly in the command refKronfit.
}
  \item{Kpar}{Parameter vectors for use in model simplification
}
  \item{seKpar}{Standard errors of the parameter estimates for use in model simplification
}
  \item{prelim}{A logical switch for a preliminary estimation.
}
  \item{details}{A logical switch to control output.
}
  \item{thres}{A threshold for t-ratios in setting parameter to zero. 
Default is 1.
}
}
\value{
\item{data}{The observed time series data}
\item{Kindex}{Kronecker indices}
\item{ARid}{Specification of AR parameters: 0 denotes fixing to zero, 1 denotes fixing to 1, and 2 denoting estimation}
\item{MAid}{Specification of MA parameters}
\item{cnst}{A logical variable: include.mean}
\item{coef}{Parameter estimates}
\item{se.coef}{Standard errors of the estimates}
\item{residuals}{Residual series}
\item{Sigma}{Residual covariance matrix}
\item{aic,bic}{Information criteria of the fitted model}
\item{Ph0}{Constant vector}
\item{Phi}{AR coefficient matrices}
\item{Theta}{MA coefficient matrices}
}
\references{Tsay (2014, Chapter 4). Multivariate Time Series Analysis with 
R and Financial Applications. John Wiley. Hoboken, NJ.
}
\author{Ruey S. Tsay
}
\seealso{refKronfit, Kronspec
}



