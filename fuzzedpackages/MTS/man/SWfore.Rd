\name{SWfore}
\alias{SWfore}
\title{Stock-Watson Diffusion Index Forecasts
}
\description{Uses the diffusion index approach of Stock and Watson to 
compute out-of-sample forecasts
}
\usage{
SWfore(y, x, orig, m)
}
\arguments{
  \item{y}{The scalar variable of interest 
}
  \item{x}{The data matrix (T-by-k) of the observed explanatory variables
}
  \item{orig}{Forecast origin
}
  \item{m}{The number of diffusion index used
}
}
\details{Performs PCA on X at the forecast origin. Then, fit a linear 
regression model to obtain the coefficients of prediction equation. 
Use the prediction equation to produce forecasts and compute forecast 
errors, if any. No recursive estimation is used.
}
\value{
\item{coef}{Regression coefficients of the prediction equation}
\item{yhat}{Predictions at the forecast origin}
\item{MSE}{Mean squared errors, if available}
\item{loadings}{Loading matrix}
\item{DFindex}{Diffusion indices}
}
\references{Tsay (2014, Chapter 6). Multivariate Time Series Analysis with R and 
Financial Applications. John Wiley. Hoboken, NJ.
}
\author{Ruey S. Tsay
}




