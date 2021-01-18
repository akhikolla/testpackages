\name{Mlm}
\alias{Mlm}
\title{Multivariate Linear Model
}
\description{Fit a multivariate multiple linear regression model via the least 
squares method
}
\usage{
Mlm(y, z, constant=TRUE, output=TRUE)
}
\arguments{
  \item{y}{data matrix of dependent variable. Each column contains one 
variable.
}
  \item{z}{data matrix of the explanatory variables. Each column contains one 
variable.
}
  \item{constant}{A logical switch for including the constant term
}
  \item{output}{A logical switch to print the output
}
}
\value{
\item{beta}{coefficient matrix}
\item{se.beta}{standard errors of the coefficient matrix}
\item{residuals}{The residual series} 
\item{sigma}{Residual covariance matrix}
}
\author{Ruey S. Tsay
}





