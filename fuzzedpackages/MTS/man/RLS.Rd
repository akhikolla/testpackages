\name{RLS}
\alias{RLS}
\title{Recursive Least Squares
}
\description{Compute recursive least squares estimation
}
\usage{
RLS(y, x, ist = 30, xpxi = NULL, xpy0 = NULL)
}
\arguments{
  \item{y}{data of dependent variable
}
  \item{x}{data matrix of regressors
}
  \item{ist}{initial number of data points used to start the estimation
}
  \item{xpxi}{Inverse of the X'X matrix
}
  \item{xpy0}{Initial value of X'y.
}
}
\value{
\item{beta}{Time-varying regression coefficient estimates}
\item{resi}{The residual series of recursive least squares estimation} 
}
\author{Ruey S. Tsay
}
\note{This function is used internally, but can also be used as a command.
}




