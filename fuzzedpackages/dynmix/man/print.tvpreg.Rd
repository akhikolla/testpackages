\name{print.tvpreg}
\alias{print.tvpreg}

\title{Prints \code{tvpreg} Object.}

\description{The function prints selected outcomes obtained from object \code{tvpreg}.}

\usage{\method{print}{tvpreg}(x, ...)}

\arguments{
\item{x}{an object of \code{tvpreg} class}
\item{...}{not used}
}

\details{
The function prints the general structure of the model, i.e., names of predictors. It also prints the number of observations (length of time-series) and the regression coefficients as estimated in the last period.
}

\examples{
data(oil)
t1<- tvp.reg(y=oil[,1,drop=FALSE],x=oil[,-1,drop=FALSE],lambda=0.99,V=100,W=100)
print(t1)
}
