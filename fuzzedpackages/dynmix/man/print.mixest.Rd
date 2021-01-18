\name{print.mixest}
\alias{print.mixest}
\alias{print}

\title{Prints \code{mixest} Object.}

\description{The function prints selected outcomes obtained from object \code{mixest}.}

\usage{\method{print}{mixest}(x, ...)}

\arguments{
\item{x}{an object of \code{mixest} class}
\item{...}{not used}
}

\details{
The function prints the general structure of the model, i.e., names of predictors. It also prints the number of observations (length of time-series) and the number of component models used in estimations (mixing). Additionally it prints the model's parameters (i.e., forecasting method, values of the initial parameters, etc.).
}

\examples{
\donttest{
data(oil) 
m1 <- mixest1(y=oil[,1,drop=FALSE],x=oil[,-1,drop=FALSE],ftype=2,V=100,W=100)
print(m1)
}
}
