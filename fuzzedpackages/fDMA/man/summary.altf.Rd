\name{summary.altf}
\alias{summary.altf}

\title{Summarizes Outcomes from \code{altf} Object.}

\description{The function summarizes selected outcomes obtained from \code{\link{altf}}.}

\usage{\method{summary}{altf}(object, ...)}

\arguments{
\item{object}{an object of \code{altf} class}
\item{...}{not used}
}

\details{
The function produces the outcomes as \code{\link{print.altf}}.

Additionally, it provides mean values of coefficients and how often p-values for t-test of statistical significance for each independent variable in the model are below 1\%, 5\% and 10\%, respectively. 
}

\examples{
\dontrun{
wti <- crudeoil[-1,1]
drivers <- (lag(crudeoil[,-1],k=1))[-1,]
a <- altf(y=wti,x=drivers)

summary(a)
}
}