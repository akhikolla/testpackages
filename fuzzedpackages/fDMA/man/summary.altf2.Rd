\name{summary.altf2}
\alias{summary.altf2}

\title{Summarizes Outcomes from \code{altf2} Object.}

\description{The function summarizes selected outcomes obtained from \code{\link{altf2}}.}

\usage{\method{summary}{altf2}(object, ...)}

\arguments{
\item{object}{an object of \code{altf2} class}
\item{...}{not used}
}

\details{
The function produces the outcomes as \code{\link{print.altf2}}.

Additionally, it provides mean values of coefficients, min, max and mean relative variable importance for each independent variable in the model, frequency when relative variable importance is over 0.5 for each independent variable in the model, and how often p-values (averaged over selected models) for t-test of statistical significance for each independent variable in the model are below 1\%, 5\% and 10\%, respectively. 
}

\examples{
\dontrun{
wti <- crudeoil[-1,1]
drivers <- (lag(crudeoil[,-1],k=1))[-1,]
a <- altf2(y=wti,x=drivers)

summary(a)
}
}