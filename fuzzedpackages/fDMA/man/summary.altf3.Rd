\name{summary.altf3}
\alias{summary.altf3}

\title{Summarizes Outcomes from \code{altf3} Object.}

\description{The function summarizes selected outcomes obtained from \code{\link{altf3}}.}

\usage{\method{summary}{altf3}(object, ...)}

\arguments{
\item{object}{an object of \code{altf3} class}
\item{...}{not used}
}

\details{
The function produces the outcomes as \code{\link{print.altf3}}.

Additionally, it provides mean values of coefficients and how often p-values (averaged over selected window sizes) for t-test of statistical significance for each independent variable in the model are below 1\%, 5\% and 10\%, respectively. 
}

\examples{
\dontrun{
wti <- crudeoil[-1,1]
drivers <- (lag(crudeoil[,-1],k=1))[-1,]
a <- altf3(y=wti,x=drivers,windows=c(36,100,150))

summary(a)
}
}