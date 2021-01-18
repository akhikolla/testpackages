\name{summary.altf4}
\alias{summary.altf4}

\title{Summarizes Outcomes from \code{altf4} Object.}

\description{The function summarizes selected outcomes obtained from \code{\link{altf4}}.}

\usage{\method{summary}{altf4}(object, ...)}

\arguments{
\item{object}{an object of \code{altf4} class}
\item{...}{not used}
}

\details{
The function produces the outcomes as \code{\link{print.altf4}}.

Additionally, it provides mean values of coefficients.
}

\examples{
\dontrun{
wti <- crudeoil[-1,1]
drivers <- (lag(crudeoil[,-1],k=1))[-1,]
a <- altf4(y=wti,x=drivers,windows=c(36,100,150))

summary(a)
}
}