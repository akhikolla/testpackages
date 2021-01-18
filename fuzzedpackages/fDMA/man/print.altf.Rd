\name{print.altf}
\alias{print.altf}

\title{Prints \code{altf} Object.}

\description{The function prints selected outcomes obtained from \code{\link{altf}}.}

\usage{\method{print}{altf}(x, ...)}

\arguments{
\item{x}{an object of \code{altf} class}
\item{...}{not used}
}

\details{
The function prints forecast quality measures from \code{x}. For details see \code{\link[forecast]{accuracy}}.
}


\examples{
\dontrun{
wti <- crudeoil[-1,1]
drivers <- (lag(crudeoil[,-1],k=1))[-1,]
a <- altf(y=wti,x=drivers)

print(a)
}
}