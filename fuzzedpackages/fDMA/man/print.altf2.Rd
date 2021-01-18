\name{print.altf2}
\alias{print.altf2}

\title{Prints \code{altf2} Object.}

\description{The function prints selected outcomes obtained from \code{\link{altf2}}.}

\usage{\method{print}{altf2}(x, ...)}

\arguments{
\item{x}{an object of \code{altf2} class}
\item{...}{not used}
}

\details{
The function prints forecast quality measures from \code{x}. For details see \code{\link[forecast]{accuracy}}.
}


\examples{
\dontrun{
wti <- crudeoil[-1,1]
drivers <- (lag(crudeoil[,-1],k=1))[-1,]
a <- altf2(y=wti,x=drivers)

print(a)
}
}