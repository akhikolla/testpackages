\name{print.altf4}
\alias{print.altf4}

\title{Prints \code{altf4} Object.}

\description{The function prints selected outcomes obtained from \code{\link{altf4}}.}

\usage{\method{print}{altf4}(x, ...)}

\arguments{
\item{x}{an object of \code{altf4} class}
\item{...}{not used}
}

\details{
The function prints forecast quality measures from \code{x}. For details see \code{\link[forecast]{accuracy}}.
}


\examples{
\dontrun{
wti <- crudeoil[-1,1]
drivers <- (lag(crudeoil[,-1],k=1))[-1,]
a <- altf4(y=wti,x=drivers,windows=c(36,100,150))

print(a)
}
}