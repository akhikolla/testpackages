\name{print.altf3}
\alias{print.altf3}

\title{Prints \code{altf3} Object.}

\description{The function prints selected outcomes obtained from \code{\link{altf3}}.}

\usage{\method{print}{altf3}(x, ...)}

\arguments{
\item{x}{an object of \code{altf3} class}
\item{...}{not used}
}

\details{
The function prints forecast quality measures from \code{x}. For details see \code{\link[forecast]{accuracy}}.
}


\examples{
\dontrun{
wti <- crudeoil[-1,1]
drivers <- (lag(crudeoil[,-1],k=1))[-1,]
a <- altf3(y=wti,x=drivers,windows=c(36,100,150))

print(a)
}
}