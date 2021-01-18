\name{convts}
\alias{convts}

\title{Renames Selected Outcomes of \code{mixest} and \code{tvpreg} Objects.}

\description{
This function renames rows of selected outcomes stored in \code{mixest} and \code{tvpreg} objects. It can be useful in generating better looking plots.
}

\usage{
convts(x,ind=NULL, ...)
}

\arguments{
\item{x}{object of class \code{mixest} or \code{tvpreg}}
\item{ind}{optional, \code{\link[base]{character}} consisting of names of time points, should have the same length as the forecasted time-series}
\item{...}{optional, alternatively, instead of providing \code{ind}, arguments of \code{\link[base]{seq.Date}} can be specified}
}

\value{object of the same class as \code{x} but with renamed rownames of selected outcomes
}

\examples{
\donttest{
data(oil)
t1 <- tvp.reg(y=oil[,1,drop=FALSE],x=oil[,-1,drop=FALSE],lambda=0.99,V=100,W=100)
plot(t1)

t1a <- convts(x=t1,from=as.Date("1990-02-15"),by="month",length.out=nrow(oil[,1,drop=FALSE]))
plot(t1a)

m1 <- mixest1(y=oil[,1,drop=FALSE],x=oil[,-1,drop=FALSE],ftype=1,V=100,W=100)
plot(m1)

ind <- as.character(seq(from=as.Date("1990-02-15"),by="month",length.out=nrow(oil[,1,drop=FALSE])))
m1a <- convts(x=m1,ind=ind)
plot(m1a)
}
}
