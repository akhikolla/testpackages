\name{summary.splmm}
\alias{summary.splmm}
\title{Summarize an 'splmm' object}
\description{Providing an elaborate summary of a \code{'splmm'} object.}
\usage{
\method{summary}{splmm}(object, ...)
}

\arguments{
  \item{object}{a \code{'splmm'} object}
  \item{\dots}{not used.}
}
\details{This functions shows a detailed summary of a \code{'splmm'} object.}

\examples{
\donttest{
data(cognitive)

x <- model.matrix(ravens ~ schoolid+treatment+relyear+sex+age_at_time0+height+weight
+head_circ+ses1+readtest+writetest+yrsofsch, cognitive)
z <- model.matrix(ravens ~ schoolid+treatment+relyear, cognitive)

fit <- splmm(x=x,y=cognitive$ravens,z=z,grp=cognitive$id,lam1=0.1,lam2=0.1,
penalty.b="scad", penalty.L="scad")
summary(fit)
}
}

\keyword{methods}
