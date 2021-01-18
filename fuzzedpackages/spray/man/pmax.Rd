\name{pmax}
\alias{pmax}
\alias{pmin}
\alias{maxpair_spray}
\alias{minpair_spray}
\alias{pmax.spray}
\alias{pmin.spray}
\alias{pmax.default}
\alias{pmin.default}
\title{
Parallel maxima and minima for sprays
}
\description{
Parallel (pairwise) maxima and minima for sprays. 
}
\usage{
maxpair_spray(S1,S2)
minpair_spray(S1,S2)
\method{pmax}{spray}(x, ...)
\method{pmin}{spray}(x, ...)
}
\arguments{
  \item{x,S1,S2}{Spray objects}
  \item{\dots}{spray objects to be compared}
}
\details{

Function \code{maxpair_spray()} finds the pairwise maximum for two
sprays.  Specifically, if \code{S3 <- maxpair_spray(S1,S2)}, then
\code{S3[v] == max(S1[v],S2[v])} for every index vector \code{v}.

Function \code{pmax.spray()} is the method for the generic
\code{pmax()}, which takes any number of arguments.  If \code{S3 <-
maxpair_spray(S1,S2,...)}, then \code{S3[v] == max(S1[v],S2[v],...)} for
every index vector \code{v}.

Function \code{pmax.spray()} operates right-associatively:

\code{pmax(S1,S2,S3,S4) == f(S1,f(S2,f(S3,S4)))} where \code{f()} is
short for \code{maxpair_spray()}.  So if performance is important, put
the smallest spray (in terms of number of nonzero entries) last.

In these functions, a scalar is interpreted as a sort of global maximum.
Thus if \code{S3 <- pmax(S,x)} we have \code{S3[v] == max(S[v],x)} for
every index \code{v}.  Observe that this operation is not defined if
\code{x>0}, for then there would be an infinity of \code{v} for which
\code{S3[v] != 0}, an impossibility (or at least counter to the
principles of a sparse array).  Note also that \code{x} cannot have
length \code{>1} as the elements of a spray object are stored in an
arbitrary order.

Functions \code{minpair_spray()} and \code{pmin.spray()} are analogous.
Note that \code{minpair_spray(S1,S2)} is algebraically equivalent to
\code{-pmax_spray(-S1,-S2)}; see the examples.

The value of \code{pmax(S)} is problematic.  Suppose
\code{all(value(S)<0)}; the current implementation returns
\code{pmax(S)==S} but there is a case for returning the null polynomial.
  
}
\value{
  Returns a spray object
}

\author{
Robin K. S. Hankin
}
\examples{

S1 <- rspray(100,vals=sample(100)-50)
S2 <- rspray(100,vals=sample(100)-50)
S3 <- rspray(100,vals=sample(100)-50)


# following comparisons should all be TRUE:

jj <- pmax(S1,S2,S3)
jj ==  maxpair_spray(S1,maxpair_spray(S2,S3))
jj ==  maxpair_spray(maxpair_spray(S1,S2),S3)

pmax(S1,S2,S3)  == -pmin(-S1,-S2,-S3)
pmin(S1,S2,S3)  == -pmax(-S1,-S2,-S3)

pmax(S1,-Inf) == S1
pmin(S1, Inf) == S2


pmax(S1,-3)

\dontrun{
pmax(S1,3) # not defined
}


}
