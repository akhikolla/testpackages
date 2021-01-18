\name{Ops.hyper2}
\alias{Ops.hyper2}
\alias{Ops}
\alias{hyper2_add}
\alias{hyper2_equality}
\alias{hyper2_prod}
\alias{hyper2_sum_numeric}
\alias{sum.hyper2}
\title{Arithmetic Ops Group Methods for hyper2 objects}
\description{

  Allows arithmetic operators \dQuote{\code{+}}, \dQuote{\code{*}} and
  comparison operators \dQuote{\code{==}} and \dQuote{\code{!=}}, to be
  used for hyper2 objects.

  Specifically, \code{H1 + H2} implements addition of two log-likelihood
  functions, corresponding to incorporation of new observational data;
  and \code{n*H1} implements \code{H1+H1+...+H1}, corresponding to
  repeated observations of the same data.

  There are no unary operations for this class.

}
\usage{
\method{Ops}{hyper2}(e1, e2 = NULL)
\method{sum}{hyper2}(x,...,na.rm=FALSE)
hyper2_add(e1,e2)
hyper2_equality(e1,e2)
hyper2_sum_numeric(H,r)
}
\arguments{
  \item{e1,e2}{Objects of class \code{hyper2}, here interpreted as
    hyperdirichlet distributions}
  \item{x,...,na.rm}{In the \code{sum()} method, objects to be summed;
    \code{na.rm} is currently ignored}
  \item{H,r}{In function \code{hyper2_sum_numeric()}, object \code{H} is
    a \code{hyper2} object and \code{r} is a length-one real vector (a
    number)}
}
\details{

  Testing for equality is not straightforward for two implementation
  reasons.  Firstly, the object itself is stored internally as a
  \code{stl} \code{map}, which does not store keys in any particular
  order; and secondly, the \code{stl} \code{set} class is used for the
  brackets.  A set does not include information about the order of its
  elements; neither does it admit repeated elements.  See examples.

  Function \code{hyper2_sum_numeric()} is defined so that idiom like

  \code{H[p] <- H[p] + 3}

  works (without this, one has to type \code{H[p] <- powers(H[p]) + 3},
  which sucks).
  
}
\value{
Returns a \code{hyper2} object or a Boolean.
}
\author{
Robin K. S. Hankin
}
\examples{

stopifnot(hyper2(list(1,1:2),1:2)==hyper2(list(1:2,1),2:1))  # stl map class
stopifnot(hyper2(list(1,1:2),1:2)==hyper2(list(1,c(1,1,2)),1:2)) # stl set class


n <- 4
H <- hyper2(d=n)
for(i in seq_len(30)){
  jj <- sample(seq_len(n),sample(n,1),replace=FALSE)
  H[jj] <- H[jj] + 1   
}


## without pnames:
x1 <- hyper2(list(1,1:2,1:3),1:3)
x2 <- hyper2(list(4,5),10:11)
x1+x2

## Same as above, but with pnames:
y1 <- hyper2(list(1,1:2,1:3),1:3,pnames=letters[1:3])
y2 <- hyper2(list(4,5),10:11,pnames=letters[1:5])
\dontrun{y1+y2}  # fails, pnames not identical 

pnames(y1) <- letters[1:5]
y1+y2




}
\keyword{symbolmath}
