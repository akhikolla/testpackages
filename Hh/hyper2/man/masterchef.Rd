\name{masterchef}
\Rdversion{1.1}
\alias{masterchef}
\alias{MasterChef}
\alias{masterchef_series6}
\alias{pmax_masterchef6}
\alias{pmax_masterchef6_constrained}
\docType{data}
\title{
Masterchef series 6
}
\description{
Data from Australian Masterchef Series 6
}
\usage{data(masterchef)}
\format{
  Object \code{masterchef_series6} is a list of \code{hyper2} objects;
\code{pmax_masterchef_series6} and
\code{pmax_masterchef_series6_constrained} are named vectors with unit
sum.
  
}
\details{

The object is created using the code in
\code{inst/masterchef_series6_allrounds.R}, which is heavily documented.
Not all the information available is included in the likelihood function
as some of the early rounds result in an unmanageably large list.
Inclusion is controlled by Boolean vector \code{doo}.
 
The definitive source is the coloured table on the wiki page.
}
\references{

  Wikipedia contributors, \dQuote{MasterChef Australia (series 6),}
  Wikipedia, The Free Encyclopedia,
  \url{https://en.wikipedia.org/w/index.php?title=MasterChef_Australia_(series_6)&oldid=758432561}
  (accessed January 5, 2017).
  
}
\seealso{\code{\link{ggol}}}
\examples{
data(masterchef)

a1 <- rep(1/13,12)                          # equal strengths
a2 <- indep(pmax_masterchef6)               # MLE
a3 <- indep(pmax_masterchef6_constrained)   # constrained MLE

like_series(a1, masterchef_series6)
like_series(a2, masterchef_series6)
like_series(a3, masterchef_series6)

}
\keyword{datasets}
