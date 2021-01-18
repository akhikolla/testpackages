\name{preprocessbboptions}
\alias{preprocessbboptions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Set controls for most functiosn in the package
}
\description{
 Preprocesses a list of argument. The return value of this function serves as argument to \code{\link{blackbox.options}} (see Examples). Providing in this way the information described in the Details section of \code{blackbox.options} is essential for further usage of the package functions.
}
\usage{
preprocessbboptions(optionList)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{optionList}{
A list, with named elements, which names will (mostly) match the names of options set by this function
}
}
\value{
A list, returned invisibly
}
\examples{
\dontrun{
GP <- list(ParameterNames=c("theta_1","theta_2"))
pp <- preprocessbboptions(GP)
do.call(blackbox.options, pp) ## essential
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
