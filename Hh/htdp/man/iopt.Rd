\name{iopt}
\alias{iopt}
\title{List reference frames recognized by HTDP}
\description{
Returns an \R \code{data.frame} containing a list of reference frames recognized by HTDP. See Table 4 in the HTDP User's Guide for more information.
}
\usage{
iopt()
}
\value{
  A \code{data.frame} with 2 columns:
    \item{key}{The number that the HTDP software uses to identify a particular reference frame. Some reference frames share the same key.}
    \item{value}{The description of the reference frame.}
}
\seealso{
  \code{\link{displace}}
}
