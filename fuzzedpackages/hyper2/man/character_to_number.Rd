\name{character_to_number}
\alias{character_to_number}
\alias{char2num}
\title{
Convert a character vector to a numeric vector
}
\description{
Convert string descriptions of competitors into their number
}
\usage{
character_to_number(char, pnames)
char2num(char, pnames)
}
\arguments{
  \item{char}{Character vector to be converted}
  \item{pnames}{Names vector (usually \code{pnames(H)})}
}
\details{

The internal mechanism of functions such as \code{ggol()}, and all the
\code{C++} code, operates with the competitors labelled with a
non-negative integer; it is then natural to refer to the competitors as
\code{p1}, \code{p2}, etc.

However, sometimes the competitors have names (as in, for example, the
\code{rowing} dataset).  If so, it is more natural to refer to the
competitors using their names rather than an arbitrary integer.

Function \code{character_to_number()} converts the names to numbers.  If
an element of \code{char} is not present in \code{pnames}, an error is
returned.

Function \code{char2num()} is an easy-to-type synonym.

}
\author{Robin K. S. Hankin}
\seealso{\code{\link{order_likelihood}}}
\examples{

H <- hyper2(pnames=letters[1:9])

# Use numbers:
H + order_likelihood(sample(9))  


# Use letters:
H + order_likelihood(character_to_number(sample(letters[1:9]),pnames(H)))



}