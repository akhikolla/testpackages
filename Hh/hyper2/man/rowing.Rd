\name{rowing}
\Rdversion{1.1}
\alias{rowing}
\alias{sculling}
\alias{sculls2016}
\alias{allrowers}
\docType{data}
\title{Rowing dataset, sculling}
\description{Data from Men's single sculls, 2016 Summer Olympics}
\usage{data(rowing)}
\format{
\code{sculls2016} is a  \code{hyper2} object that gives a likelihood
function for the 2016 men's sculls.
}
\details{
Object \code{sculls2016} is created by the code in
\code{inst/rowing_analysis.R}.  This reads file \code{inst/rowing.txt},
each line of which is a heat showing the finishing order.  Character vector
\code{allrowers} holds the names of all the rowers.

File \code{inst/rowing_minimal.txt} has the same data but with dominated
players (that is, any group of players none of whom have beaten any
player not in the group) have been removed.  This is because dominated
players have a ML strength of zero.

}
\references{

  Wikipedia contributors, \dQuote{Rowing at the 2016 Summer
  Olympics---Men's single sculls}, \emph{Wikipedia, The Free
  Encyclopedia},
  \url{https://en.wikipedia.org/w/index.php?title=Rowing_at_the_2016_Summer_Olympics_\%E2\%80\%93_Men\%27s_single_sculls&oldid=753517240}
  (accessed December 7, 2016).
  
}
\seealso{\code{\link{ggol}}}
\examples{
data(rowing)
dotchart(maxp(sculls2016))
}
\keyword{datasets}
