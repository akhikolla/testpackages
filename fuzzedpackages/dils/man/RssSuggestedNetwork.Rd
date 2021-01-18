\name{RssSuggestedNetwork}
\alias{RssSuggestedNetwork}
\title{Suggest a network with imputed links}
\usage{
  RssSuggestedNetwork(g, rss, q.impute.above = 0.8)
}
\arguments{
  \item{g}{Object type, then description of \code{arg1}.}

  \item{rss}{Object type, then description of \code{arg2}.}

  \item{q.impute.above}{Object type, then description of
  \code{arg3}.}
}
\value{
  list \item{g.imputed}{\code{igraph} contatining the
  original and the new links} \item{g.new}{\code{igraph}
  containing only the new links} \item{g.original}{original
  graph} \item{q.impute.above}{quantile of RSS scores above
  which links should be imputed}
  \item{frac.filled}{fraction of potential links that were
  actually filled with a new link}
}
\description{
  A longer description of the function.  This can be
  perhaps a paragraph, perhaps more than one.
}
\examples{
g <- graph.atlas(128)
\dontrun{plot(g)}

suggested <- RssSuggestedNetwork(g, q.impute.above=.6)
\dontrun{plot(suggested$g.imputed)}
suggested$frac.filled
}
\author{
  Stephen R. Haptonstahl \email{srh@haptonstahl.org}
}
\references{
  \url{http://www.haptonstahl.org/R}
}
\seealso{
  \code{\link{RelationStrengthSimilarity}}
}

