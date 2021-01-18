\name{interzonal}
\alias{interzonal}
\alias{curacao}
\alias{stockholm1962}
\alias{collusion}
\docType{data}
\title{1963 World Chess Championships}
\description{
Likelihood functions for players' strengths in the fith Interzonal
tournament which occurred as part of the 1963 Chess world Championships
in Stockholm, 1962.
}
\details{

  The 1963 World Chess Championships was notable for allegations of
  Soviet collusion.  Specifically, Fischer publicly alledged that
  certain Soviet players had agreed in advance to draw all their games.

  Likelihood functions \code{interzonal} and \code{interzonal_collusion}
  are created by files \code{inst/collusion_stockholm.R}, which are
  heavily documented and include some analysis.  Object
  \code{interzonal} includes a term for drawing, (\dQuote{\code{draw}}),
  assumed to be the same for all players; object \code{collusion}
  includes in addition to \code{draw}, a term for the drawing in
  Soviet-Soviet matches, \dQuote{\code{coll}}. 
}
\seealso{\code{\link{chess}},\code{\link{karpov_kasparov_anand}}}
\examples{
data(interzonal)
dotchart(maxp(interzonal))

maxp(collusion)
}
\keyword{datasets}
