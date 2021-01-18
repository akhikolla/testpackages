\name{infimalRank}
\alias{infimalRank}
\title{Adjusted Ranking of Functional Data Based on the Infimal Depth}
\usage{
infimalRank(ID, IA, ties.method = "max")
}
\arguments{
\item{ID}{The vector of infimal depths of the curves of length \code{n}.}

\item{IA}{The vector of the infimal areas corresponding to the infimal depths from \code{ID}
of length \code{n}.}

\item{ties.method}{Parameter for breaking ties in infimal area index. By default \code{max}, see 
\code{rank}.}
}
\value{
A vector of length \code{n}. Low depth values mean high ranks, i.e. potential outlyingness. 
If some of the infimal depths are identical, the ranking of these functions is made according to the 
values of the infimal area. There, higher infimal area index means higher rank, i.e. non-centrality.
}
\description{
Returns a vector of adjusted depth-based ranks for infimal depth for functional data.
}
\details{
Infimal depths for functional data tend to give to many functional observations the same 
value of depth. Using this function, the data whose depth is the same is ranked according
to the infimal area indicator. This indicator is provided in functions \code{depthf.fd1} along
the value of the infimal depth.
}
\examples{
datafA = dataf.population()$dataf[1:20]
datafB = dataf.population()$dataf[21:50]
D = depthf.fd1(datafA,datafB)
infimalRank(D$Half_ID,D$Half_IA) 

ID = c(0,1,0,0,0,1,1)
IA = c(2,3,1,0,2,4,1)
infimalRank(ID,IA)
}
\references{
Nagy, S., Gijbels, I. and Hlubinka, D.  (2017).
Depth-based recognition of shape outlying functions. 
\emph{Journal of Computational and Graphical Statistics}, \bold{26} (4), 883--893.
}
\author{
Stanislav Nagy, \email{nagy at karlin.mff.cuni.cz}
}
\keyword{depth}
\keyword{functional}
\keyword{rank}
