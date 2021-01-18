\name{reality_mining_84}
\alias{reality_mining_84}
\docType{data}

\title{Reality Mining small dataset}

\description{
Human contact data among \code{96} students of the Massachusetts Institute of Technology (MIT), collected by the Reality Mining experiment performed in \code{2004} as part of the Reality Commons project. The dataset is a subset and an adaptation of the proximity dataset extracted and made public by KONECT (url provided below).
}

\format{
The object is a binary adjacency cube of size \code{[96x96x84]}. The time frames correspond to intervals of \code{4} hours. The entry \code{[i,j,t]} is equal to \code{1} if students \code{i} and \code{j} were physically close to each other in the time interval \code{[t-1,t]}. The dataset contains only the first two weeks of interactions. The full dataset is available in this package under the name \code{reality_mining_1392}.
}

\usage{data(reality_mining_84)}

\source{\url{	http://realitycommons.media.mit.edu/realitymining.html}

\url{http://konect.uni-koblenz.de/networks/mit}}

\references{
  Nathan Eagle and Alex (Sandy) Pentland. Reality Mining: Sensing complex social systems. \emph{Personal Ubiquitous Computing, 10(4):255--268}, 2006. 
}

\keyword{datasets}

