\name{chess}
\alias{chess}
\docType{data}
\title{Chess playing dataset}
\description{
A tally of wins and losses for games between three chess players:
Topalov, Anand, Karpov
}
\usage{data(chess)}
\details{

  This is a very simple dataset that can be used for illustration of
  \code{hyper2} idiom.
  
  The players are:
  \itemize{
    \item Grandmaster Veselin Topalov. FIDE world champion 2005-2006;
    peak rating 2813
    \item Grandmaster Viswanathan Anand.  FIDE world champion 2000-2002,
    2008; peak rating 2799
    \item Grandmaster Anatoly Karpov.  FIDE world champion 1993-1999;
    peak rating 2780
  }

  Observe that Topalov beats Anand, Anand beats Karpov, and Karpov beats
  Topalov (where \dQuote{beats} means \dQuote{wins more games than}).

  The games thus resemble a noisy version of \dQuote{rock paper scissors}.

  The likelihood function does not record who played white; see
  \code{karpov_kasparov_anand} for such a dataset.

}
\references{
  \itemize{
    \item \url{http://www.chessbase.com/}
  }
}
\seealso{\code{\link{karpov_kasparov_anand}}}
\examples{
data(chess)
maxp(chess)
}
\keyword{datasets}
