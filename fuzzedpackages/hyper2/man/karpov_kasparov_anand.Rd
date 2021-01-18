\name{karpov_kasparov_anand}
\alias{karpov_kasparov_anand}
\alias{kka}
\alias{kka_3draws}
\alias{kka_3whites}
\alias{white_wins}
\alias{drawn_games}
\alias{black_wins}
\alias{plays_white_wins}
\alias{plays_white_draws}
\alias{plays_white_loses}
\alias{kka_array}
\docType{data}
\title{Karpov, Kasparov, Anand}
\description{
  Data of three chess players: Karpov, Kasparov, and Anand.  Includes
two likelihood functions for the strengths of the players, and matrices
of game results}

\details{
  The strengths of chess players may be assessed using the generalized
  Bradley-Terry model.  The \code{karpov_kasparov_anand} likelihood
  function allows one to estimate the players' strengths, propensity to
  draw, and also the additional strength conferred by playing white.

  Likelihood functions \code{karpov_kasparov_karpov}, \code{kka_3draws}
  and \code{kka_3whites} are created by files
  \code{inst/karpov_kasparov_anand.R}, \code{inst/kka_3draws} and
  \code{inst/kka_3whites}, which are heavily documented and include some
  analysis.  Object \code{karpov_kasparov_anand} assumes that the draw
  potential is the same for all three players; likelihood function
  \code{kka_3draws} allows the propensity to draw to differ between the
  three players.

  The reason that the players are different from those in the
  \code{chess} dataset is that the original data does not seem to be
  available any more.

  Dataset \code{kka} refers to scorelines of matches between three chess
  players (Kasparov, Karpov, Anand).  It is a list with names such as
  \sQuote{\code{karpov_plays_white_beats_kasparov}} which has value 18:
  we have a total of 18 games between Karpov and Kasparov in which
  Karpov played white and beat Kasparov.

  The three matrices \code{plays_white_wins}, \code{plays_white_draws},
  and \code{plays_white_loses} tabulate this information in a coherent
  way; and array \code{kka_array} presents the same information in a 3D
  array (but the names of the dimnames are lost).

  All data drawn from chessgames.com, specifically

 \url{http://www.chessgames.com/perl/ezsearch.pl?search=karpov+vs+kasparov}

 Note that the database allows one to sort by white wins or black wins
 (there is a \sQuote{refine search} tab at the bottom).  Some searches
 have more than one page of results.
 
 Numbers here downloaded 17 February 2019.  Note that only
 \sQuote{classical games} are considered here (rapid and exhibition
 games being ignored).

}
\seealso{\code{\link{chess}}}
\examples{
data(karpov_kasparov_anand)
maxp(karpov_kasparov_anand)
pie(maxp(karpov_kasparov_anand))
}
\keyword{datasets}
