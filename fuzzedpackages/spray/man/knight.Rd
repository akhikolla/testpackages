\name{knight}
\alias{knight}
\alias{king}
\alias{chess_knight}
\title{
Generating function for a chess knight and king
}
\description{
Generating function for a chess knight and king on an
arbitrarily-dimensioned chessboard
}
\usage{
knight(d=2)
king(d=2)
}
\arguments{
  \item{d}{Dimensionality of the board, defaulting to 2}
}
\value{
Returns the generating function of the piece in question.
}
\author{
Robin K. S. Hankin
}
\note{
The pieces are forced to move; if they have the option of not moving,
add 1 to the returned spray.  The vignette contains a short discussion.
}

\examples{

## How many ways can a knight return to its starting square in 6 moves?
constant(knight()^6)

## How many in 6 or fewer?
constant((1+knight())^6)

## Where does a randomly-moving knight end up?
d <- xyz(2)
kt <- (1+knight())*d^2/9
persp(1:25,1:25,as.array(d*kt^6))


## what is the probability that a 4D king is a knight's move from
##   (0,0,0,0) after 6 moves?

sum(value(((king(4)/80)^4)[knight(4)]))

}
