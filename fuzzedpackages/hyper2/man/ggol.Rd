\name{ggol}
\alias{ggol}
\alias{choose_losers}
\alias{choose_winners}
\alias{general_grouped_order_likelihood}
\alias{order_likelihood}
\title{Order statistics}
\description{
  Various functions for calculating the likelihood function for order
statistics
}
\usage{
ggol(H, ...)
general_grouped_order_likelihood(H, ...)
choose_losers(H, all_players, losers)
choose_winners(H, all_players, winners)
order_likelihood(M,times=1)
}

\arguments{
  \item{H}{Object of class \code{hyper2}}
  \item{\dots}{Numeric or character vectors specifying groups of
    players with  equal rank, with higher-ranking groups coming earlier
    in the argument list}
  \item{all_players,winners,losers}{Numeric or character vectors
    specifying (all,winning,losing) competitors in randomly chosen
    groups}
  \item{M}{In function \code{order_likelihood()}, a matrix with each
    row corresponding to a race.  The columns
    correspond to the finishing order; thus \code{a=M[i,j]} means that
    competitor \code{a} finished in place \code{j} in race \code{i}}
  \item{times}{Vector specifiying the number of times each row is
  observed}
}
\details{

  These functions are designed to return likelihood functions, in the
  form of lists of \code{hyper2()} objects, for typical order statistics
  such as the results of rowing heats or MasterChef tournaments.

  Function \code{ggol()} is an easily-typed alias for
  \code{general_grouped_order_likelihood()}.

  Function \code{addrank()} is now defunct.
  
}
\author{Robin K. S. Hankin}
\seealso{\code{\link{rrank}}}
\examples{

W <- hyper2(pnames=letters[1:7])
W1 <- ggol(W, 'a', letters[2:6],'g')  # 24-element list
W2 <- ggol(W, 'b', letters[3:7],'a')  # 24-element list

like_single_list(rep(1/7,6),W1)
like_series(rep(1/7,6),list(W1,W2))

# Run 10 races:
r1 <- rrank(10,p=(7:1)/28)
colnames(r1) <- letters[1:7]

# Likelihood function for r1:
W <- order_likelihood(r1)

# Incorporate information from another (independent) observation:
W <- W + order_likelihood(sample(7))


# consistency check:
H <- hyper2(pnames=letters[1:5])
H <- H + order_likelihood(rrank(100,5:1)) # probs = 5/15,4/15,...,1/15
maxp(H)  # should be close to (5:1)/15




}
