\name{rrank}
\alias{rrank}
\alias{print.rrank}
\title{Random ranks}
\description{
A function for producing ranks randomly, consistent with a specified
strength vector
}
\usage{
rrank(n = 1, p, pnames=NULL, fill = FALSE)
\method{print}{rrank}(x, ...)
}
\arguments{
  \item{n}{Number of observations}
  \item{p}{Strength vector}
  \item{pnames}{Character vector specifying names of the columns}
  \item{fill}{Boolean, with default \code{FALSE} meaning to interpret
    the elements of \code{p} as strengths, notionally summing to one;
    and \code{TRUE} meaning to augment \code{p} with a fillup value}
  \item{x, ...}{Arguments passed to the print method}
}

\value{

If \code{n=1}, return a vector; if \code{n>1} return a matrix with
\code{n} rows, each corresponding to a ranking.  The canonical example
is a race in which the probability of competitor \eqn{i} coming first is
\eqn{p_i/\sum p_j}{p_i/sum p}, where the summation is over the
competitors who have not already finished.

If, say, the first row of \code{rrank()} is \code{c(2,5,1,3,4)}, then
competitor 2 came first, competitor 5 came second, and so on.

Note that function \code{rrank()} returns an object of class
\code{rrank}, which has its own special print method.  The column names
appear as \dQuote{\code{c1, c2, ...}} which is intended to be read
\dQuote{came first}, \dQuote{came second}, and so on.  The difference
between \emph{rank} and \emph{order} can be confusing.

\preformatted{
R> x <- c( 3, 1, 2, 4.1)
R> rank(x)
[1] 3 1 2 4
R> order(x)
[1] 2 3 1 4
R> 
}

In the current context, we have:

\preformatted{
R> rrank(2,ptrue)
     c1 c2 c3 c4
[1,]  3  1  2  4
[2,]  2  3  4  1
R> 
}

In the above, each row is a race; we have four runners.  In the first
race, runner number 3 came first, runner 1 came second, runner 2 came
third, and so on.  In the second race, runner 2 came first, etc.  Now
consider the following R idiom which uses the first race:

\preformatted{
R> p <- c(3,1,2,4)
R> order(p)    # slicker to say p[p] <- 1:4 as in inverse_word_single()
[1] 2 3 1 4
R> 
}

This shows that runner 1 came second, runner 2 came third, runner 3 came
first, and runner 4 came fourth.

Take the first race as an example.

Rank: who came first? runner 3.  Who came second? runner 1.  Who came
third? runner 2.  Who came fourth? runner 4.

Order: where did runner 1 come?  second.  Where did runner 2 come?
third.  Where did runner 3 come?  first.  Where did runner 4 come?
fourth.  For \code{order_likelihood()}, we need rank data, not order
data.  See vignette \dQuote{\code{skating_analysis}} for more
discussion.

}
\author{Robin K. S. Hankin}
\seealso{\code{\link{ggol}},\code{\link{order_likelihood}},\code{\link{skating}}}
\examples{
ptrue <- (4:1)/10
rrank(100,p=ptrue)

r1 <- rrank(100,p=ptrue)
r2 <- rrank(10,p=ptrue)

H <- order_likelihood(r1)
maxp(H)   # should be close to ptrue

H <- H + order_likelihood(rrank(100,maxp(H)))
maxp(H)



}