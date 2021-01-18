\name{increment}
\alias{increment}
\alias{decrement}
\alias{trial}
\alias{inc}
\alias{dec}
\title{Increment and decrement operators}
\description{
Syntactic sugar for incrementing and decrementing likelihood functions
}
\usage{
inc(H, val = 1)
dec(H, val = 1)
trial(H,winners,players,val=1)
}
\arguments{
  \item{H}{A hyper2 object}
  \item{winners,players}{Numeric or character vectors specifying the
    winning team and the losing team}
  \item{val}{Numeric}
}
\details{
A very frequent operation is to increment a single term in a hyper2
object.  If

\preformatted{
R> H <- hyper2(list(1,2,3,1:2,2:3),1:5)
R> pnames(H) <- letters[1:3]
}

Then we would have

\preformatted{
R> H
a * (a + b)^4 * b^2 * (b + c)^5 * c^3
}

Suppose we wish to increment the power of \code{a+b}.  We could do:

\preformatted{
H[c("a","b")] <- H[c("a","b")] + 1
}

Or use \code{magrittr} pipes:

\preformatted{
H[c("a","b")] \%<>\% `+`(1)
}


But \code{inc} and \code{dec} furnish convenient idiom to accomplish the
same thing:

\preformatted{
H[c("a","b")] \%<>\% inc
}

Functions \code{inc} and \code{dec} default to adding or subtracting 1,
but other values can  be supplied:

\preformatted{
H[c("a","b")] \%<>\% inc(3)
}


Or even 

\preformatted{
H[c("a","b")] \%<>\% inc(H["a"])
}


The convenience function \code{trial()} takes this one step further and
increments the \sQuote{winning team} and decrements the \sQuote{losing
team}.  This function is intended to be used with magrittr pipes:


\preformatted{
> H <- hyper2(pnames=letters[1:6])
> H
(a + b + c + d + e + f)^0
> H \%<>\% trial(c("a","b"),c("a","b","c"))
> H
(a + b) * (a + b + c)^-1
}

Using \code{trial()} in this way ensures that the powers sum to zero.

The \code{inc} and \code{dec} operators are used in
\code{inst/rowing_analysis.R}; and the \code{trial()} function is used
in \code{inst/kka_draws.R}. 

}
\author{Robin K. S. Hankin}
\examples{
data(chess)

## Now suppose Topalov beats Anand. To  incoroprate this into the LF:
chess["Topalov"] \%<>\% inc
chess[c("Topalov","Anand")] \%<>\% dec
chess
}