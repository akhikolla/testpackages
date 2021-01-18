\name{head.hyper2}
\alias{head.hyper2}
\title{First few terms of a  distribution}
\description{First few terms in a hyperdirichlet distribution}
\usage{\method{head}{hyper2}(x, ...)}
\arguments{
  \item{x}{Object of class \code{hyper2}}
  \item{\dots}{Further argments, passed to \code{head()}}
}
\details{Function is  \code{x[head(brackets(x), ...)]}}
\value{Returns a \code{hyper2} object}
\author{Robin K. S. Hankin}
\examples{

H <- order_likelihood(rrank(100,p=(10:1)/55))
head(H)

}
