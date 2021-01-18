\name{cplusplus}
\alias{cplusplus}
\alias{overwrite}
\alias{accessor}
\alias{assigner}
\alias{addL}
\alias{equal}
\alias{equality}
\alias{identityL}
\alias{evaluate}
\alias{differentiate}
\alias{hyper2_accessor}
\alias{hyper2_equal}
\alias{hyper2_addL}
\alias{hyper2_assigner}
\alias{hyper2_differentiate}
\alias{hyper2_evaluate}
\alias{hyper2_identityL}
\alias{hyper2_overwrite}
\title{Wrappers to c calls}
\description{
Various low-level wrappers to c calls, courtesy of Rcpp
}
\usage{
overwrite(L1, powers1, L2, powers2)
accessor(L,powers,Lwanted)
assigner(L,p,L2,value)
addL(L1,p1,L2,p2)
identityL(L,p)
evaluate(L,powers,probs)
differentiate(L, powers, probs, n)
}

\arguments{
  \item{L,L1,L2,Lwanted}{
    Lists with integer vector elements, used to specify the brackets of
    the hyperdirichlet distribution}
  \item{p,p1,p2,powers,powers1,powers2}{
    A numeric vector specifying the powers to which the brackets are raised
  }
  \item{value}{RHS in assignment, a numeric vector}
  \item{probs}{Vector of probabilities for evaluation of log-likelihood}
  \item{n}{Integer specifying component to differentiate with respect to}
}
\details{
These functions are not really intended for the end-user, as
out-of-scope calls may cause crashes.

}
\value{These functions return a named List}
\author{Robin K. S. Hankin}
