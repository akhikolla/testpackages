\name{change_pnames}
\alias{change_pnames}
\alias{all_pnames}
\title{Change pnames of a hyper2 object}
\description{Change the pnames attribute of a hyper2 object}
\usage{
change_pnames(H, new_pnames)
all_pnames(L)
}
\arguments{
  \item{H}{Object of class hyper2}
  \item{L}{List of hyper2 objects}
  \item{new_pnames}{Character vector of the new pnames}
}
\author{Robin K. S. Hankin}
\note{

  Function \code{change_pnames()} is not a straightforward function.
  The intent is to enable one to add hyper2 objects with non-identical
  pnames attributes in a sensible way.

  The standard use for \code{change_names()} is to add new players to a
  hyper2 object: so the new names must be a superset of the old,
  although the order may differ.

  Function \code{all_pnames()} returns the union of the pnames of the
  hyper2 objects sent to it
  
}
\examples{

x1 <- hyper2(list(1,1:2,1:3),1:3,pnames=letters[1:3])
x2 <- x1
pnames(x2) <- letters[3:1]

x1
x2
# At this point, note that x1 and x2 are algebraically identical; but:
x1==x2

## Variables x1 and x2 are different!  This is by design. The pnames
## attribute is just syntatic sugar, mathematically meaningless.

## Also note that you cannot add x1 and x2 because to do so is bad weirdness:
\dontrun{
x1 + x2
}

## Here is a sensible use-case:

x1 <- hyper2(list(1,1:2,1:3),1:3,pnames=letters[1:3])
x2 <- hyper2(list(1,2),10:11,pnames=letters[4:5])

pnames(x1) <- letters[1:5]
pnames(x2) <- letters[1:5]
x1+x2


## NB:

pnames(x1) <- letters[5:1]
x1


}