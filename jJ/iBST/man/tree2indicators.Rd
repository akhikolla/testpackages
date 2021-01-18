\name{tree2indicators}
\alias{tree2indicators}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
From a tree to indicators (or dummy variables)
}
\description{
Coerces a given tree structure inheriting from rpart to binary covariates.
}
\usage{
tree2indicators(fit)
}
\arguments{
  \item{fit}{
a tree structure inheriting to the rpart method
}
}

\value{
a list of indicators defining the leaf nodes of the fitted tree from left to right
}

\author{
Cyprien Mbogning
}

\examples{
fit <- rpart(Kyphosis ~ Age + Number + Start, data = kyphosis)
tree2indicators(fit)
}

\keyword{ documentation }% use one of  RShowDoc("KEYWORDS")
\keyword{ tree }% __ONLY ONE__ keyword per line
