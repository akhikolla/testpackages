\name{calculateassociations}
\alias{calculateassociations}
\title{Generates social network based on xy spatial coordinates of individuals.}
\usage{
calculateassociations(x, y, ir)
}
\arguments{
  \item{x}{a list of x coordinates for all the individuals}

  \item{y}{a list of y coordinates for all the individuals}

  \item{ir}{a list of interaction radii for all the
  individuals}
}
\description{
This method uses the pairwise distances between each pair
of individuals. For a spatial point pattern X, the
association of individual j on individual i, Aij, is
calculated using the distances between points X_{i} and
X_{j} using a smooth interaction function first introduced
by Illian (2009). If d represents the distance between
points X_{i} and X_{j}, and the interaction radius for
individual i is R_{i}, then the association of j on i is
calculated as: ((1-(d/R_{i})^{2}))^{2} if d>0 and d<=R_{i},
and 0 otherwise. This function has been described as a
smooth interaction function because the value of the
association calculated decreases smoothly as a function of
distance.  This is in contrast to
associations calcuated using a binary function where the
association =1 if d<=R and 0 otherwise.  Such a function is
based on the assumption that the association is constant
(1) at all distances less than R_{i}, and 0 for distances
greater than R_{i}.
The figure below shows a plot of both functions (the binary interaction function is represented by the dotted line) for an example where R=1.

\figure{smooth.jpeg}
}
\references{
Illian, Moller, Waagepetersen, 2009. Hierarchical spatial
point process analysis for a plant community of high
biodiversity.Environ. Ecol. Stat. vol 16, pp 389-405
}

