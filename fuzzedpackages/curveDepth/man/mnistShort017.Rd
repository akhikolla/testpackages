\name{mnistShort017}
\alias{mnistShort017}
\title{A short version of the MNIST data set}
\usage{
data("mnistShort017")
}
\description{First 100 digits from the MNIST data set belonging to classes
  '0', '1', and '7' representable as curves agter skeletonization.}
\format{A list of three elements `0`, `1`, and `7`, each being an array with
  dimensions \code{(28, 28, 100)}. Each slice of each of these arrays, in third
  dimension, contains a single monochrome digit image represented by a
  \code{28x28} indicator matrix. See
  Lafaye De Micheaux, Mozharovksyi and Vimond (2018) and accompanying codes for
  details on preprocessing.}
\source{\url{http://yann.lecun.com/exdb/mnist/}}
\references{LeCun, Y., Bottou, L.,Bengio, Y., and Haffner, P. (1998).
Gradient-based learning applied to document recognition.
\emph{Proceedings of the IEEE}, 86(11), 2278-2324.

Lafaye De Micheaux, P., Mozharovskyi, P. and Vimond, M. (2018).
Depth for curve data and applications.}
\author{Yann LeCun (Courant Institute, NYU),
Corinna Cortes (Google Labs, New York),
Christopher J.C. Burges (Microsoft Research, Redmond),

preprocessing performed by
Myriam Vimond (CREST, Ensai, University of Bretagne Loire).}
