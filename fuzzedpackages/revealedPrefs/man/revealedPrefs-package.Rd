\name{revealedPrefs-package}
\alias{revealedPrefs-package}
\alias{revealedPrefs}
\docType{package}
\title{
\packageTitle{revealedPrefs}
}
\description{
\packageDescription{revealedPrefs}
}
\details{
This package is meant for the analysis of (quantity, price) data, eg. of bundles of goods and corresponding prices. It features fast algorithms that make the analysis of large datasets feasible.

Functions \code{directPrefs} and \code{indirectPrefs} compute revealed preferences.

Functions \code{checkWarp}, \code{checkSarp}, \code{checkGarp} perform fast non-parametric tests of rationality using the corresponding rationality axioms.

Functions \code{simWarp}, \code{simSarp}, \code{simGarp} and \code{simPrefs} generate simulated data consistent with a rationality axiom or with a given preference matrix.

Functions \code{cpLower} and \code{cpUpper} generate Crawford-Pendakur type bounds on the number of subpopulations and provide the corresponding clusterings.
}
\author{
\packageAuthor{revealedPrefs}

Maintainer: \packageMaintainer{revealedPrefs}
}
\references{
Varian, H. R. (1982) The Nonparametric Approach to Demand Analysis, 
\emph{Econometrica}, 50(4):945-973.

Varian, H. R. (1984) \emph{Microeconomic Analysis}. New York/London:
Norton, 2nd edition, pp 141-143.

Crawford, I. and Pendakur, K. (2013). How many types are there?
\emph{The Economic Journal}, 123(567):77-95.
}
\keyword{ package }
\seealso{
See \link{directPrefs} for computation of preferences, \link{checkGarp} for rationality tests, \code{simGarp} for data generation, and \link{cpUpper} for clustering of data into non-violating subsets.
}
\examples{
# Compute preferences and check rationality on a GARP-violating dataset:
data(noGarp)
indirectPrefs(noGarp$x, noGarp$p)
checkGarp(noGarp$x, noGarp$p)

# Cluster dataset into GARP-consistent subpopulations:
cpUpper(noGarp$x, noGarp$p)
}
