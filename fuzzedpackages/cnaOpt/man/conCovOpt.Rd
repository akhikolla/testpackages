\name{conCovOpt}
\alias{conCovOpt}
\alias{print.conCovOpt}
\alias{plot.conCovOpt}

\title{
Find consistency and coverage optima for configurational data
}

\description{
\code{conCovOpt} issues pairs of optimal consistency and coverage scores that atomic solution formulas (asf) of an outcome inferred from configurational data can possibly reach.
}

\usage{
conCovOpt(x, outcome = NULL, 
          type = if (inherits(x, "configTable")) attr(x, "type") else "cs", 
          maxCombs = 1e+07, approx = FALSE, allConCov = FALSE)
\method{print}{conCovOpt}(x, ...)
\method{plot}{conCovOpt}(x, con = 1, cov = 1, ...)
}

\arguments{
  \item{x}{
A \code{\link{data.frame}} or \code{\link{configTable}}. In the \code{print}- and \code{plot}-method: The output of \code{conCovOpt}.
}
  \item{outcome}{
A character vector of one or several factors in \code{x}.
}
  \item{type}{
A character string specifying the data type; one of \code{"cs"}, \code{"mv"}, \code{"fs"}. Must be specified if the data type is not \code{"cs"}.
}
  \item{maxCombs}{
Maximal number of combinations that will be tested for optimality.
}
  \item{approx}{
Logical; if TRUE, an exhaustive search is only approximated; if FALSE, an exhaustive search is conducted. 
}
  \item{allConCov}{
Logical; if TRUE, all possible con-cov scores are stored as attribute "allConCov"; if FALSE, only the optimal con-cov scores are stored.
}
  \item{con, cov}{Numeric scalars between 0 and 1 indicating consistency and coverage thresholds marking the area of "good" models in a square drawn in the plot. Points within the square correspond to models reaching these thresholds.
}
  \item{\dots}{
Currently not used.
}
}

\details{
Prior to actual CNA or QCA analyses, \code{conCovOpt} calculates consistency and coverage optima for models (i.e. atomic solution formulas, asf) of an \code{outcome} inferred from data \code{x}.

An ordered pair (con, cov) of consistency and coverage scores is a \strong{con-cov optimum} for outcome Y=k in data \code{x} iff it is not excluded (on principled grounds) for an asf of Y=k inferred from \code{x} to reach (con, cov) but excluded to score better on one element of the pair and at least as well on the other.

\code{conCovOpt} calculates con-cov optima by executing the following steps:\enumerate{\item if \code{x} is a data frame, aggregate \code{x} in a \code{configTable}, \item build exo-groups with constant values in all factors other than the \code{outcome}, \item assign output values to each exo-group that reproduce the behavior of \code{outcome} as closely as possible, \item  calculate con-cov scores for each assignment resulting in step 3, \item eliminate all non-optimal scores.}

In case of \code{"cs"} and \code{"mv"} data, at least one actual model (asf) inferrable from \code{x} and reaching an optimum's consistency and coverage scores is guaranteed to exist for every con-cov optimum. The function \code{\link{DNFbuild}} can be used to build these optimal models. The same does not hold for \code{"fs"} data. In case of \code{"fs"} data, it merely holds that the existence of a model reaching an optimum's consistency and coverage scores cannot be excluded prior to an actual application of \code{\link{cna}}. 

}

\value{
An object of class 'conCovOpt'. The exo-groups resulting from step 2 are stored as attribute \code{"exoGroups"}, the lists of output values resulting from step 3 are stored as attribute \code{"reprodList"} (reproduction list), and all possible con-cov scores are stored as attribute \code{"allConCov"}.
}

%\references{
%Evtl. paper erwähnen?}

\seealso{
\code{\link{configTable}}, \code{\link{selectMax}}, \code{\link{DNFbuild}}
}

\examples{
(cco.irrigate <- conCovOpt(d.irrigate))
conCovOpt(d.irrigate, outcome = c("R","W"))
# Plot method.
plot(cco.irrigate)
plot(cco.irrigate, con = .8, cov = .8)


dat1 <- d.autonomy[15:30, c("EM","SP","CO","AU")]
(cco1 <- conCovOpt(dat1, type = "fs", outcome = "AU"))

print(cco1, digits = 3, row.names = TRUE)
plot(cco1)

# Exo-groups (configurations with constant values in all factors other than the outcome).
attr(cco1$A, "exoGroups")

# Rep-list (list of values optimally reproducing the outcome).
attr(cco1$A, "reprodList")

# allConCov (add all possible con-cov scores, not just optimal ones).
cco1_acc <- conCovOpt(dat1, type = "fs", outcome="AU", allConCov = TRUE) 
attr(cco1_acc$A, "allConCov")
# If the allConCov table has been built, it is passed to the output of selectMax().
sm1 <- selectMax(cco1_acc)
attr(sm1$A, "allConCov")


dat2 <- d.pacts
# Maximal number of combinations exceeds maxCombs.
(cco2 <- conCovOpt(dat2, type = "fs", outcome = "PACT")) # Generates a warning
# Increase maxCombs.
\donttest{(cco2_full <- try(conCovOpt(dat2, type = "fs", outcome = "PACT", 
  maxCombs=1e+08))) # Takes long or fails due to memory shortage}
# Approximate an exhaustive search.
(cco2_approx <- conCovOpt(dat2, type = "fs", outcome = "PACT", approx = TRUE))

}
