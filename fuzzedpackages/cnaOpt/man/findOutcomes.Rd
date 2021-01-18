\name{findOutcomes}
\alias{findOutcomes}

\title{
Identify the factors that can possibly be modeled as outcomes prior to running CNA
}
\description{
Prior to running CNA (or any other configurational comparative method), \code{findOutcomes} identifies those factors in data \code{x} that can be modeled as outcomes relative to specified consistency and coverage thresholds \code{con} and \code{cov}. 
%Determine potential outcomes in a \code{configTable} fulfilling given con and cov thresholds.
}

\usage{
findOutcomes(x, con = 1, cov = 1, ...)
}

\arguments{
  \item{x}{
A \code{\link{data.frame}} or \code{\link{configTable}}; if \code{x} is of type multi-value (\code{"mv"}) or fuzzy-set (\code{"fs"}), \code{x} must be a \code{configTable}.
}
  \item{con, cov}{
Numeric scalars between 0 and 1 specifying consistency and coverage thresholds.
}
  \item{\dots}{
Additional arguments passed to \code{\link{configTable}}, for instance \code{type}, \code{rm.dup.factors}, \code{rm.dup.factors}, or \code{case.cutoff}.
}
}

\details{
\code{findOutcomes} first runs \code{\link{conCovOpt}} to find the con-cov optima for all factors in \code{x} and then applies \code{\link{selectMax}} to select those factors with con-cov optima meeting the consistency and coverage thresholds specified in \code{con} and \code{cov}.

In case of \code{"cs"} and \code{"mv"} data, an actual model (asf) meeting the specified \code{con} and \code{cov} thresholds is guaranteed to exist for every factor value with an entry \code{TRUE} in the \code{outcome} column. The function \code{\link{DNFbuild}} can be used to build these models. The same does not hold for \code{"fs"} data. In case of \code{"fs"} data, an entry \code{TRUE} in the \code{outcome} column simply means that the existence of a model reaching the specified \code{con} and \code{cov} thresholds cannot be excluded prior to an actual application of \code{\link{cna}}. 
}

\value{
A data.frame.
}

\seealso{
\code{\link{conCovOpt}}, \code{\link{selectMax}}, \code{\link{selectCases}}, \code{\link{DNFbuild}}, \code{\link{full.ct}}
}

\examples{
# CS data
findOutcomes(d.educate)
findOutcomes(d.educate, con = 0.75, cov = 0.75)

# A causal chain.
target1 <- "(A + B <-> C)*(C + D <-> E)"
dat1 <- selectCases(target1, full.ct(target1))
findOutcomes(dat1)

# A causal cycle.
target2 <- "(A + Y1 <-> B)*(B + Y2 <-> A)*(A + Y3 <-> C)" 
dat2 <- selectCases(target2, full.ct(target2))
findOutcomes(dat2)


# MV data
findOutcomes(mvct(d.pban)) # no possible outcomes at con = cov = 1 
findOutcomes(mvct(d.pban), con = 0.8) 
findOutcomes(d.pban, type = "mv", con = 0.8, cov= 0.8) 


# FS data
findOutcomes(fsct(d.jobsecurity)) # no possible outcomes at con = cov = 1 
findOutcomes(d.jobsecurity, type = "fs", con = 0.86) 


}
