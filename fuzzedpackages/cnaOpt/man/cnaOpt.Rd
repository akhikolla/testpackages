\name{cnaOpt}
\alias{cnaOpt}

\title{Find atomic solution formulas with optimal consistency and coverage
}

\description{
\code{cnaOpt} attempts to find atomic solution formulas (asfs) for a given \code{outcome} (inferred from crisp-set, \code{"cs"}, or multi-value, \code{"mv"}, data) that are optimal with respect to consistency and coverage. 
}

\usage{
cnaOpt(x, outcome, ..., crit = quote(con * cov), cond = quote(TRUE))
}

\arguments{
  \item{x}{
A \code{\link{data.frame}} or \code{\link{configTable}} of type \code{"cs"} or \code{"mv"}; if \code{x} is of type \code{"mv"}, \code{x} must be a \code{configTable}.
}
  \item{outcome}{
A character string specifying one outcome, i.e. one factor in \code{x}.
}
  \item{\dots}{
Additional arguments passed to \code{\link{configTable}}, for instance \code{type}, \code{rm.dup.factors}, \code{rm.dup.factors}, or \code{case.cutoff}.
}
  \item{crit, cond}{
A quoted expression specifying additional criteria for selecting optimal solutions, passed to \code{\link{selectMax}}.
}
}

\details{
\code{cnaOpt} infers causal models (atomic solution formulas, asf) for the \code{outcome} from data \code{x} that are optimal with respect to the optimality criterion \code{crit} and complying with conditions \code{cond}. Data \code{x} may be crisp-set (\code{"cs"}) or multi-value (\code{"mv"}), but not fuzzy-set (\code{"fs"}). The function first calculates consistency and coverage optima (con-cov optima) for \code{x}, then selects the optimum that is best according to \code{crit} and \code{cond}, builds the canonical disjunctive normal form (DNF) realizing the best optimum and, finally, generates all minimal forms of that canonical DNF. 

Roughly speaking, running \code{cnaOpt} amounts to sequentially executing \code{\link{configTable}}, \code{\link{conCovOpt}}, \code{\link{selectMax}}, \code{\link{DNFbuild}} and \code{\link{condTbl}}.

}

\value{
\code{cnaOpt} returns a \code{\link{data.frame}} with additional class "condTbl". See the "Value" section in \code{?\link{condTbl}} for details.
}

\seealso{
\code{\link{cna}}, \code{\link{conCovOpt}}
}

\examples{
# Example 1: Real-life crisp-set data, d.educate.
(res_opt1 <- cnaOpt(d.educate, "E"))

# Using the pipe operator (\%>\%), the steps processed by cnaOpt in the 
# call above can be reproduced as follows:
library(dplyr)
conCovOpt(d.educate, "E") \%>\% selectMax \%>\% DNFbuild("E", reduce = "ereduce") \%>\% 
  paste("<-> E") \%>\% condTbl(d.educate)


# Example 2: Simulated crisp-set data.
dat1 <- data.frame(
  A = c(1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0), 
  B = c(0, 0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0), 
  C = c(0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0), 
  D = c(1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1), 
  E = c(1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1), 
  F = c(0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1)
)

(res_opt2 <- cnaOpt(dat1, "E"))

# Change the optimality criterion.
cnaOpt(dat1, "E", crit = quote(pmin(con, cov)))
# Impose an additional condition.
\donttest{cnaOpt(dat1, "E", cond = quote(con >= 0.9))}


# Example 3: All logically possible configurations.
(res_opt3 <- cnaOpt(full.ct(4), "D"))  # All combinations are equally bad.


# Example 4: Real-life multi-value data, d.pban.
cnaOpt(mvct(d.pban), outcome = "PB=1")
\donttest{cnaOpt(mvct(d.pban), outcome = "PB=1", crit = quote(pmin(con, cov)))}
\donttest{cnaOpt(mvct(d.pban), outcome = "PB=1", cond = quote(con > 0.93))}
cnaOpt(mvct(d.pban), outcome = "PB=0")
cnaOpt(mvct(d.pban), outcome = "PB=0", cond = quote(con > 0.93))
cnaOpt(d.pban, type = "mv", outcome = "F=2")
cnaOpt(d.pban, type = "mv", outcome = "F=2", cond = quote(con > 0.75))

}
