\name{selectMax}
\alias{selectMax}
\alias{multipleMax}

\title{
Select the con-cov optimum from a '\code{conCovOpt}' object that is best according to a specified optimality criterion
}

\description{
\code{selectMax} selects a con-cov optimum from a '\code{conCovOpt}' object that is best according to some specified optimality criterion. \code{multipleMax} checks a '\code{selectMax}' object for multiple solutions with identical values in the optimality criterion. 
}

\usage{
selectMax(x, crit = quote(con * cov), cond = quote(TRUE))
multipleMax(x, outcome)
}

\arguments{
  \item{x}{
In \code{selectMax}: an object output by \code{\link{conCovOpt}}.
In \code{multipleMax}: an object output by \code{selectMax}.
}
  \item{crit}{
A quoted expression specifying the optimality criterion (see examples). 
}
  \item{cond}{
A quoted expression specifying additional constraints imposed on the optimality criterion (see examples).
}
  \item{outcome}{
A character string specifying a single outcome in the original data.
}
}

\details{
While \code{\link{conCovOpt}} identifies \emph{all} con-cov optima in an analyzed data set, \code{selectMax} selects \emph{one} con-cov optimum from a '\code{conCovOpt}'  object that is best according to the optimality criterion specified in the argument \code{crit}. The default is to select a con-cov maximum: An ordered pair (con, cov) of consistency and coverage scores is a \strong{con-cov maximum} for outcome Y=k in data \eqn{\delta} iff (con, cov) is a con-cov optimum for Y=k in \eqn{\delta} with highest product of consistency and coverage (con-cov product). However, the argument \code{crit} allows for specifying any other optimality criterion, e.g. \code{pmin(con, cov), pmax(con, cov)}, etc. If the '\code{conCovOpt}'  object contains multiple outcomes, the selection of a best con-cov optimum is done separately for each outcome.

Whereas \code{selectMax} selects only one con-cov optimum satisfying \code{crit}, \code{multipleMax} selects all elements in an \code{allConCov} list contained in the '\code{conCovOpt}' object reaching identical scores on the optimality criterion. It is executed for \emph{one} \code{outcome} only (see the examples below).

Via the column \code{id} in the output of \code{selectMax} it is possible to select one among many equally good maxima, for instance, by means of \code{\link{reprodAssign}} (see the examples below).
}

\value{
\code{selectMax} returns an object of class 'selectMax'. 

\code{multipleMax} returns a \code{data.frame}.
}

\seealso{
\code{\link{conCovOpt}}, \code{\link{reprodAssign}}

See also examples in \code{\link{conCovOpt}}.


}


\examples{
dat1 <- d.autonomy[15:30, c("EM","SP","CO","AU")]
(cco1 <- conCovOpt(dat1, type = "fs", outcome = "AU"))
selectMax(cco1)
selectMax(cco1, cond = quote(con > 0.95))
selectMax(cco1, cond = quote(cov > 0.98))
selectMax(cco1, crit = quote(pmin(con, cov)))
selectMax(cco1, crit = quote(pmax(con, cov)), cond = quote(cov > 0.9))


# Multiple equally good maxima.
(cco2 <- conCovOpt(dat1, type = "fs", outcome = "AU", allConCov = TRUE)) 
(sm2 <- selectMax(cco2, cond = quote(con > 0.93)))
multipleMax(sm2, "AU")
# Each maximum corresponds to a different rep-assignment, which can be selected
# using the id argument.
reprodAssign(sm2, "AU", id = 10)
reprodAssign(sm2, "AU", id = 11)
reprodAssign(sm2, "AU", id = 13)


}