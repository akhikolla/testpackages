\name{reprodAssign}
\alias{reprodAssign}
\alias{DNFbuild}

\title{
Build disjunctive normal forms realizing con-cov optima
}

\description{
\code{reprodAssign} generates the output values of a disjunctive normal form (DNF) reaching a con-cov optimum. \code{DNFbuild} builds a DNF realizing a targeted con-cov optimum; it only works for crisp-set and multi-value data.
}

\usage{
reprodAssign(x, outcome, id = xi$id)
DNFbuild(x, outcome, reduce = c("rreduce", "ereduce", "none"), id = xi$id)
}

\arguments{
  \item{x}{
An object produced by \code{\link{selectMax}}.  
}
  \item{outcome}{
A character string specifying \emph{one} outcome in \code{attr(x, "configTable")}.
}
  \item{id}{
An integer referring to the identifier of the targeted con-cov optimum.
}
  \item{reduce}{
A character string: if \code{"rreduce"} or \code{"ereduce"}, the canonical DNF realizing the con-cov optimum is freed of redundancies using \code{\link{rreduce}} or \code{\link{ereduce}}, respectively; if \code{"none"}, the unreduced canonical DNF is returned. \code{reduce=TRUE} is interpreted as \code{"rreduce"}, \code{reduce=FALSE} and \code{reduce=NULL} as \code{"none"}.
}
}

\details{An atomic CNA model (asf) accounts for the behavior of the \code{outcome} in terms of a redundancy-free DNF. \code{reprodAssign} generates the output values such a DNF has to return in order to reach a con-cov optimum stored in an object of class '\code{selectMax}'. If the data stored in \code{attr(x, "configTable")} is of type \code{"cs"} or \code{"mv"}, \code{DNFbuild} builds a concrete DNF realizing the targeted con-cov optimum. (For data of type \code{"fs"} an error is returned.) If \code{reduce = "rreduce"} (default), \emph{one} (randomly selected) redundancy-free DNF is built using \code{\link{rreduce}}; if \code{reduce = "ereduce"}, \emph{all} redundancy-free DNFs are built using \code{\link{ereduce}}; if \code{reduce = "none"}, the non-reduced canonical DNF is returned.
The argument \code{id} allows for selecting a targeted con-cov optimum via its identifier (see examples below).
}

\value{
\code{reprodAssign}: A vector of scores.
\code{DNFbuild}: A Boolean formula in disjunctive normal form (DNF).
}


\seealso{
\code{\link{conCovOpt}}, \code{\link{selectMax}}, \code{\link{condTbl}}, \code{\link{reprodAssign}}
}

\examples{
# CS data, d.educate
cco1 <- conCovOpt(d.educate)
best1 <- selectMax(cco1)
reprodAssign(best1, outcome =  "E")
DNFbuild(best1, outcome =  "E")
DNFbuild(best1, outcome =  "E", reduce = FALSE) # canonical DNF
DNFbuild(best1, outcome =  "E", reduce = "ereduce") # all redundancy-free DNFs
DNFbuild(best1, outcome =  "E", reduce = "rreduce") # one redundancy-free DNF
DNFbuild(best1, outcome =  "E", reduce = "none") # canonical DNF


# Simulated mv data
datMV <- data.frame(
  A = c(3,2,1,1,2,3,2,2,2,1,1,2,3,2,2,2,1,2,3,3,3,1,1,1,3,1,2,1,2,3,3,2,2,2,1,2,2,3,2,1,2,1,3,3),
  B = c(1,2,3,2,1,1,2,1,2,2,3,1,1,1,2,3,1,3,3,3,1,1,3,2,2,1,1,3,3,2,3,1,2,1,2,2,1,1,2,2,3,3,3,3),
  C = c(1,3,3,3,1,1,1,2,2,3,3,1,1,2,2,2,3,1,1,2,1,2,2,3,3,1,2,2,2,3,2,1,1,2,2,2,1,1,1,2,2,1,1,2),
  D = c(3,1,2,2,1,1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,1,1,1,1,1,2,2,2,2,2,3,1,1,1,1,1,2,2,2,2,2,3,3,3),
  E = c(3,2,2,3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3)
)

# Apply conCovOpt and selectMax.
cco2 <- conCovOpt(mvct(datMV))
best2 <- selectMax(cco2)

# Apply DNFbuild to build one redundancy-free DNF reaching best2.
(formula1 <- DNFbuild(best2, outcome = "D=3"))
# formula1 reaches the con-cov score stored in best2 for outcome "D=3".
condTbl(paste0(formula1, "<-> D=3"), mvct(datMV))
# Build all redundancy-free DNFs reaching best2.
DNFbuild(best2, outcome = "D=3", reduce = "ereduce")
# Any factor value in datMV can be treated as outcome.
(formula2 <- DNFbuild(best2, outcome =  "E=3"))
condTbl(paste0(formula2, "<-> E=3"), mvct(datMV))
# Any con-cov optimum in cco2 can be targeted via its identifier.
(formula3 <- DNFbuild(best2, outcome =  "E=3", id = 508))
condTbl(paste0(formula3, "<-> E=3"), mvct(datMV))


# Simulated fs data
datFS <- data.frame(
  A = c(.73, .85, .94, .36, .73, .79, .39, .82, .15, .12, .67, .27, .3), 
  B = c(.21, .03, .91, .64, .39, .12, .06, .7,  .73, .15, .88, .73, .36), 
  C = c(.61,  0,  .61,  1,  .94, .15, .88, .27, .12, .12, .27, .15, .15), 
  D = c(.64, .67, .3,  .06, .33, .03, .76, .94, .67, .76, .18, .27, .36), 
  E = c(.91, .94, .67, .85, .73, .79, .24, .09, .03, .21, .33, .36, .27)
)

# Apply conCovOpt and selectMax.
cco3 <- conCovOpt(datFS, outcome = "E", type = "fs", allConCov = TRUE)
best3 <- selectMax(cco3)

# Apply reprodAssign.
reprodAssign(best3, outcome = "E")
# Select a con-cov optimum in cco3 via its identifier.
reprodAssign(best3, outcome = "E", id = 252)

# DNFbuild does not work for fs data; it generates an error.
try(DNFbuild(best3, outcome = "E"))

}
