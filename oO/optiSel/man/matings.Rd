\name{matings}
\Rdversion{1.1}
\alias{matings}
\title{Mate Allocation}

\description{
Males and females are allocated for mating such that all breeding animals have the desired number of matings. The mean inbreeding coefficient in the offspring is minimized if matrix \code{Kin} contains pairwise kinships of the selection candidates.
}
\usage{
matings(phen, Kin,  alpha=1, 
    ub.n=NA, max=FALSE, solver="default", ...)}

\arguments{
\item{phen}{Data frame with desired number of matings (column \code{n}), sexes (column \code{Sex}), and IDs (column \code{Indiv}) of the selection candidates. The data frame may also contain column \code{herd} containing the names of the herds to which the females belong (\code{NA} for males).}
\item{Kin}{Kinship matrix (or an other similarity matrix) for selection candidates.}
\item{alpha}{If \code{alpha<1} then the proportion of matings with the same male is at most \code{alpha} in each herd. A value \code{alpha<1} increases genetic connectedness between herds and enables to estimate more accurate breeding values.}
\item{ub.n}{Maximum number of matings of the same individuals. Without this constraint (i.e. \code{ub.n=NA}), some superior animals may always be mated to the same inferior animal, so their offspring would likely not be suitable for breeding.}
\item{max}{The default \code{max=FALSE} means that the objective function is minimized.}
\item{solver}{Either \code{solver="default"}, or \code{solver=Rsymphony_solve_LP}. The latter is possible only if package \code{Rsymbhony} is loaded, which is not available for all platforms. }
\item{...}{Further optimization parameters. By default, they are passed to function \link[ECOSolveR]{ecos.control}}.
}

\details{
Males and females are allocated for mating such that all breeding animals have the desired number of matings. If \code{Kin} is a kinship matrix, then the mean inbreeding coefficient in the offspring is minimized. In general, the mean similarity of the parents is minimized.

The maximum number of matings of the same individuals can be constrained. For each herd, the proportion \code{alpha} of matings with the same male can be constrained as well, but this increases computation time. 

}

\value{
Data frame with columns \code{Sire}, \code{Dam}, \code{n}, and possibly \code{herd}, whereby column \code{n} contains the desired number of matings, and column \code{herd} contains the herd of the dam.

The data frame has attributes \code{objval} with the value of the objective function (usually the mean inbreeding coefficient), and attribute \code{info} describing the solution as reported by the solver.

}


\examples{

data("map")
data("Cattle")
dir   <- system.file("extdata", package = "optiSel")
files <- paste(dir, "/Chr", 1:2, ".phased", sep="")

sKin  <- segIBD(files, map, minSNP=20, minL=2.0)
Phen  <- Cattle[Cattle$Breed=="Angler", ]

cont  <- data.frame(
  age   = c(   1,    2,    3,    4,    5,    6), 
  male  = c(0.11, 0.11, 0.10, 0.08, 0.06, 0.04),
  female= c(0.11, 0.11, 0.10, 0.08, 0.06, 0.04))

cand  <- candes(phen=Phen, sKin = sKin, cont=cont)
con   <- list(uniform="female", ub.sKin = 0.047)
Offspring <- opticont("max.BV", cand, con, trace=FALSE)


#####  Minimize inbreeding   #####
Candidate   <- Offspring$parent
Candidate$n <- noffspring(Candidate, N=20)$nOff
Mating      <- matings(Candidate, sKin)
Mating 
attributes(Mating)$objval


\dontrun{
## This is faster but not available on all platforms:
library("Rsymphony")
Mating <- matings(Candidate, sKin, alpha=0.30, solver=Rsymphony_solve_LP)
Mating
attributes(Mating)$objval

attributes(Mating)$info
#[1] "Optimum solution found"
}

}


\author{Robin Wellmann}
