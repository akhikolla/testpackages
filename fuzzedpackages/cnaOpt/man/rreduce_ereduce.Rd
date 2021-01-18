\name{rreduce_ereduce}
\alias{ereduce}
\alias{rreduce}

\title{Eliminate redundancies from disjunctive normal forms (DNF)
}

\description{\code{ereduce} and \code{rreduce} implement different algorithmic approaches to eliminate redundancies from disjunctive normal forms (DNF), i.e. disjunctions of conjunctions of literals. If there are several minimal DNF, \code{ereduce} will return them all, while \code{rreduce} selects one at random.
}

\usage{
ereduce(cond, x = full.ct(cond), full = !missing(x), simplify2constant = TRUE)
rreduce(cond, x = full.ct(cond), full = !missing(x), verbose = FALSE, maxiter = 1000,
        simplify2constant = TRUE) 
}

\arguments{
  \item{cond}{
A character string specifying a disjunctive normal form; can be either crisp-set or multi-value.
}
  \item{x}{A 
\code{truthTab} or \code{data.frame}; can be either crisp-set or multi-value.
}
  \item{full}{
Logical; if \code{TRUE} (the default), redundancies are eliminated relative to \code{full.ct(x)}, otherwise relative to \code{x}.
}
  \item{simplify2constant}{
Logical; if \code{TRUE} (the default), a tautologous or contradictory \code{cond} is reduced to a constant \code{"1"} or \code{"0"}, respectively. If \code{FALSE}, a minimal tautology or contradiction, i.e. \code{"A+a"} or \code{"A*a"}, will result.
}
  \item{verbose}{
Logical; if TRUE, the reduction process will be traced in the console. 
}
  \item{maxiter}{
Maximal number of iterations.  
}
}

\details{

\code{ereduce} and \code{rreduce} eliminate conjuncts and disjuncts from a DNF \code{cond} as long as the consistency and coverage of \code{cond} in data \code{x} does not change, that is, as long as the result of \code{condition(cond, x)} remains the same. The only required argument is \code{cond}. If \code{x} is not provided, redundancies are eliminated relative to \code{full.ct(cond)}. If \code{x} is provided and \code{full = TRUE}, redundancies are eliminated relative to \code{full.ct(x)}. If \code{x} is provided and \code{full = FALSE}, redundancies are eliminated relative to \code{x}.

While \code{ereduce} generates all redundancy-free forms of \code{cond}, \code{rreduce} only returns one randomly chosen one. \code{rreduce} is faster than \code{ereduce}, but often incomplete. \code{ereduce}, in a nutshell, searches for minimal hitting sets in \code{cond} preventing \code{cond} from being false in data \code{x}. %the negation of \code{cond} from becoming true.

}

\value{
Redundancy-free disjunctive normal form (DNF).
}

\seealso{
\code{\link{full.ct}}, \code{\link{conCovOpt}}, \code{\link{DNFbuild}}.
}


\examples{
# Logical redundancies.
cond1 <- "A*b + a*B + A*C + B*C"
ereduce(cond1)
rreduce(cond1)
cond2 <- "A*b + a*B + A*B + a*b"
ereduce(cond2)
ereduce(cond2, simplify2constant = FALSE)


# Redundancy elimination relative to simulated cs data.
dat1 <- data.frame(
  A = c(0, 0, 0, 0, 1, 1, 0, 1), 
  B = c(0, 1, 0, 1, 1, 0, 0, 0), 
  C = c(1, 1, 0, 1, 1, 0, 1, 1), 
  D = c(0, 0, 0, 0, 0, 1, 1, 1))
cco1 <- conCovOpt(dat1, "D")
best1 <- selectMax(cco1)
formula1 <- DNFbuild(best1, outcome = "D", reduce = FALSE)
# ereduce
ereduce(formula1, dat1, full = FALSE)
# rreduce
rreduce(formula1, dat1, full = FALSE)
rreduce(formula1, dat1, full = FALSE, verbose = TRUE)


# Redundancy elimination relative to simulated mv data.
dat2 <- data.frame(
  A = c(3,2,1,1,2,3,2,2,2,1,1,2,3,2,2,2,1,2,3,3,3,1,1,1,3,1,2,1,2,3,3,2,2,2,1,2,2,3,2,1,2,1,3,3),
  B = c(1,2,3,2,1,1,2,1,2,2,3,1,1,1,2,3,1,3,3,3,1,1,3,2,2,1,1,3,3,2,3,1,2,1,2,2,1,1,2,2,3,3,3,3),
  C = c(1,3,3,3,1,1,1,2,2,3,3,1,1,2,2,2,3,1,1,2,1,2,2,3,3,1,2,2,2,3,2,1,1,2,2,2,1,1,1,2,2,1,1,2),
  D = c(3,1,2,2,1,1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,1,1,1,1,1,2,2,2,2,2,3,1,1,1,1,1,2,2,2,2,2,3,3,3),
  E = c(3,2,2,3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3,3)
)
cco2 <- conCovOpt(dat2, "D=3", type="mv")
best2 <- selectMax(cco2)
formula2 <- DNFbuild(best2, outcome = "D=3", reduce = FALSE)
# ereduce
ereduce(formula2, mvct(dat2), full = FALSE)
# rreduce
rreduce(formula2, mvct(dat2), full = FALSE)


# Any Boolean expressions.
cond <- "!(A*B*C)*!(a*b*c)" # or "A + B*!(D + e) <-> C" 
x <- selectCases(cond) 
cond <- cna:::getCond(x) # Returns a DNF equivalent to cond, but with many redundancies.
ereduce(cond)
rreduce(cond)
}
