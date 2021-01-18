\name{fwsim}
\alias{fwsim}
\alias{fwsim_fixed}
\alias{print.fwsim}
\alias{summary.fwsim}
\alias{plot.fwsim}
\title{Fisher-Wright Population Simulation}
\description{This package provides tools to simulate a population under the Fisher-Wright model with a stepwise neutral mutation process on \eqn{r} loci, where mutations on loci happen independently. The population sizes are either fixed (traditional/original Fisher-Wright model) or random Poisson distributed with exponential growth supported. Intermediate generations can be saved in order to study e.g. drift.

For stochastic population sizes:
Model described in detail at \url{http://arxiv.org/abs/1210.1773}. Let \eqn{M} be the population size at generation \eqn{i} and \eqn{N} the population size at generation \eqn{i + 1}.
Then we assume that \eqn{N} conditionally on \eqn{M} is \eqn{\mbox{Poisson}(\alpha M)}{Poisson(\alpha*M)} distributed for \eqn{\alpha > 0} (\eqn{\alpha > 1} gives expected growth and \eqn{0 < \alpha < 1} gives expected decrease).

For each haplotype \eqn{x} occuring \eqn{m} times in the \eqn{i}'th generation, the number
of children \eqn{n} is \eqn{\mbox{Poisson}(\alpha m)}{Poisson(\alpha*m)} distributed independently of other haplotypes. 
It then follows that the sum of the number of haplotypes follows a \eqn{\mbox{Poisson}(\alpha M)}{Poisson(\alpha*M)} distribution (as just stated in the previous paragraph) and that \eqn{n} conditionally on \eqn{N} follows a \eqn{\mbox{Binomial}(N, m/M)}{Binomial(N, m/M)} as expected.

The mutation model can be e.g. the stepwise neutral mutation model. See \code{\link{init_mutmodel}} for details.
}
\usage{
fwsim(G, H0, N0, mutmodel, alpha = 1.0, SNP = FALSE, 
  save_generations = NULL, progress = TRUE, trace = FALSE, ensure_children = FALSE, ...)
fwsim_fixed(G, H0, N0, mutmodel, SNP = FALSE, 
  save_generations = NULL, progress = TRUE, trace = FALSE, ...)
\method{print}{fwsim}(x, ...)
\method{summary}{fwsim}(object, ...)
\method{plot}{fwsim}(x, which = 1L, ...)
}
\arguments{
\item{G}{ number of generations to evolve (integer, remember postfix L).}
\item{H0}{ haplotypes of the initial population. Must be a vector or matrix (if more than one initial haplotype). The number of loci is the length or number of columns of \code{H0}.}
\item{N0}{ count of the \code{H0} haplotypes. The i'th element is the count of the haplotype \code{H0[i, ]}. \code{sum(N0)} is the size of initial population.}
\item{mutmodel}{ a \code{mutmodel} object created with \code{\link{init_mutmodel}}. Alternatively, a numeric vector of length r of mutation probabilities (this will create a stepwise mutation model with r loci divide the mutation probabilities evenly between upwards and downwards mutation).}
\item{alpha}{ vector of length 1 or \code{G} of growth factors (1 correspond to expected constant population size). If length 1, the value is reused in creating a vector of length \code{G}.}
\item{SNP}{ to make alleles modulus 2 to immitate SNPs.}
\item{save_generations}{ to save intermediate populations. \code{NULL} means that no intermediate population will be saved. Else, a vector of the generation numbers to save.}
\item{progress}{ whether to print progress of the evolution.}
\item{trace}{ whether to print trace of the evolution (more verbose than \code{progress}).}
\item{ensure_children}{ Ensures that every generation has at least one child; implemented by getting \eqn{\mbox{Poisson}(\alpha M) + 1}{Poisson(\alpha*M) + 1} children.}
\item{x}{A \code{fwsim} object.}
\item{object}{A \code{fwsim} object.}
\item{which}{A number specifying the plot (currently only 1: the actual population sizes vs the expected sizes).}
\item{\dots}{ not used.}
}
%\details{
%}
\value{
A \code{fwsim} object with elements
\item{pars}{ the parameters used for the simulation}
\item{saved_populations}{ a list of haplotypes in the intermediate populations}
\item{population}{ haplotypes in the end population after \code{G} generations}
\item{pop_sizes}{ the population size for each generation}
\item{expected_pop_sizes}{ the expected population size for each generation}
}
\author{Mikkel Meyer Andersen <mikl@math.aau.dk> and Poul Svante Eriksen}
\examples{
# SMM (stepwise mutation model) example
set.seed(1)
fit <- fwsim(G = 100L, H0 = c(0L, 0L, 0L), N0 = 10000L, 
  mutmodel = c(Loc1 = 0.001, Loc2 = 0.002, Loc3 = 0.003))
summary(fit)
fit

# SMM (stepwise mutation model) example
H0 <- matrix(c(0L, 0L, 0L), 1L, 3L, byrow = TRUE)

mutmodel <- init_mutmodel(modeltype = 1L, 
  mutpars = matrix(c(c(0.003, 0.001), rep(0.004, 2), rep(0.001, 2)), 
                   ncol = 3, 
                   dimnames = list(NULL, c("DYS19", "DYS389I", "DYS391"))))
mutmodel

set.seed(1)
fit <- fwsim(G = 100L, H0 = H0, N0 = 10000L, mutmodel = mutmodel)

xtabs(N ~ DYS19 + DYS389I, fit$population)
plot(1L:fit$pars$G, fit$pop_sizes, type = "l", 
  ylim = range(range(fit$pop_sizes), range(fit$expected_pop_sizes)))
points(1L:fit$pars$G, fit$expected_pop_sizes, type = "l", col = "red")

set.seed(1)
fit_fixed <- fwsim_fixed(G = 100L, H0 = H0, N0 = 10000L, mutmodel = mutmodel)

# LMM (logistic mutation model) example
mutpars.locus1 <- c(0.149,   2.08,    18.3,   0.149,   0.374,   27.4) # DYS19
mutpars.locus2 <- c(0.500,   1.18,    18.0,   0.500,   0.0183,  349)  # DYS389I
mutpars.locus3 <- c(0.0163,  17.7,    11.1,   0.0163,  0.592,   14.1) # DYS391
mutpars <- matrix(c(mutpars.locus1, mutpars.locus2, mutpars.locus3), ncol = 3)
colnames(mutpars) <- c("DYS19", "DYS389I", "DYS391")
mutmodel <- init_mutmodel(modeltype = 2L, mutpars = mutpars)
mutmodel

set.seed(1)
H0_LMM <- matrix(c(15L, 13L, 10L), 1L, 3L, byrow = TRUE)
fit_LMM <- fwsim(G = 100L, H0 = H0_LMM, N0 = 10000L, mutmodel = mutmodel)
xtabs(N ~ DYS19 + DYS389I, fit_LMM$population)
}
\keyword{Fisher-Wright}

