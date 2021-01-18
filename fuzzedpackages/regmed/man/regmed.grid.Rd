\name{regmed.grid}
\alias{regmed.grid}
\alias{print.regmed.grid}
\title{
Regularized mediation models over a vector grid of lambda penalty values.
}
\description{
Fit regularized mediation models over a vector grid of lambda penalty
values. Structural equation models for analysis of multiple mediators
are extended by creating a sparse group lasso penalized model such that
the penalty considers the natural groupings of the pair of parameters
that determine mediation, as well as encourages sparseness of the model
parameters. The model is x[alpha] -> mediator -> [beta]outcome, where
alpha and beta are the parameters for the indirect effect of x on y,
through the mediator. The model also allows a direct effect of x on y: x[delta] -> y.
}
\usage{
regmed.grid(x, mediator, y, lambda.vec, frac.lasso=0.8, max.outer=5000,
max.inner=100, x.std=TRUE, med.std=TRUE, step.multiplier = 0.5,
wt.delta = .5, print.iter=FALSE)
}
\arguments{
  \item{x}{
vector representing "exposure" variable (sometimes called instrumental variable) 
}
  \item{mediator}{
matrix of mediators, rows are observations, columns are different mediators 
}
  \item{y}{
vector representing outcome
}
  \item{lambda.vec}{
vector of lambda penalty parameters
}
  \item{frac.lasso}{
fraction of penalty (lambda) that is allocated to L1 penalty
(lasso). The remaining fraction, (1-frac.lasso) is allocated to
group-lasso penalty, where the group is the pair of parameters alpha and
beta that determine mediation (x [alpha] -> mediator -> [beta] y).
}
  \item{max.outer}{
maximum number of outer loop iterations. The outer loop cycles over
several inner loops.
}
  \item{max.inner}{
maximum number of iterations for each inner loop. There is an inner loop
for each pair of alpha-beta parameters for each mediator, an inner loop
for direct effect (delta), and inner loops for residual variances for x
and for y.
}
  \item{x.std}{
logical (TRUE/FALSE) whether to standardize x by dividing by standard
devation of x. Note that x will be centered on its mean.
}
  \item{med.std}{
logical (TRUE/FALSE) whether to standardize mediators by dividing each mediator by its standard
deviation. Note that mediators will be centered on their means.
}
\item{step.multiplier}{
a value between 0 and 1 for backtracking, to shrink step size. Value of
0.5 is typical default.
}
\item{wt.delta}{
  a weight >=0 for how much weight should be given to shrinking delta
  parameter, by penalaty lambda*wt.delta.
}
  \item{print.iter}{
print iteration history during fitting routine 
}
}
\details{
Altough outcome y is not required to be scaled by its standard
deviation, it can be beneficial to scale y. This helps with setting
range of lambda penalty parameters, because when all x, y, and mediators
are scaled, it is reasonable to consider lambda values within the range
of 0 to 1. See reference for details of algorithm.
}
\value{
regmed.grid object
}
\references{
Schaid, DJ, Sinnwell JP. (2020) Penalized Models for Analysis of
Multiple Mediators. To appear in Genetic Epidemiology
}
\author{
 Dan Schaid, Greg Jenkins, Jason Sinnwell
}
\seealso{
\code{\link{plot.regmed.grid}}
\code{\link{regmed.fit}}
}
\examples{
 data(regmed_example)

 y <- regmed_example$y
 x <- regmed_example$x
 med <- regmed_example[, -c(1,2)]
 fit.grid <- regmed.grid(x, med, y, lambda.vec= c(seq(from=1, to=0, by = -.1)), frac.lasso=.8)
 fit.grid
}

\keyword{ models }
