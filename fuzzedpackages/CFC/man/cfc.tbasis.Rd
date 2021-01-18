\name{cfc.tbasis}
\alias{cfc.tbasis}
\title{
Cause-specific competing-risk survival analysis in time denomination
}
\description{
Constructing cumulative incidence and event-free probability functions from cause-specific survival probabilities evaluated at fixed time points.
}
\usage{
cfc.tbasis(p1, p2, unity.tol = 1e-06, diff.tol = 0.01,
  diff.tol.policy = c("mean", "all"), check = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p1}{Multi-dimensional array containing survival probabilities for cause 1 (i.e. exponential of the negative integral of hazard function). First dimension must correspond to time points at which probabilities are calculated. Elements with same time, but distributed in the space of remaining dimensions, are treated independently. These diemensions can correspond, e.g., to observations or samples (in Bayesian frameworks).}
  \item{p2}{Multi-dimensional array containing survival probabilities for cause 2. See note for \code{p1}.}
  \item{unity.tol}{Tolerance for difference of survival probabilities from 1.0 at time=0.0, which is the first 'row' of arrays \code{p1} and \code{p2}. For example, for two-dimensional arrays, we need \code{all(abs(p1[1,]-1.0) < unity.tol)}, and a similar condition for \code{p2}.}
  \item{diff.tol}{Tolerance for change in survival probabilities from one time point to the next. Large changes lead to higher errors during numerical integration.}
  \item{diff.tol.policy}{If \code{"mean"}, then average change in survival probabilities are compared to \code{diff.tol}. If \code{"all"}, each values is compared. The latter is more strict.}
  \item{check}{Boolean flag indicating whether or not to check probability arrays for validity. Current validity checks are: 1) ensuring all probabilities are between 0.0 and 1.0, 2) all probabilities at time=0.0 are equal to 1.0 (see \code{unity.tol}), 3) No changes in probabilities from one time point to next are too large (see \code{diff.tol}). Dimensional consistency between \code{p1} and \code{p2} is always checked regardless of the value of this flag.}
}
\details{
Assuming one-dimensional \code{p1} and \code{p2} for clarity, the algorithm calculates cumulative incidence function for cuase 1 using a recursive formula: \code{ci1[n+1] = ci1[n] + dci1[n]}, where \code{dci1[n] = 0.5*(p2[n] + p2[n+1])*(p1[n] - p1[n+1])}. The increment in cumulative incidence function for cause 2 is similarly calculated, \code{dci2[n] = 0.5*(p1[n] + p1[n+1])*(p2[n] - p2[n+1])}. These equations guarantee that \code{dci1[n] + dci2[n] = p1[n]*p2[n] - p1[n+1]*p2[n+1]}. Event-free probability is simply calculated as code{efp[n] = p1[n]*p2[n]}. Taken together, this numerical integration ensures that \code{efp[n+1] - efp[n] + dci1[n] + dci2[n] = 0}.
}
\value{
If \code{p1} and \code{p2} are one-dimensional arrays (i.e. vectors), a matrix with columns named \code{"ci1"}, \code{"ci2"} and \code{"efp"} is returned, representing the cummulative incidence functions for cause 1 and cause 2 and the event-free probability, evaluated at same time points as \code{p1} and \code{p2} are provided. If \code{p1} and \code{p2} are multi-dimensional arrays, a list is returned with elements \code{"ci1"}, \code{"ci2"} and \code{"efp"}, each one with the same interpretation, and all of the same dimensions as \code{p1} and \code{p2}.
}
\references{
Mahani A.S. and Sharabiani M.T.A. (2019). Bayesian, and Non-Bayesian, Cause-Specific Competing-Risk Analysis for Parametric and Nonparametric Survival Functions: The R Package CFC. Journal of Statistical Software, 89(9), 1-29. doi:10.18637/jss.v089.i09

Prentice et al (1978). The analysis of failure times in the presence of competing risks. Biometrics, 541-554.
}
\author{
Mansour T.A. Sharabiani, Alireza S. Mahani
}
\note{
The integration algorithm described above does not require knowledge of time step. (Alternatively, using hazard functions for integration would have required specification of time step.) Since \code{p1} and \code{p2} are integrals (followed by exponentiation) of cause-specific hazard functions, using them directly adds to robustness of numerical integration and avoids error accumulation. The returned cumulative incidence and event-free probabilities correspond to the same time points assumed for input cause-specific probabilities.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{cfc.pbasis}}
}
\examples{
\dontrun{

# prepare data for cause-specific competing-risk analysis
data(bmt)
bmt$status1 <- 1*(bmt$cause==1)
bmt$status2 <- 1*(bmt$cause==2)
f1 <- Surv(time, status1) ~ platelet + age + tcell
f2 <- Surv(time, status2) ~ platelet + age + tcell

# sample-based bayesian weibull regression
library(BSGW)
reg1 <- bsgw(f1, bmt, ordweib = TRUE, control = bsgw.control(iter = 500, burnin = 100, nskip = 50))
reg2 <- bsgw(f2, bmt, ordweib = TRUE, control = bsgw.control(iter = 500, burnin = 100, nskip = 50))

# prediction on a uniform grid of 100 time points
# (use first 50 observations for speed)
pred1 <- predict(reg1, newdata = bmt[1:50,], tvec = 100)
pred2 <- predict(reg2, newdata = bmt[1:50,], tvec = 100)

# permuting dimensions of survival objects to conform with cfc
S1 <- aperm(pred1$smp$S, c(2,1,3))
S2 <- aperm(pred2$smp$S, c(2,1,3))

# cause-specific competing risk analysis - time mode
my.cfc <- cfc.tbasis(S1, S2)

# calculating averages across observations (e.g. patients in the study)
my.summ <- summary(my.cfc, MARGIN = c(1,2))

# plotting mean CI and event-free functions
# as well as their sampled-based confidence intervals
plot(my.summ, t = pred1$tvec)


}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
