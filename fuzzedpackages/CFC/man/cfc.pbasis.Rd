\name{cfc.pbasis}
\alias{cfc.pbasis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cause-specific competing-risk survival analysis in probability denomination
}
\description{
Constructing cumulative incidence and event-free probability functions from cause-specific survival times give for a fixed set of probabilities.
}
\usage{
cfc.pbasis(t1, t2, probs, unity.tol = 1e-06, diff.tol = 0.01,
  diff.tol.policy = c("all", "mean"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t1}{Multi-dimensional array containing survival times for cause 1 (i.e. exponential of the negative integral of hazard function). First dimension must correspond to probabilities at which times are calculated. Elements with same time, but distributed in the space of remaining dimensions, are treated independently. These diemensions can correspond, e.g., to observations or samples (in Bayesian frameworks). Survival times must be increasing along the first dimension.}
  \item{t2}{Multi-dimensional array containing survival times for cause 2. See note for \code{t1}.}
  \item{probs}{Probabilities for which survival times are provided in \code{t1} and \code{t2}. Must begin at 1.0 and be decreasing.}
  \item{unity.tol}{Tolerance for difference of survival probability from 1.0 at time=0.0. In other words, we require that \code{abs(probs[1] - 1.0) < unity.tol}.}
  \item{diff.tol}{Tolerance for change in survival probabilities from one time point to the next. Large changes lead to higher errors during numerical integration.}
  \item{diff.tol.policy}{If \code{"mean"}, then average change in survival probabilities are compared to \code{diff.tol}. If \code{"all"}, each values is compared. The latter is more strict.}
}
\details{
For each 'row' of \code{t1}, and \code{t2}, all elements are processed independently. To combine the survival curves from corresponding elements of \code{t1} and \code{t2}, we first form a 'comon denominator' time vector by combining the two time vectors and sorting the results (after removing duplicates). We limit the maximum value in the combined time vector to minimum of the the two maxima from each cause. Next, we use interpolation to find survival probabilities of each cause at all the time points in the combined time vector. Finally, we call the function \code{\link{cfc.tbasis}}.
}
\value{
If \code{t1} and \code{t2} are one-dimensional, a matrix with columns named \code{"time"}, \code{"ci1"}, \code{"ci2"} and \code{"efp"} is returned. For multi-dimensional arrays, a list is returned with one such matrix for each element of the consolidated dimension representing all but the first dimension of \code{t1} and \code{t2}.
}
\references{
Mahani A.S. and Sharabiani M.T.A. (2019). Bayesian, and Non-Bayesian, Cause-Specific Competing-Risk Analysis for Parametric and Nonparametric Survival Functions: The R Package CFC. Journal of Statistical Software, 89(9), 1-29. doi:10.18637/jss.v089.i09

Prentice et al (1978). The analysis of failure times in the presence of competing risks. Biometrics, 541-554.
}
\author{
Mansour T.A. Sharabiani, Alireza S. Mahani
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{cfc.tbasis}}
}
\examples{
\dontrun{

# prepare data for cause-specific competing-risk analysis
data(bmt)
bmt$status1 <- 1*(bmt$cause==1)
bmt$status2 <- 1*(bmt$cause==2)
f1 <- Surv(time, status1) ~ platelet + age + tcell
f2 <- Surv(time, status2) ~ platelet + age + tcell

# perform weibull regression on each cause independently
library(survival)
reg1 <- survreg(f1, bmt)
reg2 <- survreg(f2, bmt)

# predict times for given probabilities
# transpose predictions so that first dimension
# is time/probability (use first 50 observations for speed)
pvec <- seq(from=1.0, to = 0.1, length.out = 100)
pred1 <- t(predict(reg1, newdata = bmt[1:50,], p = 1-pvec, type = "quantile"))
pred2 <- t(predict(reg2, newdata = bmt[1:50,], p = 1-pvec, type = "quantile"))

# cause-specific competing risk analysis - probability mode
my.cfc <- cfc.pbasis(pred1, pred2, probs = pvec)

# calculating averages across observations (e.g. patients in the study)
my.summ <- summary(my.cfc)

# plotting average CI and event-free probability curves
plot(my.summ)

}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
