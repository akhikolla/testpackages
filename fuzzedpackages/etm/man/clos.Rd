\name{clos}
\alias{clos}
\alias{clos.etm}
\alias{clos.msfit}
\title{Change in Length of Stay}
\description{
  The function estimates the expected change in length of stay (LOS)
  associated with an intermediate event.
}
\usage{
clos(x, aw, ratio, ...)
\S3method{clos}{etm}(x, aw = FALSE, ratio = FALSE, ...)
\S3method{clos}{msfit}(x, aw = FALSE, ratio = FALSE, cox_model, ...)
}
\arguments{
  \item{x}{An object of class \code{etm}. Argument \code{delta.na} in
    \code{\link{etm}} must be set to \code{TRUE} in order to use this
    function.}
  \item{aw}{Logical. Whether to compute the expected change of LOS using
    alternative weighting. Default is \code{FALSE}.}
  \item{ratio}{Logical. Compute the ratio of the expected length-of-stay
    given instermediate event status instead of a difference. Default
    value is \code{FALSE}}
  \item{cox_model}{TODO}
  \item{...}{Further arguments}
}

\details{
  The approach for evaluating the impact of an intermediate
  event on the expected change in length of stay is based on Schulgen
  and Schumacher (1996). They suggested to consider the difference of
  the expected subsequent stay given infectious status at time s. 

  Extensions to the methods of Schulgen and Schumacher and the earlier
  implementation in the \pkg{changeLOS} include the possibility to
  compute the extra length of stay both for competing endpoints and the
  more simple case of one absorbing state, as well as the possibility to
  compute this quantity for left-truncated data. 
}

\value{
  An object of class \code{clos.etm} with the following components:
  \item{e.phi}{Change in length of stay}
  \item{phi.case}{Estimates of \eqn{E(\mbox{LOS} | X_s =
      \mbox{intermediate event})}{E(LOS | X_s = intermediate event)} for
    all observed transition times \eqn{s}{s}, where
    \eqn{X_s}{X_s}denotes the state by time \eqn{s}{s}}
  \item{phi.control}{Estimates of \eqn{E(\mbox{LOS} | X_s =
      \mbox{initial state})}{E(LOS|X_s = initial state)} for
    all observed transition times \eqn{s}{s}.}
  \item{e.phi2}{Weighted average of the difference between
    \code{phi2.case} and \code{phi2.control}.}
  \item{phi2.case}{Estimates of \eqn{E(\mbox{LOS}
      \mathbf{1}(X_{\mbox{LOS}} = \mbox{discharge}) | X_s =
      \mbox{intermediate event})}{E(LOS \strong{1}(X_LOS = discharge)|X_s =
      intermediate event)}, where \eqn{\mathbf{1}}{\strong{1}} denotes
    the indicator function.}
  \item{phi2.control}{\eqn{E(\mbox{LOS}
      \mathbf{1}(X_{\mbox{LOS}} = \mbox{discharge}) | X_s =
      \mbox{initial state})}{E(LOS \strong{1}(X_LOS = discharge)|X_s =
      initial state)}.}
  \item{e.phi3}{Weighted average of the difference between
    \code{phi3.case} and \code{phi3.control}.}
  \item{phi3.case}{Estimates of \eqn{E(\mbox{LOS}
      \mathbf{1}(X_{\mbox{LOS}} = \mbox{death}) | X_s =
      \mbox{intermediate event})}{E(LOS \strong{1}(X_LOS = death)|X_s =
      intermediate event)}.}
  \item{phi3.control}{\eqn{E(\mbox{LOS}
      \mathbf{1}(X_{\mbox{LOS}} = \mbox{death}) | X_s =
      \mbox{initial state})}{E(LOS \strong{1}(X_LOS = death)|X_s =
      initial state)}.}
  \item{weights}{Weights used to compute the weighted averages.}
  \item{w.time}{Times at which the weights are computed.}
  \item{time}{All transition times.}
  \item{e.phi.weights.1}{Expected change in LOS using \code{weights.1}}
  \item{e.phi.weights.other}{Expected change in LOS using
    \code{weights.other}}
  \item{weights.1}{Weights corresponding to the conditional waiting
    time in the intial state given one experiences the intermediate event.}
  \item{weights.other}{Weights corresponding to the conditional
    waiting time given one does not experience the intermediate event.} 
}

\references{
  G Schulgen and M Schumacher (1996). Estimation of prolongation of
  hospital stay attributable to nosocomial infections.  \emph{Lifetime
  Data Analysis} 2, 219-240.

  J Beyersmann, P Gastmeier, H Grundmann, S Baerwolf, C Geffers,
  M Behnke, H Rueden, and M Schumacher (2006). Use of Multistate
  Models to Assess Prolongation of Intensive Care Unit Stay Due to
  Nosocomial Infection. \emph{Infection Control and Hospital
    Epidemiology} 27, 493-499.

  Allignol A, Schumacher M, Beyersmann J: Estimating summary functionals
  in multistate models with an application to hospital infection
  data. \emph{Computation Stat}, 2011; 26: 181-197.

  M Wrangler, J Beyersmann and M Schumacher (2006). changeLOS: An
  R-package for change in length of hospital stay based on the
  Aalen-Johansen estimator. \emph{R News} 6(2), 31--35.
}
  
\author{Arthur Allignol \email{arthur.allignol@gmail.com},
  Matthias Wangler, Jan Beyersmann}
\seealso{\code{\link{etm}}}
\examples{
data(los.data)

## putting los.data in the long format
my.observ <- prepare.los.data(x=los.data)

tra <- matrix(FALSE, 4, 4)
tra[1, 2:4] <- TRUE
tra[2, 3:4] <- TRUE

tr.prob <- etm(my.observ, c("0","1","2","3"), tra, NULL, 0)

cLOS <- etm::clos(tr.prob)
plot(cLOS)


### Compute bootstrapped SE

## function that performs the bootstrap
## nboot: number of bootstrap samples. Other arguments are as in etm()
boot.clos <- function(data, state.names, tra, cens.name, s = 0, nboot) {
    res <- double(nboot)
    for (i in seq_len(nboot)) {
        index <- sample(unique(data$id), replace = TRUE)
        inds <- new.id <- NULL
        for (j in seq_along(index)){
            ind <- which(data$id == index[j])
            new.id <- c(new.id, rep(j, length(ind)))
            inds <- c(inds, ind)
        }
        dboot <- cbind(data[inds, ], new.id)
        dboot[, which(names(dboot) == "id")]
        dboot$id <- dboot$new.id
        tr.prob <- etm(dboot, state.names, tra, cens.name, s, cova = FALSE)
        res[i] <- etm::clos(tr.prob)$e.phi
    }
    res
}

## bootstrap
se <- sqrt(var(boot.clos(my.observ, c("0","1","2","3"), tra, NULL, 0,
                         nboot = 10)))
}

\keyword{survival}
