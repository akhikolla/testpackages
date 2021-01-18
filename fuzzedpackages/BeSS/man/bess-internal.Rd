\name{BeSS-internal}
\title{Internal BeSS functions}
\alias{plot_solution}
\alias{plot_loss}
\alias{bess_lm}
\alias{gbess_lm}
\alias{get_A}
\alias{getcox_A}
\alias{gget_A}
\alias{bess.lm}
\alias{bess.glm}
\alias{bess.cox}
\alias{gbess.lm}
\alias{gbess.glm}
\alias{gbess.cox}
\alias{EigenR}
\alias{uniqueR}
\description{Internal BeSS functions}
\usage{
bess.lm(x, y, beta0, s, max.steps=20, factor = NULL,
        weights=rep(1,nrow(x)), normalize=FALSE)
bess.glm(x, y, beta0, intercept=0, s, max.steps=20,
         glm.max=1e6, factor = NULL,
         weights=rep(1,nrow(x)), normalize=FALSE)
bess.cox(x, y, beta0, s, cox.max=20, max.steps=20, factor=NULL,
         weights=rep(1,nrow(x)), normalize=FALSE)
gbess.lm(x, y, Gi, beta0, s, max.steps = 20,
         weights=rep(1,nrow(x)), normalize=FALSE)
gbess.glm(x, y, Gi, beta0, intercept=0, s, max.steps = 10, glm.max=1e6,
          weights=rep(1,nrow(x)), normalize=FALSE)	
gbess.cox(x, y, Gi, beta0, s, cox.max=20, max.steps=10,
          weights=rep(1,nrow(x)), normalize=FALSE)		  
}
\author{
  Canhong Wen, Aijun Zhang, Shijie Quan, and Xueqin Wang.
}
\references{
 Wen, C., Zhang, A., Quan, S. and Wang, X. (2020). BeSS: An R Package for Best Subset Selection in Linear, Logistic and Cox Proportional Hazards Models, \emph{Journal of Statistical Software}, Vol. 94(4).  doi:10.18637/jss.v094.i04.
}
\details{These are not intended for use by users. \code{bess.lm}fit a linear regression model. \code{bess.glm}fit a logistic model. \code{bess.cox}fit a cox model.}
\keyword{internal}
