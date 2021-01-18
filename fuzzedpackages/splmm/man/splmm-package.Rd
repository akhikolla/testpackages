\name{splmm-package}
\alias{splmm-package}
\docType{package}
\title{\packageTitle{splmm}}
\description{\packageDescription{splmm}}
\details{
  The DESCRIPTION file: \packageDESCRIPTION{splmm}
  \packageIndices{splmm}

 Contains functions that fit linear mixed-effects models for high-dimensional data (p>>n) with penalty for both the fixed effects and random effects for variable selection.
}
\author{
  \packageAuthor{splmm}
  
  Maintainer: \packageMaintainer{splmm}
}
\references{

Luoying Yang PhD thesis

      SCHELLDORFER, J., BUHLMANN, P. and DE GEER, S.V. (2011), Estimation for High-Dimensional Linear Mixed-Effects Models Using L1-Penalization. Scandinavian Journal of Statistics, 38: 197-214. doi:10.1111/j.1467-9469.2011.00740.x
  
}
% Optionally other standard keywords, one per line,
% from the file KEYWORDS in the R documentation.
\keyword{package}

\examples{
\donttest{
## Use splmm on the Kenya school cognitive data set


data(cognitive)

x <- model.matrix(ravens ~ schoolid+treatment+relyear+sex+age_at_time0+height+weight+head_circ+ses1
      +readtest+writetest+yrsofsch, cognitive)
z <- model.matrix(ravens ~ schoolid+treatment+relyear, cognitive)

fit <- splmm(x=x,y=cognitive$ravens,z=z,grp=cognitive$id,lam1=0.01,
lam2=0.01,penalty.b="scad", penalty.L="scad")
summary(fit)
}

}
