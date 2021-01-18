\name{splmmTuning}
\alias{splmmTuning}
\title{Tuning funtion of \code{'splmm'} object}
\usage{
splmmTuning(x, y, z, grp, lam1.seq, lam2.seq, nonpen.b=1,nonpen.L=1,
                          penalty.b=c("lasso","scad"),
                          penalty.L=c("lasso","scad"),
                          CovOpt=c("nlminb","optimize"),
                          standardize=TRUE,control=splmmControl())}
                          
\arguments{
  \item{x}{matrix of dimension N x p including the fixed-effects covariables. An intercept has to be included in the first column as (1,...,1).}
  \item{y}{response variable of length N.}
  \item{z}{random effects matrix of dimension N x q. It has to be a matrix, even if q=1.}
  \item{grp}{grouping variable of length N}
  \item{lam1.seq}{a grid of regularization parameter for fixed effects penalization, could be a scalar if no need to tune. }
  \item{lam2.seq}{a grid of regularization parameter for random effects penalization, could be a scalar if no need to tune. }
  \item{nonpen.b}{Index of indices of fixed effects not penalized. The default value is 1, which means the fixed intercept is not penalized. }
  \item{nonpen.L}{Index of indices of random effects not penalized. The default value is 1, which means the random intercept is not penalized. }
  \item{penalty.b}{The penalty method for fixed effects penalization. Currently available options include LASSO penalty and SCAD penalty. }
  \item{penalty.L}{The penalty method for fixed effects penalization. Currently available options include LASSO penalty and SCAD penalty. }
  \item{CovOpt}{which optimization routine should be used for updating the variance parameter. The available options include optimize and nlminb. nlminb uses the estimate                     of the last iteration as a starting value. nlminb is faster if there are many Gauss-Seidel iterations.}
  \item{standardize}{A logical parameter specifying whether the fixed effects matrix x and random effects matrix z should  be standardized such that each column has mean 0 and standard deviation 1. The default value is \code{TRUE}.}
  \item{control}{control parameters for the algorithm and the Armijo Rule, see \code{'splmmControl'} for the details. }
}

\description{
This function fits \code{'splmm'} function over grids of lambda1 and/or lambda2 and determine the best fit model based on model selection information criterion. 
The function takes a scalar or a grid of lambda1 and/or lambda2 and determine the optimal tuning parameter value for the best model fit. If both lambda1 and lambda2 are inputted as scalars, an \code{'splmm'} object is returned; if either or both lambda1 and lambda2 are inputted as grids, an \code{'splmm.tuning'} object is returned. Currently the model selection criterion include AIC and BIC, and BIC is used to determine the optimal model. 
}


\value{A \code{'splmm.tuning'} object is returned, for which \code{plot} method exist.


  \item{lam1.seq}{lambda1 grid used for tuning. Only available when lambda1 is inputted as a vector.}
  \item{lam2.seq}{lambda2 grid used for tuning. Only available when lambda2 is inputted as a vector.}
  \item{BIC.lam1}{A vector of BIC values of splmm models fitting over a lambda1 grid.}
  \item{BIC.lam2}{A vector of BIC values of splmm models fitting over a lambda2 grid.}
  \item{fit.BIC}{An array of BIC values of splmm models fitting over lambda 1 grid x lambda2 grid.}
  \item{AIC.lam1}{A vector of AIC values of splmm models fitting over a lambda1 grid.}
  \item{AIC.lam2}{A vector of AIC values of splmm models fitting over a lambda2 grid.}
  \item{fit.AIC}{An array of AIC values of splmm models fitting over lambda 1 grid x lambda2 grid.}
  \item{BBIC.lam1}{A vector of BBIC values of splmm models fitting over a lambda1 grid.}
  \item{BBIC.lam2}{A vector of BBIC values of splmm models fitting over a lambda2 grid.}
  \item{fit.BBIC}{An array of BBIC values of splmm models fitting over lambda 1 grid x lambda2 grid.}
  \item{EBIC.lam1}{A vector of EBIC values of splmm models fitting over a lambda1 grid.}
  \item{EBIC.lam2}{A vector of EBIC values of splmm models fitting over a lambda2 grid.}
  \item{fit.EBIC}{An array of EBIC values of splmm models fitting over lambda 1 grid x lambda2 grid.}
  \item{min.BIC}{The minimum BIC value from tuning over a grid. This is only available when either lambda1 or lambda2 is a scalar. }
  \item{min.AIC}{The minimum AIC value from tuning over a grid. This is only available when either lambda1 or lambda2 is a scalar. }
  \item{min.BBIC}{The minimum BBIC value from tuning over a grid. This is only available when either lambda1 or lambda2 is a scalar. }
  \item{min.EBIC}{The minimum EBIC value from tuning over a grid. This is only available when either lambda1 or lambda2 is a scalar. }
  \item{best.model}{The index of the optimal model. This is only available when either lambda1 or lambda2 is a scalar. }
  \item{best.fit}{The optimal model chosen by the minimum BIC as an \code{splmm} object. }
  \item{min.lam1}{lambda1 value that results in the optimal model. This is only available when input lambda1 is a vector.}
  \item{min.lam2}{lambda2 value that results in the optimal model. This is only available when input lambda2 is a vector.}
  \item{lam1.tuning}{A \code{logical} parameter specifying if tuning is performed over lamdbda1 grid. \code{lam1.tuning=TRUE} if input lambda1 is a vector. }
  \item{lam2.tuning}{A \code{logical} parameter specifying if tuning is performed over lamdbda2 grid. \code{lam1.tuning=TRUE} if input lambda2 is a vector. }

}
\examples{

\donttest{
data(cognitive)

x <- model.matrix(ravens ~schoolid+treatment+year+sex+age_at_time0
                  +height+weight+head_circ+ses+mom_read+mom_write
                  +mom_edu, cognitive)
z <- x

## Tuning over lambda1 grid
lam1 = seq(0.1,0.5,0.1)
lam2 = 0.1
fit1 <-splmmTuning(x=x,y=cognitive$ravens,z=z,grp=cognitive$id,lam1.seq=lam1,
lam2.seq=lam2,penalty.b="scad", penalty.L="scad")
plot.splmm(fit1)

## Tuning over lambda2 grid
lam1 = 0.1
lam2 = seq(0.1,0.5,0.1)
fit2 <-splmmTuning(x=x,y=cognitive$ravens,z=z,grp=cognitive$id,lam1.seq=lam1,
lam2.seq=lam2,penalty.b="scad", penalty.L="scad")
plot.splmm(fit2)

## Tuning over both lambda1 and lambda2 grid
lam1 = seq(0.1,0.5,0.2)
lam2 = seq(0.1,0.5,0.2)
fit3 <-splmmTuning(x=x,y=cognitive$ravens,z=z,grp=cognitive$id,lam1.seq=lam1,
lam2.seq=lam2,penalty.b="scad", penalty.L="scad")
plot.splmm(fit3)
}

}


