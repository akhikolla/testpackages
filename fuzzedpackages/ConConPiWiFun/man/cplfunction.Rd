\name{cplfunction}
\alias{cplfunction}
\alias{Suml}
\alias{cplfunction}
\alias{Rcpp_cplfunction-class}
\alias{plot-methods}
\alias{plot,ANY-method}
\alias{plot,Rcpp_cplfunction-method}
\alias{InfConvl}
%\alias{show-methods}
%\alias{show,ANY-method}
\alias{show,Rcpp_cplfunction-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
This class implements continuous convex piecewise linear functions
}
\description{
This includes functions that are ccpl (resp. ccpq) on a convex set (i.e. an interval or a point) and infinite out of the domain. These functions can be very usefull for a large class of optimisation problems. Efficient manipulation (such as log(N) insertion) of such data structure is obtained with map standard template library of C++ (that hides balanced trees). This package is a wrapper on such a class based on Rcpp modules. 
}

\author{
Robin Girard
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 to See Also as \code{\link{cplfunction}},
}
\examples{
##
#Construction of a piecewise linear function
##


Slopes=c(-1,2,Inf) # increasing ! convexity is required
Breakpoints=c(-Inf,2,4) # increasing. length is number of slopes +1
FirstNonInfBreakpointVal=3
CCPWLfunc1=new(cplfunction,Slopes,Breakpoints,FirstNonInfBreakpointVal)
plot(CCPWLfunc1) #visualisation method


###Etoile transformation (legendre transform of f)
# Changes f no return value
CCPWLfunc1$Etoile()
plot(CCPWLfunc1)  #if f = CCPWLfunc1 CCPWLfunc1 becomes is f^*(y) =inf_x {xy-f(x)}
CCPWLfunc1$Etoile()
plot(CCPWLfunc1)   ## (f^*)^* is f !


###Squeeze function 
# Changes f, no return value
left=-Inf; right=3
CCPWLfunc1$Squeeze(left,right) # CCPWLfunc1 is now infinite (or not definite) out of [left,right]
# i.e. all breakpoints out of [left,right] removed


###Swap function 
# Changes f  no return value !
y=2;
CCPWLfunc1$Swap(y)
plot(CCPWLfunc1); #now f = CCPWLfunc1 is replaced by x -> f(y-x) 


### Sum function (uses fast insertion) do not affect operands 
CCPWLfunc1=new(cplfunction,c(-1,2,Inf) ,c(-Inf,2,4),0)
CCPWLfunc2=new(cplfunction,c(-1,2,Inf),c(-Inf,1,3),0)
CCPWLfunc1plus2=Suml(CCPWLfunc1,CCPWLfunc2)
CCPWLfunc1plus2


par(mfrow=c(1,3))
plot(CCPWLfunc2,col='red'); 
plot(CCPWLfunc1,col='blue'); 
plot(CCPWLfunc1plus2);


rm(list=ls())
gc()

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
