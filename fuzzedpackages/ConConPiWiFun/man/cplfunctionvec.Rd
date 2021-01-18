\name{cplfunctionvec}
\alias{cplfunctionvec}
\alias{Rcpp_cplfunctionvec-class}
%\alias{show-methods}
%\alias{show,ANY-method}
\alias{SerialOptimPriceStorage}
\alias{show,Rcpp_cplfunctionvec-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
This class implements "optimized list" of continuous convex piecewise linear functions
}
\description{
This is a wrapper to stl vector of convex piecewise linear functions. Allows to loop efficiently on such list. 
}

\author{
Robin Girard
}

\seealso{
 to See Also as \code{\link{cplfunction}}, \code{\link{cpqfunctionvec}}
}
\examples{



####
# construction of a vector of 
# continuous convex piecewise linear functions

CCPWLfuncList=new(cplfunctionvec) 
CCPWLfuncList$push_back(new(cplfunction,c(-1,1) ,c(-Inf,0),0))
CCPWLfuncList$push_back(new(cplfunction,c(-1,1) ,c(-Inf,0),0))

CCPWLfuncList=new(cplfunctionvec) 
n=1000; Y=rnorm(n); S1=array(-1,n);S2=array(1,n); B0=array(-Inf,n); B1=rnorm(n); 
for (i in 1:n){
  CCPWLfuncList$push_back(new(cplfunction,c(S1[i],S2[i]) ,c(B0[i],B1[i]),0))
}
CCPWLfuncList$size() ## gives the size 
## The same but faster
CCPWLfuncList=new(cplfunctionvec) 
CCPWLfuncList$SerialPush_2Breaks_Functions(S1,S2,B0,B1);

#### method OptimMargInt solves 
#         		min_x sum_i=1^n C_i(x_i)
#                   Pmoins_i<=	x_i				<=Pplus_i 		i=1,...,n
#					Cmoins_i<=	sum_j=1^i x_j	<=Cplus_i 	i=1,...,n

Pmoins=array(-1,n);Pplus=array(1,n);Cmoins=array(0,n);Cplus=array(5,n);
res=CCPWLfuncList$OptimMargInt(Pmoins,Pplus,Cmoins,Cplus)

par(mfrow=c(1,2))
plot(Y,type='l',ylim=range(res$xEtoile))
lines(y=Pmoins,x=1:n,col='blue'); lines(y=Pplus,x=1:n,col='blue');
lines(y=res$xEtoile,x=1:n,col='red')
text(x=800,y=3,paste("Optimum=",signif(sum(abs(res$xEtoile-Y)),digits=6)))
plot(Y,type='l',ylim=c(min(Y),max(diffinv(res$xEtoile)[1:n+1])))
lines(y=Cmoins,x=1:n,col='blue'); lines(y=Cplus,x=1:n,col='blue');
lines(y=diffinv(res$xEtoile)[1:n+1],x=1:n,col='red')

rm(list=ls())
gc()



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
