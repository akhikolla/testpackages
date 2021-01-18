\name{SAMM}
\alias{SAMM}
\title{Some Algoritms for Mixed Models}
\description{
The LMM model of focus can be expressed as \deqn{Y=XB+\sum_{j=1}^k Z_j G_j+WE,} where \eqn{Y} is the \eqn{n \times d} response variable, \eqn{X} is the \eqn{n \times q} design matrix of \eqn{q \times d} the fixed effects \eqn{B}, \eqn{Z_j} for \eqn{j=1,2,\ldots,k}  (\eqn{k\geq 1}) are the \eqn{n \times q_j} design matrices of the \eqn{q_j \times d} random effects \eqn{G_j,}  and \eqn{W} is the \eqn{n \times t} design matrix of \eqn{t \times d} residual effecs \eqn{E}. The random effects and the residual are independently distributed, and have matrix variate distributions (\eqn{G_j\sim N_{q_j \times d}(0_{q_j \times d}, K_j,\Sigma_j)} for \eqn{j=1,2,\ldots,k} and \eqn{E\sim N_{t \times d}(0_{t \times d}, R,\Sigma_E)}). The matrices  \eqn{K_j, R, \Sigma_j, \Sigma_E} might be further parametrized.  When the response is univariate \eqn{(d=1)} there is no need to specify a structure for \eqn{\Sigma_j}'s and\eqn{\Sigma_E.}  On the other hand, the \eqn{K_j} and \eqn{R} are kernel matrices (covariance matrices with some standardization over the diagonals) and they need to be provided by the user. 

Please refer to the examples below and the other help files for more details about the kernel and covariance structures. 

}
\usage{
SAMM(Y, X, Zlist, Klist, lambda, W, R, Siglist, corfunc,
                 corfuncfixed, sigfunc, mmalg, tolparconv = 1e-10,
                 tolparinv = 1e-10, maxiter = 1000L, geterrors = FALSE,
                 Hinv = FALSE)
}
\arguments{
  \item{Y}{\eqn{Y} is the \eqn{n \times d} response variable}
  \item{X}{\eqn{n \times q} design matrix of \eqn{q \times d} the fixed effects \eqn{B}}
  \item{Zlist}{a list object containing \eqn{Z_j} for \eqn{j=1,2,\ldots,k}  (\eqn{k\geq 1}), the \eqn{n \times q_j} design matrices of the \eqn{q_j \times d} random effects \eqn{G_j,}}
  \item{Klist}{a list to specify the kernel matrices \eqn{K_j} for \eqn{j=1,2,\ldots,k.} For each \eqn{j}, the user needs to provide a constant matrix, or a list specifying the kernel structure}
  \item{lambda}{a scalar shrinkage parameter for shrinkage of variance components (only works with the choice mmalg=''mmmk_ml'').}
  \item{W}{\eqn{W} is the \eqn{n \times t} design matrix of \eqn{t \times d} residual effecs \eqn{E}}
  \item{R}{a list to specify the kernel matrix \eqn{R},the user needs to provide a constant matrix, or a list specifying the kernel structure}
  \item{Siglist}{a list to specify the covariance structures \eqn{\Sigma_j}'s and\eqn{\Sigma_E}}
  \item{corfunc}{a boolian vector specifying whether \eqn{K_j} for \eqn{j=1,2,\ldots,k}  and \eqn{R} are functions or given matrices (TRUE for functions)}
  \item{corfuncfixed}{a boolian vector specifying whether \eqn{K_j} for \eqn{j=1,2,\ldots,k} and \eqn{R} are fixed at the initial parameter values specified }
  \item{sigfunc}{a boolian vector specifying whether \eqn{\Sigma_j} for \eqn{j=1,2,\ldots,k, E} are functions or unstructured (TRUE for functions)}
  \item{mmalg}{The mixed model solving algorithm}
  \item{tolparconv}{convergence criteria}
  \item{tolparinv}{a small scalar to add to the diagonals of positive semidefinite matrices for inversion or for calculating the Cholesky decompositions}
  \item{maxiter}{Maximum number of iterations}
  \item{geterrors}{TRUE or FALSE, if true prediction error variances for the random effects are supplied in the output}
\item{Hinv}{TRUE or FALSE, if TRUE the inverse of the}
}
\value{This might change with respect to the algorithm or analysis.
}
\details{
This might change with respect to the algorithm or analysis.
}

\examples{
\dontrun{
library(SAMM)
#set.seed(12345)

n=120
ntrain=100
M1<-matrix(rbinom(n*180,2,.3)-1, nrow=n)
K<-relmatcov_cppforR(c(0), M1)
K[1:5,1:5]
det(K)
K=K+1e-3*diag(n)
mean(diag(K))

covY<-2*K+1*diag(n)

Y<-10+crossprod(chol(covY),rnorm(n))


#training set
Trainset<-sample(1:n,ntrain,replace=(ntrain>n))
y=Y[Trainset]
X=matrix(rep(1, n)[Trainset], ncol=1)
Z=diag(n)[Trainset,]

X=X
y=y



########1
samout1<-SAMM(Y=matrix(y,ncol=1),X=X, Zlist=list(Z), Klist=list(K),
lambda=0, W=diag(ntrain),R=list(diag(ntrain)),Siglist=list(),
corfunc=c(F,F), corfuncfixed=c(F,F), sigfunc=c(F),mmalg="emm_ml",
tolparconv=1e-10, tolparinv=1e-10,maxiter=1000,geterrors=F)

samout2<-SAMM(Y=matrix(y,ncol=1),X=X, Zlist=list(Z), Klist=list(K),
lambda=0, W=diag(ntrain),R=list(diag(ntrain)),Siglist=list(),
corfunc=c(F,F), corfuncfixed=c(F,F), sigfunc=c(F),mmalg="emm_reml",
tolparconv=1e-10, tolparinv=1e-10,maxiter=1000,geterrors=F)

samout3<-SAMM(Y=matrix(y,ncol=1),X=X, Zlist=list(Z),
Klist=list(list("Const",c(0),K)), lambda=0, 
W=diag(ntrain),R=list(diag(ntrain)),Siglist=list(),
corfunc=c(F,F), corfuncfixed=c(F,F), sigfunc=c(F),
mmalg="dmm_ml", tolparconv=1e-10, tolparinv=1e-10,
maxiter=1000,geterrors=F)

samout4<-SAMM(Y=matrix(y,ncol=1),X=X, Zlist=list(Z), 
Klist=list(list("Const",c(0),K)), lambda=0, 
W=diag(ntrain),R=list(diag(ntrain)),Siglist=list(), 
corfunc=c(F,F), corfuncfixed=c(F,F), sigfunc=c(F),
mmalg="dmm_reml", tolparconv=1e-10, tolparinv=1e-10,
maxiter=1000,geterrors=F)

samout5<-SAMM(Y=matrix(y,ncol=1),X=X, Zlist=list(Z),
Klist=list(K), lambda=0, W=diag(ntrain),
R=list(diag(ntrain)),Siglist=list(), corfunc=c(F,F), 
corfuncfixed=c(F,F), sigfunc=c(F),mmalg="mm_ml", 
tolparconv=1e-10, tolparinv=1e-10,maxiter=1000,geterrors=F)

samout6<-SAMM(Y=matrix(y,ncol=1),X=X, Zlist=list(Z), 
Klist=list(K), lambda=0, W=diag(ntrain),R=list(diag(ntrain)),
Siglist=list(), corfunc=c(F,F), corfuncfixed=c(F,F), 
sigfunc=c(F),mmalg="dermm_reml1", tolparconv=1e-10, 
tolparinv=1e-10,maxiter=1000,geterrors=F)

samout7<-SAMM(Y=matrix(y,ncol=1),X=X, Zlist=list(Z),
Klist=list(K), lambda=0, W=diag(ntrain),R=list(diag(ntrain)),
Siglist=list(), corfunc=c(F,F), corfuncfixed=c(F,F), 
sigfunc=c(F),mmalg="dermm_reml2", tolparconv=1e-10, 
tolparinv=1e-10,maxiter=1000,geterrors=F)

samout8<-SAMM(Y=matrix(y,ncol=1),X=X, Zlist=list(Z),
Klist=list(K), lambda=0, W=diag(ntrain),R=list(diag(ntrain)),
Siglist=list(), corfunc=c(F,F), corfuncfixed=c(F,F), 
sigfunc=c(F),mmalg="mmmk_ml", tolparconv=1e-10, 
tolparinv=1e-10,maxiter=1000,geterrors=F)

samout9<-SAMM(Y=matrix(y,ncol=1),X=X, Zlist=list(Z), 
Klist=list(K), lambda=0, W=diag(ntrain),R=list(diag(ntrain)),
Siglist=list(), corfunc=c(F,F), corfuncfixed=c(F,F), 
sigfunc=c(F),mmalg="emmmk_reml", tolparconv=1e-10, 
tolparinv=1e-10,maxiter=1000,geterrors=F)

samout10<-SAMM(Y=matrix(y,ncol=1),X=X, Zlist=list(Z), 
Klist=list(K), lambda=0, W=diag(ntrain),R=list(diag(ntrain)),
Siglist=list(), corfunc=c(F,F), corfuncfixed=c(F,F), 
sigfunc=c(F),mmalg="emmmk_ml", tolparconv=1e-10, 
tolparinv=1e-10,maxiter=1000,geterrors=F)

samout11<-SAMM(Y=matrix(y,ncol=1),X=X, Zlist=list(Z),
Klist=list(K), lambda=0, W=diag(ntrain),R=list(diag(ntrain)),
Siglist=list(), corfunc=c(F,F), corfuncfixed=c(F,F), sigfunc=c(F),
mmalg="emmmv_ml", tolparconv=1e-10, tolparinv=1e-10,maxiter=1000,geterrors=F)

###########2
R<-diag(ntrain)
diag(R)<-diag(R)+.01*rnorm(nrow(R))

samout12<-SAMM(Y=matrix(y,ncol=1),X=X, Zlist=list(Z),
Klist=list(K), lambda=0, W=diag(ntrain),R=list(R),
Siglist=list(), corfunc=c(F,F), corfuncfixed=c(F,F),
sigfunc=c(F),mmalg="mmmk_ml", tolparconv=1e-10, 
tolparinv=1e-10,maxiter=1000,geterrors=F)

samout13<-SAMM(Y=matrix(y,ncol=1),X=X, Zlist=list(Z), 
Klist=list(K), lambda=0, W=diag(ntrain),R=list(R),
Siglist=list(), corfunc=c(F,F), corfuncfixed=c(F,F), 
sigfunc=c(F),mmalg="mm_ml", tolparconv=1e-10, 
tolparinv=1e-10,maxiter=1000,geterrors=F)

samout14<-SAMM(Y=matrix(y,ncol=1),X=X, Zlist=list(Z),
Klist=list(K), lambda=0, W=diag(ntrain),R=list(R),
Siglist=list(), corfunc=c(F,F), corfuncfixed=c(F,F), 
sigfunc=c(F),mmalg="mmmv_ml", tolparconv=1e-10, 
tolparinv=1e-10,maxiter=1000,geterrors=F)

samout15<-SAMM(Y=matrix(y,ncol=1),X=X, Zlist=list(Z), 
Klist=list(K), lambda=0, W=diag(ntrain),R=list(R),
Siglist=list(), corfunc=c(F,F), corfuncfixed=c(F,F), 
sigfunc=c(F),mmalg="mmmkmv_ml", tolparconv=1e-10, 
tolparinv=1e-10,maxiter=1000,geterrors=F)

samout16<-SAMM(Y=matrix(y,ncol=1),X=X, Zlist=list(Z), 
Klist=list(K), lambda=0, W=diag(ntrain),R=list(R),Siglist=list(),
corfunc=c(F,F), corfuncfixed=c(F,F), sigfunc=c(F),
mmalg="dermm_reml1", tolparconv=1e-10,
tolparinv=1e-10,maxiter=1000,geterrors=F)

samout17<-SAMM(Y=matrix(y,ncol=1),X=X, Zlist=list(Z), 
Klist=list(K), lambda=0, W=diag(ntrain),R=list(R),
Siglist=list(), corfunc=c(F,F), corfuncfixed=c(F,F),
sigfunc=c(F),mmalg="dermm_reml2", tolparconv=1e-10, 
tolparinv=1e-10,maxiter=1000,geterrors=F)

samout18<-SAMM(Y=matrix(y,ncol=1),X=X, Zlist=list(Z),
Klist=list(list("Const",c(0),K)), lambda=0, W=diag(ntrain),
R=list(R),Siglist=list(), corfunc=c(F,F), 
corfuncfixed=c(F,F), sigfunc=c(F),mmalg="dmm_ml",
tolparconv=1e-10, tolparinv=1e-10,maxiter=1000,geterrors=F)

samout19<-SAMM(Y=matrix(y,ncol=1),X=X, Zlist=list(Z), 
Klist=list(list("Const",c(0),K)), lambda=0, W=diag(ntrain),
R=list(R),Siglist=list(), corfunc=c(F,F), corfuncfixed=c(F,F), 
sigfunc=c(F),mmalg="dmm_reml", tolparconv=1e-10,
tolparinv=1e-10,maxiter=1000,geterrors=F)




##3
samout20<-SAMM(Y=matrix(y,ncol=1),X=X, Zlist=list(Z), 
Klist=list(list("Const",c(0),K)), lambda=0, W=diag(ntrain),
R=list(R),Siglist=list(), corfunc=c(T,F), corfuncfixed=c(F,F),
sigfunc=c(F),mmalg="dmm_reml", tolparconv=1e-10, 
tolparinv=1e-10,maxiter=1000,geterrors=F)

samout21<-SAMM(Y=matrix(y,ncol=1),X=X, Zlist=list(Z),
Klist=list(list("ar1",c(0),matrix(n,1,1))), lambda=0, 
W=diag(ntrain),R=list(R),Siglist=list(), corfunc=c(T,F),
corfuncfixed=c(F,F), sigfunc=c(F),mmalg="dmm_reml",
tolparconv=1e-10, tolparinv=1e-10,maxiter=1000,geterrors=F)


samout22<-SAMM(Y=matrix(y,ncol=1),X=X, Zlist=list(Z),
Klist=list(list("ar1het",rep(0,n),matrix(n,1,1))), lambda=0,
W=diag(ntrain),R=list(R),Siglist=list(), corfunc=c(T,F),
corfuncfixed=c(F,F), sigfunc=c(F),mmalg="dmm_reml",
tolparconv=1e-10, tolparinv=1e-10,maxiter=1000,geterrors=F)

samout1$Vu
samout2$Vu
samout3$Vu
samout4$Vu
samout5$Vu
samout6$Vu
samout7$Vu
samout8$Vu
samout9$Vu
samout10$Vu
samout11$Vu

samout12$Vu
samout13$Vu
samout14$sigmahatlist[[1]]
samout15$sigmahatlist[[1]]
samout16$Vu
samout17$Vu
samout18$Vu
samout19$Vu

samout20$Ve
str(samout19)


str(samout20)
str(samout21)
str(samout22)




###

samout1$Vu
samout2$Vu
samout3$Vu
samout4$Vu
samout5$Vu
samout6$Vu
samout7$Vu
samout8$Vu
samout9$Vu
samout10$Vu
samout11$Vgt

samout12$Vu
samout13$Vu
samout14$sigmahatlist[[1]]
samout15$sigmahatlist[[1]]
samout16$Vu
samout17$Vu
samout18$Vu
samout19$Vu

samout20$Ve
str(samout19)


str(samout20)
str(samout21)
str(samout22)



###############################


n=100
nsample=80
rhotrans=5
rho=tan((pi/2)*print((2/pi)*atan(rhotrans)))
rho
M1<-matrix(rbinom(n*300, 2, .2)-1, nrow=n)
K1<-relmatcov_cppforR(c(.01), M1)

M2<-matrix(rbinom(n*300, 2, .2)-1, nrow=n)
K2<-relmatcov_cppforR(c(-0.3), M2)
##K2<-ar1_cppforR(c(20), matrix(n))
W=(diag(5)[sample(1:5,n, replace=T),])
covY<-3*K1+5*K2+10*(W\%*\%ar1cov_cppforR(c(rhotrans),matrix(5))\%*\%t(W))
K1[1:5,1:5]
dim(W)
dim(ar1cov_cppforR(c(6),matrix(5)))
Y<-10+crossprod(chol(covY),rnorm(n))


#training set
Trainset<-sample(1:n,nsample)
ytrain=Y[Trainset]
Xtrain=matrix(rep(1, n)[Trainset], ncol=1)
Ztrain=diag(n)[Trainset,]
Wtrain=W[Trainset,]



samout27<-SAMM(Y=matrix(ytrain,ncol=1),X=Xtrain,
Zlist=list(Ztrain, Ztrain), Klist=list(K1,K2), 
lambda=0, W=Wtrain,R=list(list("ar1",c(0),matrix(5,1,1))),
Siglist=list("","",""), corfunc=c(F,F,T), 
corfuncfixed=c(F,F,F), sigfunc=c(F,F,F),mmalg="mmmk_ml",
tolparconv=1e-10, tolparinv=1e-10,maxiter=1000,geterrors=F)
str(samout27)


samout29<-SAMM(Y=matrix(ytrain,ncol=1),X=Xtrain,
Zlist=list(Ztrain, Ztrain), Klist=list(K1,K2),
lambda=0, W=Wtrain,R=list("ar1",c(0),matrix(5,1,1)),
Siglist=list("","",""), corfunc=c(F,F,T),
corfuncfixed=c(F,F,F), sigfunc=c(F,F,F),
mmalg="dermm_reml1", tolparconv=1e-10,
tolparinv=1e-10,maxiter=1000,geterrors=F)

str(samout29)

samout30<-SAMM(Y=matrix(ytrain,ncol=1),X=Xtrain,
Zlist=list(Ztrain, Ztrain), Klist=list(K1,K2),
lambda=0, W=Wtrain,R=list("ar1",c(0),matrix(5,1,1)),
Siglist=list("","",""), corfunc=c(F,F,T), corfuncfixed=c(F,F,F),
sigfunc=c(F,F,F),mmalg="dermm_reml2", tolparconv=1e-10,
tolparinv=1e-10,maxiter=1000,geterrors=F)
str(samout30)


###########################

n=100
M1<-matrix(rbinom(n*150, 2, .2)-1, nrow=n)
K1<-relmatcov_cppforR(c(0), M1)
M2<-matrix(rbinom(n*150, 2, .2)-1, nrow=n)
K2<-relmatcov_cppforR(c(0), M2)
M3<-matrix(rbinom(n*100, 2, .2)-1, nrow=n)
K3<-relmatcov_cppforR(c(0), M3)
M4<-matrix(rbinom(n*100, 2, .2)-1, nrow=n)
K4<-relmatcov_cppforR(c(0), M4)


covY<-2*K1+3*K2+1*K3+2*K4+3*diag(n)

Y<-10+crossprod(chol(covY),rnorm(n))


#training set
Trainset<-sample(1:n,80)
y=Y[Trainset]
X=matrix(rep(1, n)[Trainset], ncol=1)
Z=diag(n)[Trainset,]

X=X
y=y

samout35<-SAMM(Y=matrix(y,ncol=1),X=X,Zlist=list(Z,Z,Z,Z),
lambda=0,Klist=list(K1,K2,K3,K4), W=Z,R=list(diag(n)),
Siglist=list("","","","",""), corfunc=c(F,F,F,F,F), 
corfuncfixed=c(T,T,T,T,T),sigfunc=c(F,F,F,F,F),
mmalg="mmmk_ml", tolparconv=1e-10,
tolparinv=1e-10,maxiter=1000,geterrors=F)
str(samout35)

samout36<-SAMM(Y=matrix(y,ncol=1),X=X,Zlist=list(Z,Z,Z,Z),
lambda=0.99999,Klist=list(K1,K2,K3,K4), W=Z,
R=list(diag(n)),Siglist=list("","","","",""),
corfunc=c(F,F,F,F,F), corfuncfixed=c(T,T,T,T,T),
sigfunc=c(F,F,F,F,F),mmalg="mmmk_ml", tolparconv=1e-10,
tolparinv=1e-10,maxiter=1000,geterrors=F)
str(samout36)


outmat<-c()
for (lambda in seq(0,.999999, length=30)){
  
  
  samout37<-SAMM(Y=matrix(y,ncol=1),X=X,
  Zlist=list(Z,Z,Z,Z),lambda=lambda,Klist=list(K1,K2,K3,K4),
  W=Z,R=list(diag(n)),Siglist=list("","","","",""),
  corfunc=c(F,F,F,F,F), corfuncfixed=c(T,T,T,T,T), 
  sigfunc=c(F,F,F,F,F),mmalg="mmmk_ml",
  tolparconv=1e-10, tolparinv=1e-10,
  maxiter=1000,geterrors=F)
  
  
  outmat<-cbind(outmat,c(samout37$Vu*samout37$weights, samout37$Ve))
  
}
str(samout37)
colnames(outmat)<-seq(0,.999999, length=30)

maxmat<-max(c(outmat))
minmat<-min(c(outmat))

plot(seq(0,.999999, length=30),outmat[1,],
ylim=c(minmat-1, maxmat+1), col=2, type="l")
for (i in 2:5){
  par(new=T)
  plot(seq(0,.999999, length=30), outmat[i,],
  axes=F, ylim=c(minmat-1, maxmat+1), col=i+1, type="l", xlab="", ylab="")
}

} 
}
\author{Deniz Akdemir}

