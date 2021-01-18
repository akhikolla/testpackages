\name{sigcovfuncs_cppforR}

\alias{FA1hetSig_cppforR}
\alias{FA1homSig_cppforR}
\alias{FAhetSig_cppforR}
\alias{FAhomSig_cppforR}
\alias{UnstrKronIdentSig_cppforR}
\alias{IdentKronUnstrSig_cppforR}
\alias{IdentSig_cppforR}
\alias{ar1hetcov_cppforR}
\alias{arma11cov_cppforR}
\alias{compsymmcov_cppforR}
\alias{compsymmhetSig_cppforR}
\alias{compsymmhetcov_cppforR}
\alias{compsymmhomSig_cppforR}
\alias{diagSig_cppforR}
\alias{diagcov_cppforR}
\alias{expcov_cppforR}
\alias{lincombcov_cppforR}
\alias{rbfcov_cppforR}
\alias{ar1cov_cppforR}
\alias{unstrKronKcov_cppforR}
\alias{relmatcov_cppforR}
\alias{UnstrKronUnstrcov_cppforR}
\alias{unstrcov_cppforR}
\alias{ConstMatcov_cppforR}
\alias{KKronar1cov_cppforR}
\alias{KKronar1hetcov_cppforR}
\alias{KKronarma11cov_cppforR}
\alias{KKroncompsymmcov_cppforR}
\alias{KKroncompsymmhetcov_cppforR}
\alias{KKrondiagcov_cppforR}
\alias{KKronunstrcov_cppforR}
\alias{ar1KronKcov_cppforR}
\alias{ar1hetKronKcov_cppforR}
\alias{arma11KronKcov_cppforR}
\alias{compsymmKronKcov_cppforR}
\alias{compsymmhetKronKcov_cppforR}
\alias{diagKronKcov_cppforR}
\alias{expdistcov_cppforR}
\alias{rbfdistcov_cppforR}
\alias{splincov_cppforR}
\alias{splinlogcov_cppforR}
\alias{sppowcov_cppforR}

\title{Covariance and Sigma Functions}
\description{
The "kernel" functions end with "cov_cppforR" and sigma functions end with "Sig_cppforR". Check table below and the examples for details and usage. Documentation for some of these functions is missing. Let \eqn{K} be a given covariance matrix, \eqn{D} be a given Euclidean distance matrix, let \eqn{q} be the dimension of desired sigma or covariance matrices. 
}

\arguments{
  \item{params}{a numeric vector for the parameters (a mapping of the original parameters)}
  \item{data}{a numeric matrix (see examples and details)}
  
}
\value{
A kernel or sigma matrix. 
}
\details{
\tabular{rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr}{
  \strong{function name}\tab \strong{function ref} \tab \strong{params} \tab \strong{data} \tab \strong{formula} \cr
 

   IdentSig_cppforR \tab Ident \tab 1 \tab  matrix(q,1,1) \tab \eqn{A=\sigma I_q}\cr
   \tab   \tab   \tab   \tab  \cr

  ar1hetcov_cppforR \tab ar1het \tab q \tab matrix(q,1,1) \tab \eqn{A_{ij}=\sigma_i\sigma_j\rho^{|i-j|}},\cr
   \tab   \tab   \tab   \tab \eqn{\sigma_1=1} \cr

  arma11cov_cppforR\tab arma11 \tab 2 \tab matrix(q,1,1) \tab \eqn{A_{ij}=(\lambda\rho^{|i-j|-1}*}\cr
   \tab   \tab   \tab   \tab \eqn{1(i\neq j)+1(i=j))} \cr
 FA1hetSig_cppforR \tab FA1het \tab 2q-2 \tab matrix(q,1,1) \tab \eqn{f_if_j+\sigma_{ii}1(i=j) }\cr
   \tab   \tab   \tab   \tab \eqn{f_1=1,} \eqn{\sigma_1=1} \cr

  FA1homSig_cppforR \tab FA1hom \tab q \tab  matrix(q,1,1) \tab \eqn{f_if_j+\sigma 1(i=j)} \cr
    \tab   \tab   \tab   \tab \eqn{f_1=1} \cr
    
  compsymmcov_cppforR \tab compsymm \tab 1 \tab matrix(q,1,1) \tab  \eqn{A_{ij}=\rho}, \eqn{i\neq j;} \cr
   \tab   \tab   \tab   \tab \eqn{A_{ij}=1}, \eqn{i=j} \cr

   compsymmhetSig_cppforR \tab compsymmhet \tab q+1 \tab matrix(q,1,1) \tab \eqn{A_{ij}=\sigma_i\sigma_j\rho,} \eqn{i\neq j}\cr
   \tab   \tab   \tab   \tab \eqn{A_{ij}=\sigma_i\sigma_j,} \eqn{i= j} \cr

  compsymmhetcov_cppforR \tab compsymmhet \tab q \tab matrix(q,1,1) \tab \eqn{A_{ij}=\sigma_i\sigma_j\rho,} \eqn{i\neq j}\cr
   \tab   \tab   \tab   \tab \eqn{A_{ij}=\sigma_i\sigma_j,} \eqn{i= j} \cr
   \tab   \tab   \tab   \tab \eqn{\sigma_1=1}  \cr

  compsymmhomSig_cppforR \tab compsymmhom \tab 2 \tab matrix(q,1,1) \tab \eqn{A_{ij}=\sigma\rho}, \eqn{i\neq j;} \cr
   \tab   \tab   \tab   \tab \eqn{A_{ij}=\sigma}, \eqn{i=j} \cr

  diagSig_cppforR \tab diag \tab q \tab  matrix(q,1,1) \tab \eqn{diag(\sigma_1,..,\sigma_q)}\cr
   \tab   \tab   \tab   \tab  \cr

   diagcov_cppforR \tab Diag \tab q-1 \tab matrix(q,1,1) \tab \eqn{diag(1,\sigma_2,..,\sigma_q)} \cr

  expcov_cppforR \tab exp \tab 1 \tab M \tab \eqn{exp(-\sigma*d_{ij})}\cr
   \tab   \tab   \tab   \tab M defines D \cr

  lincombcov_cppforR \tab lincomb \tab k \tab \eqn{(K_1;\ldots; K_k)} \tab \eqn{\sum^k_{j=1}w_jK_j} \cr

   rbfcov_cppforR \tab rbf \tab 1 \tab M \tab \eqn{exp(-\sigma*d^2_{ij})}\cr
   \tab   \tab   \tab   \tab M defines D \cr

  ar1cov_cppforR \tab ar1 \tab 1 \tab matrix(nrow, 1,1) \tab \eqn{A_{ij}=\rho^{|i-j|}}\cr

  relmatcov_cppforR \tab RelMat\tab 1 \tab M \tab Genetic Similarity+\cr
   \tab   \tab   \tab M is coded as -1, 0, 1  \tab  \eqn{\sigma I}\cr

  unstrcov_cppforR \tab unstr \tab \eqn{\frac{q(q+1)}{2}-1} \tab matrix(q,1,1) \tab \eqn{A_{ij}=\sigma_i\sigma_j\rho_{ij}}\cr
   \tab   \tab   \tab   \tab\eqn{\sigma_1=1,} \cr
   \tab   \tab   \tab   \eqn{\rho_{ii}=1} \cr

  ConstMatcov_cppforR \tab Const \tab 0 \tab K \tab \eqn{A=K}\cr
  
  expdistcov_cppforR \tab expdist \tab 1 \tab D \tab \eqn{exp(-\sigma d_{ij})}\cr
   
   rbfdistcov_cppforR \tab rbfdist \tab 1 \tab D \tab \eqn{exp(-\sigma d^2_{ij})}\cr
   
   splincov_cppforR \tab splin \tab 1 \tab D \tab \eqn{1-\rho d_{ij},} \eqn{\rho d_{ij}\leq 1} \cr
   \tab   \tab   \tab   \tab \eqn{0,} \eqn{\rho d_{ij}> 1} \cr
   
   sppowcov_cppforR \tab sppow \tab 1 \tab D \tab \eqn{\rho^{d_{ij}}}

 }
}

\examples{
\dontrun{
library(SAMM)
n=100
nsample=80
rhotrans=5
ar1cov_cppforR(c(rhotrans),matrix(5))
rho=(2/pi)*atan(rhotrans)
rho
tan((pi/2)*(rho))

M1<-matrix(rbinom(n*300, 2, .2)-1, nrow=n)
K1<-relmatcov_cppforR(c(.01), M1)

M2<-matrix(rbinom(n*300, 2, .2)-1, nrow=n)
K2<-relmatcov_cppforR(c(0.03), M2)
W=(diag(5)[sample(1:5,n, replace=TRUE),])
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

samout<-SAMM(Y=matrix(ytrain,ncol=1),X=Xtrain,
Zlist=list(Ztrain, Ztrain), Klist=list(K1,K2),
lambda=0, W=Wtrain,R=list("ar1",c(0),matrix(5,1,1)),
Siglist=list("","",""), corfunc=c(F,F,T), corfuncfixed=c(F,F,F),
sigfunc=c(F,F,F),mmalg="dermm_reml2", tolparconv=1e-10,
tolparinv=1e-10,maxiter=1000,geterrors=F)
samout$corfuncparamslist[[3]]
rhohat=(2/pi)*atan(samout$corfuncparamslist[[3]])
rhohat
ar1cov_cppforR(c(samout$corfuncparamslist[[3]]),matrix(5,1,1))
}
}

\author{Deniz Akdemir}

