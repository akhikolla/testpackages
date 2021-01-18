#
#    Copyright (C) 2020 David Preinerstorfer
#    david.preinerstorfer@ulb.ac.be
#
#    This file is a part of wbsd.
#
#    wbsd is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details. A copy may be obtained at
#    http://www.r-project.org/Licenses/

###########################################################################
###########################################################################
#
# Wrapper function for test statistic of the form (in case y is vector):
#
# ``(R\hat{\beta}(y) - rboot)' VCESTIMATOR^{-1} (R\hat{\beta}(y) - rboot)''
#
# here \hat{\beta} denotes the OLS estimator
# and VCESTIMATOR is a covariance matrix estimator that might be based
# on restricted or unrestricted residuals (restriction: R\beta = r);
# (the choice of restricted or unrestricted residuals also has an effect
# on the weights in the construction of the HC1 - HC4 estimators)
#
###########################################################################
###########################################################################

F.wrap <- function(
y,          #matrix (n rows) of observations 
R,          #restriction matrix (q times k, rank q)
r,          #right-hand side in restriction (q-vector)
rboot,      #centralisation used in test statistic (typically equal to r;
            #but in case of bootstrap with unrestricted centering required 
            #to be ``R \hat{\beta} (y)''); (q-vector)
X,          #design matrix (n times k, rank k)
n,          #sample size
k,          #number of columns of X
q,          #number of rows of R
qrX,        #qr decomposition of X
hcmethod,   #-1:4; -1 = classical F-test without df adjustment; 0 = HC0, 1 = HC1, etc.
cores,      #number of cores used in computation
restr.cov,  #Covariance Matrix estimator computed with restricted residuals (TRUE or FALSE)  
Bfac,       #R(X'X)^{-1}X' (cf. function Bfactor.matrix below)
Bfac2,      #(R(X'X)^(-1)R')^(-1) (cf. function Bfactor.matrix2 below)
qrM0lin,    #qr decomposition of M0lin (cf. function M0lin below)
RF,         #RF = (X'X)^(-1)R'(R(X'X)^(-1)R')^(-1) used in res.OLSRF (cf. below)
tol         #tolerance parameter used in checking invertibility of VCESTIMATOR in test statistic
            #nonpositive numbers are reset to tol = machine epsilon (.Machine$double.eps)
            #if VCESTIMATOR is not invertible, -1 is returned 
            #if VCESTIMATOR is invertible, the test statistic is always nonnegative
){

#computation of residuals and weights used in the construction of the covariance matrix estimators

if(restr.cov == TRUE){
 umat <- res.OLSRF(y, qrX, R, r, RF)$Res.res
 Wmat <- wvec(k-q, n, qrM0lin, hcmethod)
}

if(restr.cov ==  FALSE){
 umat <- qr.resid(qrX, y)
 Wmat <- wvec(k, n, qrX, hcmethod)
}

Rbmat <- R %*% qr.coef(qrX, y) - matrix(rboot, byrow = FALSE, nrow = q, ncol = dim(y)[2])

if(hcmethod == -1){
 if(tol <= 0){
 tol <- .Machine$double.eps
 }
 df <- n-(k - restr.cov*q)
 var.est <- apply(umat^2,2,sum)/df
 nonzero <- (var.est > tol)
 test.val <- nonzero*1
 CBFac <- chol(Bfac2)
 Wmat <- CBFac%*%Rbmat
 Wmat <- apply(Wmat^2, 2, sum)
 test.val[test.val == 1] <- (Wmat[test.val == 1]/var.est[test.val == 1])
 test.val[test.val == 0] <- -1
} else {
test.val <- .Call('_wbsd_ctest', PACKAGE = 'wbsd', umat, Rbmat, Wmat, Bfac, cores, tol)
}

return(test.val)
}

###########################################################################
###########################################################################
#
# Bootstrap Sample Generation
#
# Given observation y, generate sample 
#
#             Xbeta.est(y) + diag(res(y)) xi
#
# Here res(y) can be the restricted or unrestricted OLS residuals,
# depending on whether boot.res.restr is TRUE or FALSE, respectively,
# and beta.est can be a null-restricted or unrestricted OLS estimator
# depending on whether boot.center.restr is TRUE or FALSE, respectively.
#
###########################################################################
###########################################################################

boot.sample <- function(
#INPUT
y,                  #matrix (n times 1) of observations
xi,                 #bootstrap xis (matrix n rows, each column is a bootstrap obs)
qrX,                #qr decomposition of X (n times k, rank k)
R,                  #restriction matrix (q times k, rank q)
r,                  #right-hand side in restriction to be tested (q-vector)
boot.res.restr,     #bootstrap variances using restr. residuals (FALSE or TRUE)
boot.center.restr,  #bootstrap centers using restr. residuals (FALSE or TRUE)
RF = NULL           #RF = (X'X)^(-1)R'(R(X'X)^(-1)R')^(-1) used in res.OLSRF (cf. below)
                    #RF is only used in case boot.res.restr or boot.center.restr is TRUE
)
{

# compute restricted OLS estimator and residuals if needed

if(boot.res.restr == TRUE | boot.center.restr == TRUE){
  rer <- res.OLSRF(y, qrX, R, r, RF)
}

# obtain residuals from original sample (unrestricted or restricted)

if(boot.res.restr == FALSE){
  boot.res <- qr.resid(qrX, y)
} else {
  boot.res <- rer$Res.res
}

# compute the center of the bootstrap samples (unrestricted or restricted)

if(boot.center.restr == FALSE){
  boot.cen <- qr.fitted(qrX, y)
} else {
  boot.cen <- qr.X(qrX)%*%rer$Coef.restr
}

boot.cen <- matrix(rep(boot.cen, times = dim(xi)[2]), 
            nrow = dim(xi)[1], ncol = dim(xi)[2], byrow = FALSE)

# combine and return bootstrap sample

return(boot.cen + diag(c(boot.res))%*%xi)
}

###########################################################################
###########################################################################
#
# Check if a given design matrix X satisfies Assumption 1 or 2
# whether Assumption 1 or 2 is checked depends on the input qr decomposition
#
# if the qr decomposition of X is provided, Assumption 1 is checked
# if the qr decomposition of a basis of M0lin is provided, Assumption 2 is
# checked
# 
# Bfac = R(X'X)^{-1} X' (cf. the function Bfactor.matrix below)
#
###########################################################################
###########################################################################

As.check <- function(
qrXM,   #qr decomposition X (n times k, rank k), or of a basis of 
        #M0lin (n times (k-q), rank (k-q), cf. function M0lin below)
        #in the function we abbreviate ``k-q'' by l
Bfac,   #R(X'X)^{-1}X' (cf. function Bfactor.matrix below)
q,      #number of rows of restriction matrix R
as.tol  #tolerance parameter used in checking Assumptions 1 or 2, respectively
        #nonpositive numbers are reset to as.tol = machine epsilon (.Machine$double.eps)
){

XM <- qr.X(qrXM)
n <- dim(XM)[1]
l <- dim(XM)[2]

#if l == 0, the Assumption to be checked automatically satisfied

if(l == 0){
return(TRUE)
} else {

#if l != 0, do the following:

  e0 <- matrix(0, nrow = n, ncol = 1)  
  ind <- rep(NA, length = n) 
  
  #check whether elements of the canonical basis are in the span of XM

  for(i in 1:n){
    ei <- e0
    ei[i,1] <- 1
    ind[i] <- ( rrank(cbind(XM, ei), as.tol) > l )
  }
  
#if all entries of ind are TRUE, the Assumption is automatically
#satisfied, because of the rank assumptions on R and X
#if there are entries of ind that are FALSE, the Assumption
#has to be checked via a rank computation

  if(sum(ind) == n){
   return ( TRUE )
  } else {
   return( rrank(Bfac[, ind, drop=FALSE], as.tol) == q )
  }
}
}


###########################################################################
###########################################################################
#
# Compute the 
#
#   rank of a matrix based on the Eigen FullPivLU decomposition
#
# in case input is a matrix with one column or with one row, the 
# function checks if the maximal absolute entry of the input exceeds tol 
#
############################################################################
############################################################################

rrank <- function(
A,    #input matrix
tol   #tolerance parameter passed to Eigen FullPivLU decomposition;
      #if tolerance parameter is nonpositive, tol is reset to machine epsilon
){
if(min(dim(A)) == 1){
  if(tol <= 0){
  tol <- .Machine$double.eps
  }
 (max(abs(A)) > tol)*1
} else {
 .Call('_wbsd_rrank', PACKAGE = 'wbsd', A, tol)
}
}

###########################################################################
###########################################################################
#
# Compute the 
#
#     kernel of a matrix based on the Eigen FullPivLU decomposition
#
###########################################################################
###########################################################################

rkernel <- function(
A,  #input matrix
tol #tolerance parameter passed to Eigen FullPivLU decomposition;
    #if tolerance parameter is nonpositive, tol is reset to machine epsilon
){
.Call('_wbsd_rkernel', PACKAGE = 'wbsd', A, tol)
}

###########################################################################
###########################################################################
#
# Basis of $M_0^{lin} = \{y \in \mathbb{R}^n: y = Xb, Rb = 0\}$
#
#Computation is based on function rkernel with tol = -1 
#
# In case q = k, i.e., M_0^{lin} = (0), output vector of dimension n times 0
#
###########################################################################
###########################################################################

M0lin <- function(
X,  #input matrix (n times k, rank k)
R   #restriction matrix (q times k, rank q)
){
if(dim(R)[1] == dim(X)[2]){
  return(X[,-(1:dim(X)[2])]) 
} else {
  X%*%rkernel(R, -1)
}
}

###########################################################################
###########################################################################
#
# Compute the matrix
#
#     R(X'X)^(-1)X' 
#
###########################################################################
###########################################################################

Bfactor.matrix <- function(
qrX,        #qr decomposition of X (n times k, rank k)
n,          #number of rows of X
R = FALSE   #restriction matrix R (q times k, rank q)
            #if R = FALSE, then (X'X)^(-1)X' is returned
){
A <- qr.coef(qrX, diag(n))
if(is.matrix(R) == TRUE){
return(R%*%A)
} else {
return(A)
}
}

###########################################################################
###########################################################################
#
# Compute the matrix 
#
#     (R(X'X)^(-1)R')^(-1) 
#
###########################################################################
###########################################################################

Bfactor.matrix2 <- function(
qrX,      #qr decomposition of X (n times k, rank k)
R         #restriction matrix R (q times k, rank q)
){
factor.tmp <- qr.R(qrX)
factor.tmp <- R%*% backsolve(qr.R(qrX), diag(dim(factor.tmp)[1]))
Bfactor2 <- solve(tcrossprod(factor.tmp))
return(Bfactor2)
}

###########################################################################
###########################################################################
#
# Generate the 
#
#     weights d_i (unrestricted case) or \tilde{d}_i (restricted case)
#
# used in the construction of the covariance estimators H0 - H4 
#
# Whether d_i or \tilde{d}_i is generated is via qrXM, which determines the
# hat matrix used
#
###########################################################################
###########################################################################

wvec <- function(
l,        #see description of qrXM
n,        #see description of qrXM
qrXM,     #qr decomposition of matrix XM (n times l, rank l)
          #weights are based on the diag of hat matrix corresp. to XM
hcmethod  #0-4 (HC0, HC1, HC2, HC3, HC4)
){
if(hcmethod == 0 | l == 0){
 return(rep(1, length = n))
}

if(hcmethod == 1){
 return(rep(n/(n-l), length = n))
}

Q <- qr.Q(qrXM)
H <- diag(tcrossprod(Q,Q))
H[H == 1] <- 0
h <- 1/(1-H)

if(hcmethod == 2){
 return(h)
}

if(hcmethod == 3){
 return(h^2)
}

if(hcmethod == 4){
 delta <-  sapply(n*H/l, function(x) {min(4, x)})
 return(h^delta)
}

}

###########################################################################
###########################################################################
#
# Function that computes 
#
#     restricted OLS estimators and restricted residuals 
#
###########################################################################
###########################################################################

res.OLSRF <- function(
y,     #matrix (n rows) of observations 
qrX,   #qr decomposition of X (n times k, rank k)
R,     #restriction matrix (q times k, rank q)
r,     #right-hand side in the restriction (q-vector)
RF     #RF = (X'X)^(-1)R'(R(X'X)^(-1)R')^(-1)
){

OLS.coef <- qr.coef(qrX, y)
OLS.coefs.restr <- OLS.coef - RF%*%R%*%OLS.coef
OLS.coefs.restr <- OLS.coefs.restr + matrix(RF%*%r, nrow = dim(R)[2], ncol = dim(y)[2])
OLS.res.restr <- y - qr.X(qrX)%*%OLS.coefs.restr

return(list("Coef.restr" = OLS.coefs.restr, 
		"Res.res" = OLS.res.restr))
}

###########################################################################
###########################################################################
#
# Generator of support points of bootstrap distribution on {-1,1}^n
# Either Rademacher or Mammen probabilities
# B.sample = number of bootstrap samples
# Returns an n x B.sample matrix
#
###########################################################################
###########################################################################

XI.generator <- function(wilddist, n, B.sample){

check.wilddist(wilddist)

if(wilddist == "rademacher"){
  XI <- matrix(sample(c(-1, 1), size = n * B.sample, 
              replace = TRUE, prob = c(.5,.5)), ncol = B.sample, nrow = n)
} 

if(wilddist == "mammen"){
  XI <- matrix(sample(c(-(sqrt(5)-1)/2, (sqrt(5)+1)/2), size = n * B.sample, 
              replace = TRUE, prob = c((sqrt(5)+1)/(2*sqrt(5)), 
              (sqrt(5)-1)/(2*sqrt(5)))), ncol = B.sample, nrow = n) 
}

return(XI)
}

###########################################################################
###########################################################################
#
# Check if a given design matrix X allows a size-controlling critical value
# Bfac = R(X'X)^{-1} X' (cf. the function Bfactor.matrix below)
#
# Returns TRUE if a size-controlling CV exists, FALSE else.
#
###########################################################################
###########################################################################

ExC.check <- function(
X,      #X
R,      #R
as.tol  #tolerance parameter used in checking rank conditions
){

n <- dim(X)[1]
k <- dim(X)[2]
q <- dim(R)[1]
e0 <- matrix(0, nrow = n, ncol = 1)  
XM <- M0lin(X, R)
l <- dim(XM)[2]
qrX <- qr(X)

Bfac <- Bfactor.matrix(qrX, n, R)
  
#check condition for every index such that ei notin span(XM)

v <- c()

for(i in 1:n){
  ei <- e0
  ei[i,1] <- 1

  if( ( rrank(cbind(XM, ei), as.tol) > l ) ){
      A <- Bfac %*% diag(c(qr.resid(qrX, ei)))
      v <- c(v, rrank(A, as.tol))
  }
}

return(!( length(v[v < q])>0 ))
}

###########################################################################
###########################################################################
#
# Input checks to function ``theta''
#
###########################################################################
###########################################################################

check.X.R <- function(X, R){

if( !is.matrix(X) ){
 stop("Invalid 'X' value - must be a matrix")
}

if( !is.matrix(R) ) {
 stop("Invalid 'R' value - must be a matrix")
}

if( dim(X)[2] >= dim(X)[1] ){
 stop("Number of columns of 'X' is not smaller than its number of rows")
}

if( dim(X)[2] == 0 ){
 stop("Number of rows of 'X' must be greater than 0")
}

if( dim(X)[2] != dim(R)[2] ) {
 stop("Matrices 'X' and 'R' have different numbers of columns")
}

if( rrank(X, -1) < dim(X)[2] ) {
 stop("The matrix 'X' is numerically of rank < k")
}

if( rrank(R, -1) < dim(R)[1] ) {
 stop("The matrix 'R' is numerically of rank < q")
}

}

check.r <- function(r, R){
if( !is.vector(r) | !is.numeric(r) | (length(r) != dim(R)[1])){
 stop("Invalid 'r' value - 'r' must be a real vector
 the length of which coincides with the number of rows of R") 
}
}

check.hcmethod <- function(hcmethod){
if(length(hcmethod) != 1 | sum(c(-1:4) == hcmethod) != 1) {
 stop("Invalid 'hcmethod' value - 'hcmethod' must be -1, 0, 1, 2, 3, or 4")
}
}

check.restr.cov <- function(restr.cov){
 if(!is.logical(restr.cov) | length(restr.cov) != 1){
 stop("Invalid 'restr.cov' value - 'restr.cov' must be TRUE or FALSE")
}
}

check.wilddist <- function(wilddist){
if(length(wilddist)!= 1 | sum(c("mammen", "rademacher") == wilddist) != 1) {
 stop("Invalid 'wilddist' value - 'wilddist' must be 'mammen' or 'rademacher' ")
}
}

check.wildmult <- function(wildmult){
if(length(wildmult) != 1 | sum(c(0:4) == wildmult) != 1) {
 stop("Invalid 'wildmult' value - 'wildmult' must be 0, 1, 2, 3, or 4")
}
}

check.wildmult.restr <- function(wildmult.restr){
if(!is.logical(wildmult.restr) | length(wildmult.restr) != 1){
 stop("Invalid 'wildmult.restr' value - 'wildmult.restr' must be TRUE or FALSE")
}
}

check.boot.res.restr <- function(boot.res.restr){
if(!is.logical(boot.res.restr) | length(boot.res.restr) != 1) {
 stop("Invalid 'boot.res.restr' value - 'boot.res.restr' must be TRUE or FALSE")
}
}

check.boot.center.restr <- function(boot.center.restr){
if(!is.logical(boot.center.restr) | length(boot.center.restr) != 1){
 stop("Invalid 'boot.center.restr' value - 'boot.center.restr' must be TRUE or FALSE")
}
}

check.tol <- function(tol){
if(!is.numeric(tol) | length(tol) != 1){
 stop("Invalid 'tol' value - 'tol' must be a real number")
}
if(tol > .1){
warning("Tolerance parameter 'tol' should typically be chosen small, e.g., 1e-07; 
your choice seems unusually large, and might yield wrong results.")
}
if(tol <= 0){
warning("Your non-positive 'tol' parameter was converted to tol = machine epsilon")
}
}

check.as.tol <- function(as.tol){
if(!is.numeric(as.tol) | length(as.tol) != 1){
 stop("Invalid 'as.tol' value - 'as.tol' must be a real number")
}
if(as.tol > .1){
warning("Tolerance parameter 'as.tol' should typically be chosen small, e.g., 1e-07; 
your choice seems unusually large, and might yield wrong results.")
}
if(as.tol <= 0){
warning("Your non-positive 'as.tol' parameter was converted to as.tol = machine epsilon")
}
}

check.in.tol <- function(in.tol){
if(!is.numeric(in.tol) | length(in.tol) != 1){
 stop("Invalid 'in.tol' value - 'in.tol' must be a real number")
}
if(in.tol > .1){
warning("Tolerance parameter 'in.tol' should typically be chosen small, e.g., 1e-07; 
your choice seems unusually large, and might yield wrong results.")
}
if(in.tol <= 0){
warning("Your non-positive 'in.tol' parameter was converted to in.tol = machine epsilon")
}
}

check.comp.meth <- function(comp.meth){
if(length(comp.meth)!= 1 | sum(c("exact", "approximation") == comp.meth) != 1){
 stop("Invalid 'comp.meth' value - 'comp.meth' must be 'exact' or 'approximation'")
}
}

check.Boot.supp <- function(Boot.supp, comp.meth, n){
if(comp.meth == "approximation"){
 if(!is.matrix(Boot.supp) | !is.numeric(Boot.supp) | dim(Boot.supp)[1] != n | dim(Boot.supp)[2] < 1){
  stop("Invalid 'Boot.supp' value - 'Boot.supp' must be a real matrix with n rows")
 }
}
}

check.checks <- function(checks){
if(!is.logical(checks) | length(checks) != 1){
 stop("Invalid 'checks' value - 'checks' must be TRUE or FALSE")
}
}

check.cores <- function(cores){
if( cores%%1 != 0 | cores <= 0) {
 stop("Invalid 'cores' value - 'cores' must be a positive integer")
}
}

#additional check in the boot.pval function

check.y <- function(y, X){
if( !is.matrix(y) ){
 stop("Invalid 'y' value - must be a matrix")
}

if( dim(y)[1] != dim(X)[1]){
 stop("Invalid 'y' value - must be a matrix with n rows")
}
}