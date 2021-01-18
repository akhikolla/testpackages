#               
#    Copyright (C) 2016  David Preinerstorfer
#    david.preinerstorfer@econ.au.dk
#
#    This file contains the internal auxiliary functions:
#    Bfactor.matrix; gen.start; wm; and several functions for input checking
#    
#    This file is a part of acrt.
#
#    acrt is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details. A copy may be obtained at
#    http://www.r-project.org/Licenses/

#Compute the matrix R(X'X)^(-1)X' given the qr decomposition of a design matrix
#(n times k with rank k) and a restriction matrix R (q times k)

Bfactor.matrix <- function(qrX, R){
factor.tmp <- qr.R(qrX)
factor.tmp <- backsolve(factor.tmp, diag(dim(factor.tmp)[1]))
Bfactor <- R %*% factor.tmp %*% t(qr.Q(qrX))
return(Bfactor)
}

#Draw uniformly from the stationarity region of AR(ar.order) processes 
#parameterized by its partial autocovariance function; 
#Based on results in: MC Jones, Applied Statistics, 1987. 

gen.start <- function(ar.order){
v <- 1:ar.order
s1 <- floor(.5 * (v + 1))
s2 <- floor(.5 * v) + 1
s.val <- c()
for(i in 1:ar.order){
s.val <- c(s.val, rbeta(1, s1[i], s2[i]))
}
return(2*s.val - 1)
}

#Generate an n x n dimensional weights matrix from a kernel function
#and a bandwidth parameter using the kweights function from the sandwich package
#by A. Zeileis (Achim Zeileis (2004). Econometric Computing with HC and HAC 
#Covariance Matrix Estimators. Journal of Statistical Software 11(10), 1-17.)

wm <- function(n, bandwidth, ker = "Bartlett"){
toeplitz(sandwich::kweights(0:(n-1)/bandwidth, ker))
}

#input checks

check.alpha <- function(alpha){
if( !is.numeric(alpha) | alpha <= 0 | alpha >= 1 ) {
 stop("Invalid 'alpha' value - 'alpha' must be in the interval (0,1)")
}
}

check.bandwidth <- function(bandwidth){
if( !is.numeric(bandwidth) | bandwidth <= 0) {
 stop("Invalid 'bandwidth' value - 'bandwidth' must be a real number > 0")
}
}

check.Eicker <- function(Eicker){
if( !is.logical(Eicker) ) {
 stop("Invalid 'Eicker' value - must be logical")
}
}

check.X.R.order <- function(X, R, ar.order.max){

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
 stop("Matrices 'X' and 'R' have different number of columns")
}

if( qr(X)$rank < dim(X)[2] ) {
 stop("The matrix 'X' is numerically of rank < k")
}

if( qr(R)$rank < dim(R)[1] ) {
 stop("The matrix 'R' is numerically of rank < q")
}

if( ar.order.max%%1 != 0 | ar.order.max < 0 ){
 stop("Invalid 'ar.order.max' value - 'ar.order.max' must be an integer >= 0")
}

}

check.N.M <- function(N0, N1, N2, Mp, M1, M2){

if( N0%%1 != 0 | N0 < 0 ){
 stop("Invalid 'N0' value - 'N0' must be a positive integer")
}

if( N1%%1 != 0 | N1 < 0 ){
 stop("Invalid 'N1' value - 'N1' must be a positive integer")
}

if( N2%%1 != 0 | N2 < 0 ){
 stop("Invalid 'N2' value - 'N2' must be a positive integer")
}

if( Mp%%1 != 0 | Mp < 0 ){
 stop("Invalid 'Mp' value - 'Mp' must be a positive integer")
}

if( M1%%1 != 0 | M1 < 0 ){
 stop("Invalid 'M1' value - 'M1' must be a positive integer")
}

if( M2%%1 != 0 | M2 < 0 ){
 stop("Invalid 'M2' value - 'M2' must be a positive integer")
}

if( N0 > N1 ) {
 warning("'N1' should be greater than 'N0'")
}

if( N1 > N2 ) {
 warning("'N2' should be greater than 'N1'")
}

if( M1 > Mp ) {
 stop("Invalid 'M1' value - 'M1' can not be greater than 'Mp'")
}

if( M2 > M1 ) {
 stop("Invalid 'M2' value - 'M2' can not be greater than 'M1'")
}
                
}

check.cores <- function(cores){
if( cores%%1 != 0 | cores <= 0) {
 stop("Invalid 'cores' value - 'cores' must be a positive integer")
} 
}

check.margin <- function(margin, ar.order.max){
if( !is.vector(margin) | length(margin) != ar.order.max |  min(margin) <= 0 | 
max(margin) > 1 ) {
 stop("Invalid 'margin' value - 'margin' must be a vector of the same length as
 'ar.order.max' and with coordinates greater than 0 and not exceeding 1")
} 
}

check.C <- function(C){
if( !is.numeric(C) | length(C)!=1 | C <= 0) {
 stop("Invalid 'C' value - 'C' must be a positive real number")
}
}

check.y <- function(y, X){
if( !is.numeric(y) | !is.matrix(y) | dim(y)[2] == 0 | dim(y)[1] != dim(X)[1]){
 stop("Invalid 'y' value - 'y' must either be a real vector of length 
 dim(X)[1], or a real matrix with dim(X)[1] rows with more than 0 columns")
}
}

check.r <- function(r, R){
if( !is.vector(r) | !is.numeric(r)){
 stop("Invalid 'r' value - 'r' must be a real vector") 
}

if( length(r) != dim(R)[1] ){
 stop("Invalid 'r' dimension - the length of 'r' must coincide with the 
 number of rows of the matrix 'R'") 
}

}

check.ker <- function(ker){
if( !(ker %in% c("Bartlett", "Parzen", "Quadratic Spectral")) ){
 stop("Invalid 'ker' value - 'ker' must be one of the following kernels:
 'Bartlett', 'Parzen', 'Quadratic Spectral'") 
} 
}

check.method <- function(opt.method.1, opt.method.2){
if( !(opt.method.1 %in% c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN")) ){
 stop("Invalid optimization method in Stage 1 - opt.method.1 must be one of
 the following methods: 'Nelder-Mead', 'BFGS', 'CG, 'L-BFGS_B', 'SANN'.")
}

if( !opt.method.2 %in% c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN") ){
 stop("Invalid optimization method in Stage 2 - opt.method.2 must be one of
 the following methods: 'Nelder-Mead', 'BFGS', 'CG, 'L-BFGS_B', 'SANN'.")
}

}