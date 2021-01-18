
## Copyright (C) 2010 - 2019  Dirk Eddelbuettel and Romain Francois
##
## This file is part of RcppClassic.
##
## RcppClassic is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
##
## RcppClassic is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with RcppClassic.  If not, see <http://www.gnu.org/licenses/>.

Rcpp::sourceCpp("cpp/RcppMatrix.cpp")

#test.RcppMatrix.int <- function() {
funx <- RcppMatrix_int
M <- matrix(1:6, 2, 3, byrow=TRUE)
expect_equal(funx(M), list(dim1=2, dim2=3, rows=2, cols=3, p22=5, m=M), info = "RcppMatrix.int")


#test.RcppMatrix.double <- function() {
funx <- RcppMatrix_double
M <- matrix(1:6,2,3,byrow=TRUE)
expect_equal(funx(M), list(dim1=2, dim2=3, rows=2, cols=3, p22=5, m=M), info = "RcppMatrix.double")


#test.RcppMatrix.double.na.nan <- function() {
funx <- RcppMatrix_double_na_nan
M <- matrix(1:6,3,2,byrow=TRUE)
M[2,1] <- NA
M[3,1] <- NaN
expect_equal(funx(M),
             list(na_21=1, na_22=0, nan_31=1, nan_32=0),
             info = "RcppMatrix.double.na.nan")

#test.RcppMatrixView.int <- function() {
funx <- RcppMatrixView_int
expect_equal(funx(matrix(1:6,2,3,byrow=TRUE)),
             list(dim1=2, dim2=3, rows=2, cols=3, p22=5),
             info = "RcppViewMatrix.int")


#test.RcppMatrixView.double <- function() {
funx <- RcppMatrixView_double
expect_equal(funx(matrix(1.0*(1:6),2,3,byrow=TRUE)),
             list(dim1=2, dim2=3, rows=2, cols=3, p22=5),
             info = "RcppMatrixView.double")


#test.RcppVector.int <- function() {
funx <- RcppVector_int
expect_equal(funx(c(1:6)), list(size=6, p2=2, v=c(1:6)), info="RcppVector.int")

#test.RcppVector.double <- function() {
funx <- RcppVector_double
expect_equal(funx(c(1:6)), list(size=6, p2=2, v=c(1:6)), info="RcppVector.double")

#test.RcppVector.double.na.nan <- function() {
funx <- RcppVector_double_na_nan
x <- 1:6
x[2] <- NA
x[4] <- NaN
expect_equal(funx(x), list(na_2=1, na_3=0, nan_4=1, nan_5=0), info = "RcppMatrix.double.na.nan")

#test.RcppVectorView.int <- function() {
funx <- RcppVectorView_int
expect_equal(funx(c(1:6)), list(size=6, p2=2), info="RcppVectorView.int")

#test.RcppVectorView.double <- function() {
funx <- RcppVectorView_double
expect_equal(funx(1.0*c(1:6)), list(size=6, p2=2), info="RcppVectorView.double")

#test.RcppStringVector.classic <- function() {
fun <- RcppStringVector_classic
sv <- c("tic", "tac", "toe")
expect_equal(fun(sv), list(string=sv), info = "RcppStringVector.classic")

#test.RcppStringVector.wrap <- function() {
fun <- RcppStringVector_wrap
sv <- c("tic", "tac", "toe")
expect_equal(fun(sv), sv, info = "RcppStringVector.wrap")

#test.RcppStringVector.begin <- function() {
fun <- RcppStringVector_begin
sv <- c("tic", "tac", "toe")
expect_equal(fun(sv), sv[1], info = "RcppStringVector.begin")

#test.RcppStringVector.end <- function() {
fun <- RcppStringVector_end
sv <- c("tic", "tac", "toe")
expect_equal(fun(sv), sv[3], info = "RcppStringVector.begin")
