
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

Rcpp::sourceCpp("cpp/RcppMisc.cpp")

#test.RcppFrame <- function() {
fun <- RcppFrameFunc
dframe <- data.frame(fun()[[1]]) ## needs a data.frame() call on first list elem
expect_equal(dframe, data.frame(A=c(1.23,4.56), B=c(42,21), C=c(FALSE,TRUE)), info = "RcppFrame")

#test.RcppList <- function() {
fun <- RcppListFunc
expect_equal(fun(), list(foo=1L, bar=2, biz="xyz"), info="RcppList")

#test.RcppParams.Double <- function() {
fun <- RcppParams_Double
expect_equal(fun(list(val=1.234)), 2*1.234, info="RcppParams.getDoubleValue")

#test.RcppParams.Int <- function() {
fun <- RcppParams_Int
expect_equal(fun(list(val=42)), 2*42, info="RcppParams.getIntValue")

#test.RcppParams.String <- function() {
fun <- RcppParams_String
expect_equal(fun(list(val="a test string")), "a test stringa test string",
             info = "RcppParams.getStringValue")

#test.RcppParams.Bool <- function() {
fun <- RcppParams_Bool
expect_equal(fun(list(val=FALSE)), FALSE, info = "RcppParams.getBoolValue")


#test.RcppParams.Date <- function() {
fun <- RcppParams_Date
expect_equal(fun(list(val=as.Date("2000-01-01")))[[1]], as.Date("2000-01-01"),
             info = "RcppParams.getDateValue")

#test.RcppParams.Datetime <- function() {
fun <- RcppParams_Datetime
posixt <- as.POSIXct(strptime("2000-01-02 03:04:05.678", "%Y-%m-%d %H:%M:%OS"))
attr(posixt, "tzone") <- NULL    ## because we don't set a tzone attribute in C++
result <- fun(list(val=posixt))[[1]]
expect_true( (result-posixt) == 0.0 , info = "RcppParams.getDatetimeValue")
