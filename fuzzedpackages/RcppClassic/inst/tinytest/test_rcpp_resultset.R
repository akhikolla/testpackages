
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

Rcpp::sourceCpp("cpp/RcppResultSet.cpp")

#test.RcppResultSet.double <- function() {
fun <- double_
expect_equal(fun()[[1]], 1.23456, info = "RcppResultRet.double")

#test.RcppResultSet.int <- function() {
fun <- int_
expect_equal(fun()[[1]], 42, info = "RcppResultSet.int")

#test.RcppResultSet.string <- function() {
fun <- string_
expect_equal(fun()[[1]], "hello unit tests", info = "RcppResultSet.string")

#test.RcppResultSet.double.vector <- function() {
fun <- double_vector
expect_equal(fun()[[1]], c(1.1, 2.2, 3.3), info = "RcppResultSet.double.vector")

#test.RcppResultSet.int.vector <- function() {
fun <- int_vector
expect_equal(fun()[[1]], c(11, 22, 33), info = "RcppResultSet.int.vector")

#test.RcppResultSet.double.matrix <- function() {
fun <- double_matrix
expect_equal(fun()[[1]], matrix(c(1.1, 2.2, 3.3, 4.4), 2, byrow=TRUE), info = "RcppResultSet.double.matrix")

#test.RcppResultSet.int.matrix <- function() {
fun <- int_matrix
expect_equal(fun()[[1]], matrix(c(11, 22, 33, 44), 2, byrow=TRUE), info = "RcppResultSet.int.matrix")

#test.RcppResultSet.RcppDate <- function() {
fun <- RcppDate_
expect_equal(fun()[[1]], as.Date("2000-01-01"), info = "RcppResultSet.RcppDate")

#test.RcppResultSet.RcppDateVector <- function() {
fun <- RcppDateVector_
v <- c(as.Date("2000-01-01"), as.Date("2001-01-01"))
expect_equal(fun(v)[[1]], v, info = "RcppResultSet.RcppDateVector")

#test.RcppResultSet.RcppDatetime <- function() {
fun <- RcppDatetime_
## setting tz = "UTC" because otherwise the format gets set as the tz
posixt <- as.POSIXct("2000-01-01 01:02:03.456", "%Y-%m-%d %H:%M:%OS", tz = "UTC" )
result <- fun(as.numeric(posixt))[[1]]
## RcppDateTime discards the timezone, so we have to set it back
## otherwise the comparison fails on the attributes
attr( result, "tzone") <- "UTC"
expect_true( (result - posixt) == 0.0 , info = "RcppResultSet.RcppDatetime")

#test.RcppResultSet.RcppDatetimeVector <- function() {
fun <- RcppDatetimeVector_
now <- Sys.time()
attr(now, "tzone") <- NULL # no attribute gets set at the C++ level
v <- now + 0:9
expect_true( sum( fun(v)[[1]] - v ) == 0.0 , info = "RcppResultSet.RcppDatetimeVector")

#test.RcppResultSet.RcppStringVector <- function() {
fun <- RcppStringVector_
v <- c("hello", "goodbye")
expect_equal(fun(v)[[1]], v, info = "RcppResultSet.RcppStringVector")

#test.RcppResultSet.std.vector.double <- function() {
fun <- std_vector_double
expect_equal(fun()[[1]], c(1.1, 2.2, 3.3), info = "RcppResultSet.std.vector.double")

#test.RcppResultSet.std.vector.int <- function() {
fun <- std_vector_int
expect_equal(fun()[[1]], c(11, 22, 33), info = "RcppResultSet.std.vector.int")

#test.RcppResultSet.std.vector.std.vector.double <- function() {
fun <- std_vector_std_vector_double
expect_equal(fun()[[1]], matrix(c(1.1, 2.2, 3.3, 1.1, 2.2, 3.3), nrow=2, ncol=3, byrow=TRUE),
             info = "RcppResultSet.std.vector.std.vector.double")

#test.RcppResultSet.std.vector.std.vector.int <- function() {
fun <- std_vector_std_vector_int
expect_equal(fun()[[1]], matrix(c(11, 22, 33, 11, 22, 33), nrow=2, ncol=3, byrow=TRUE),
             info = "RcppResultSet.std.vector.std.vector.int")

#test.RcppResultSet.std.vector.std.vector.string <- function() {
fun <- std_vector_std_vector_string
expect_equal(fun()[[1]], c("hello", "goodbye"), info = "RcppResultSet.std.vector.std.string")

#test.RcppResultSet.RcppVector.int <- function() {
fun <- RcppVector_int
x <- c(11,22,33)
expect_equal(fun(x)[[1]], x, info = "RcppResultSet.RcppVector.int")

#test.RcppResultSet.RcppVector.double <- function() {
fun <- RcppVector_double
x <- c(1.1,2.2,3.3)
expect_equal(fun(x)[[1]], x, info = "RcppResultSet.RcppVector.double")

#test.RcppResultSet.RcppMatrix.int <- function() {
fun <- RcppMatrix_int
x <- matrix(1:9, 3, 3)
expect_equal(fun(x)[[1]], x, info = "RcppResultSet.RcppMatrix.int")

#test.RcppResultSet.RcppMatrix.double <- function() {
fun <- RcppMatrix_double
x <- matrix(1.1*c(1:9), 3, 3)
expect_equal(fun(x)[[1]], x, info = "RcppResultSet.RcppMatrix.double")

#test.RcppResultSet.RcppFrame <- function() {
fun <- RcppFrame_
x <- data.frame(x=1:9, y=LETTERS[1:9], z=sample(c(TRUE,FALSE), 9, replace=TRUE))
expect_equal( as.data.frame(fun(x)[[1]]), x, info = "RcppResultSet.RcppFrame")

#test.RcppResultSet.SEXP <- function() {
fun <- SEXP_
x <- list(foo=1.23, bar=123, glim="glom")
expect_equal( fun(x)[[1]], x, info = "RcppResultSet.SEXP")

#test.RObject.asStdVectorIntResultsSet <- function(){
funx <- vector_int_rs
expect_equal(funx(2:5), 2:5*2L, info = "as<std::vector<int> >(integer) via RcppResultSet")
expect_equal(funx(2:5+.1), 2:5*2L, info = "as<std::vector<int> >(numeric) via RcppResultSet")
expect_equal(funx(as.raw(2:5)), 2:5*2L, info = "as<std::vector<int> >(raw) via RcppResultSet")
expect_error(funx("foo"), info = "as<std::vector<int> >(character) -> exception")
