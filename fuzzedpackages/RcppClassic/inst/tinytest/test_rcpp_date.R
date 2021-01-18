
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

Rcpp::sourceCpp("cpp/RcppDate.cpp")

Sys.setenv("TZ"="UTC")          # to ensure localtime is GMT

#test.RcppDate.get.functions <- function() {
fun <- get_functions
expect_equal(fun(), list(month=12, day=31, year=1999, julian=10956), info = "RcppDate.get.functions")

#test.RcppDate.operators <- function() {
fun <- RcppDate_operators
expect_equal(fun(), list(diff=1, bigger=TRUE, smaller=FALSE, equal=FALSE, ge=TRUE, le=FALSE),
            info = "RcppDate.operators")

#test.RcppDate.wrap <- function() {
fun <- RcppDate_wrap
expect_equal(fun(), as.Date("1999-12-31"), info = "RcppDate.wrap")

#test.RcppDatetime.get.functions <- function() {
#    fun <- .rcpp.RcppDate$RcppDatetime_functions
#    expect_equal(#fun(as.numeric(as.POSIXct("2001-02-03 01:02:03.123456", tz="UTC"))),
#                fun(981162123.123456),
#                list(year=2001, month=2, day=3, wday=6, hour=1, minute=2, second=3, microsec=123456),
#                info = "RcppDate.get.functions")
#}

#test.RcppDatetime.operators <- function() {
fun <- RcppDatetime_operators
expect_equal(fun(), #as.numeric(as.POSIXct("2001-02-03 01:02:03.123456", tz="UTC"))),
             list(diff=3600, bigger=TRUE, smaller=FALSE, equal=FALSE, ge=TRUE, le=FALSE),
             info = "RcppDatetime.operators")

#test.RcppDatetime.wrap <- function() {
fun <- RcppDatetime_wrap
expect_equal(as.numeric(fun()), as.numeric(as.POSIXct("2001-02-03 01:02:03.123456", tz="UTC")),
            info = "RcppDatetime.wrap")
