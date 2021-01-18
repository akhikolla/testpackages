
// Copyright (C) 2010 - 2019  Dirk Eddelbuettel and Romain Francois
//
// This file is part of RcppClassic.
//
// RcppClassic is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// RcppClassic is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with RcppClassic.  If not, see <http://www.gnu.org/licenses/>.

// NB Files updated in 2019 when switching to tinytest; usage of RcppClassic
// is now more idiomatic Rcpp use given that Rcpp (and hence Rcpp Attribute)
// are available.
//
// For more RcppClassic usage see e.g. the RcppClassicExamples packages.

#include <RcppClassic.h>

// [[Rcpp::depends(RcppClassic)]]

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List get_functions() {
    RcppDate dt = RcppDate(12,31,1999);
    RcppResultSet rs;
    rs.add("month", dt.getMonth());
    rs.add("day",   dt.getDay());
    rs.add("year",  dt.getYear());
    rs.add("julian",dt.getJulian());
    return rs.getReturnList();
}

// [[Rcpp::export]]
Rcpp::List RcppDate_operators() {
    RcppDate d1 = RcppDate(12,31,1999);
    RcppDate d2 = d1 + 1;
    RcppResultSet rs;
    rs.add("diff",    d2 - d1);
    rs.add("bigger",  d2 > d1);
    rs.add("smaller", d2 < d1);
    rs.add("equal",   d2 == d1);
    rs.add("ge",      d2 >= d1);
    rs.add("le",      d2 <= d1);
    return rs.getReturnList();
}

// [[Rcpp::export]]
RcppDate RcppDate_wrap() {
    RcppDate dt = RcppDate(12,31,1999);
    return dt;                  // wrap() implicit
}

// [[Rcpp::export]]
Rcpp::List RcppDatetime_functions(Rcpp::NumericVector x) {
    RcppDatetime dt = RcppDatetime(x);
    RcppResultSet rs;
    rs.add("year",     dt.getYear());
    rs.add("month",    dt.getMonth());
    rs.add("day",      dt.getDay());
    rs.add("wday",     dt.getWeekday());
    rs.add("hour",     dt.getHour());
    rs.add("minute",   dt.getMinute());
    rs.add("second",   dt.getSecond());
    rs.add("microsec", dt.getMicroSec());
    return rs.getReturnList();
}

// [[Rcpp::export]]
Rcpp::List RcppDatetime_operators() {
    RcppDatetime d1 = RcppDatetime(946774923.123456);
    //RcppDatetime d1 = RcppDatetime(1152338523.456789);
    // as.POSIXct("2006-07-08 01:02:03.456789")
    RcppDatetime d2 = d1 + 60*60;
    RcppResultSet rs;
    rs.add("diff",    d2 - d1);
    rs.add("bigger",  d2 > d1);
    rs.add("smaller", d2 < d1);
    rs.add("equal",   d2 == d1);
    rs.add("ge",      d2 >= d1);
    rs.add("le",      d2 <= d1);
    return rs.getReturnList();
}

// [[Rcpp::export]]
RcppDatetime RcppDatetime_wrap() {
    RcppDatetime dt = RcppDatetime(981162123.123456);
    return dt;                  // wrap() implicit
}
