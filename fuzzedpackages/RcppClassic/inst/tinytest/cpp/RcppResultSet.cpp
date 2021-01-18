
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
Rcpp::List double_() {
    double y = 1.23456;
    RcppResultSet rs;
    rs.add("foo", y);
    return rs.getReturnList();
}

// [[Rcpp::export]]
Rcpp::List int_() {
    int y = 42;
    RcppResultSet rs;
    rs.add("foo", y);
    return rs.getReturnList();
}

// [[Rcpp::export]]
Rcpp::List string_() {
    std::string y = "hello unit tests";
    RcppResultSet rs;
    rs.add("foo", y);
    return rs.getReturnList();
}

// [[Rcpp::export]]
Rcpp::List double_vector() {
    double y[3] = { 1.1, 2.2, 3.3 };
    RcppResultSet rs;
    rs.add("foo", y, 3);
    return rs.getReturnList();
}

// [[Rcpp::export]]
Rcpp::List int_vector() {
    int y[3] = { 11, 22, 33 };
    RcppResultSet rs;
    rs.add("foo", y, 3);
    return rs.getReturnList();
}

// [[Rcpp::export]]
Rcpp::List double_matrix() {
    double r1[2] = { 1.1, 2.2 };
    double r2[2] = { 3.3, 4.4 };
    double *y[2] = { r1, r2 };
    RcppResultSet rs;
    rs.add("foo", y, 2, 2);
    return rs.getReturnList();
}

// [[Rcpp::export]]
Rcpp::List int_matrix() {
    int r1[2] = { 11, 22 };
    int r2[2] = { 33, 44 };
    int *y[2] = { r1, r2 };
    RcppResultSet rs;
    rs.add("foo", y, 2, 2);
    return rs.getReturnList();
}

// [[Rcpp::export]]
Rcpp::List RcppDate_() {
    RcppDate y(01,01,2000); // silly North American mon-day-year
    RcppResultSet rs;
    rs.add("foo", y);
    return rs.getReturnList();
}

// [[Rcpp::export]]
Rcpp::List RcppDateVector_(SEXP x) {
    RcppDateVector y(x);
    RcppResultSet rs;
    rs.add("foo", y);
    return rs.getReturnList();
}

// [[Rcpp::export]]
Rcpp::List RcppDatetime_(RcppDatetime y) {
    RcppResultSet rs;
    rs.add("foo", y);
    return rs.getReturnList();
}

// [[Rcpp::export]]
Rcpp::List RcppDatetimeVector_(RcppDatetimeVector y) {
    RcppResultSet rs;
    rs.add("foo", y);
    return rs.getReturnList();
}

// [[Rcpp::export]]
Rcpp::List RcppStringVector_(RcppStringVector y) {
    RcppResultSet rs;
    rs.add("foo", y);
    return rs.getReturnList();
}

// [[Rcpp::export]]
Rcpp::List std_vector_double() {
    std::vector<double> y;
    y.push_back(1.1);
    y.push_back(2.2);
    y.push_back(3.3);
    RcppResultSet rs;
    rs.add("foo", y);
    return rs.getReturnList();
}

// [[Rcpp::export]]
Rcpp::List std_vector_int() {
    std::vector<int> y;
    y.push_back(11);
    y.push_back(22);
    y.push_back(33);
    RcppResultSet rs;
    rs.add("foo", y);
    return rs.getReturnList();
}

// [[Rcpp::export]]
Rcpp::List std_vector_std_vector_double() {
    std::vector<double> yy;
    yy.push_back(1.1);
    yy.push_back(2.2);
    yy.push_back(3.3);
    std::vector< std::vector<double> > y;
    y.push_back(yy);
    y.push_back(yy);
    RcppResultSet rs;
    rs.add("foo", y);
    return rs.getReturnList();
}

// [[Rcpp::export]]
Rcpp::List std_vector_std_vector_int() {
    std::vector<int> yy;
    yy.push_back(11);
    yy.push_back(22);
    yy.push_back(33);
    std::vector< std::vector<int> > y;
    y.push_back(yy);
    y.push_back(yy);
    RcppResultSet rs;
    rs.add("foo", y);
    return rs.getReturnList();
}

// [[Rcpp::export]]
Rcpp::List std_vector_std_vector_string() {
    std::string a("hello");
    std::string b("goodbye");
    std::vector< std::string > y;
    y.push_back(a);
    y.push_back(b);
    RcppResultSet rs;
    rs.add("foo", y);
    return rs.getReturnList();
}

// [[Rcpp::export]]
Rcpp::List RcppVector_int(RcppVector<int> y) {
    RcppResultSet rs;
    rs.add("foo", y);
    return rs.getReturnList();
}

// [[Rcpp::export]]
Rcpp::List RcppVector_double(RcppVector<double> y) {
    RcppResultSet rs;
    rs.add("foo", y);
    return rs.getReturnList();
}

// [[Rcpp::export]]
Rcpp::List RcppMatrix_int(RcppMatrix<int> y) {
    RcppResultSet rs;
    rs.add("foo", y);
    return rs.getReturnList();
}

// [[Rcpp::export]]
Rcpp::List RcppMatrix_double(RcppMatrix<double> y) {
    RcppResultSet rs;
    rs.add("foo", y);
    return rs.getReturnList();
}

// [[Rcpp::export]]
Rcpp::List RcppFrame_(RcppFrame y) {
    RcppResultSet rs;
    rs.add("", y);
    return rs.getReturnList();
}

// [[Rcpp::export]]
Rcpp::List SEXP_(SEXP x) {
    RcppResultSet rs;
    rs.add("", x, false);
    return rs.getReturnList();
}

// [[Rcpp::export]]
SEXP vector_int_rs(std::vector<int> iv) {
    for (size_t i=0; i<iv.size(); i++) {
        iv[i] = 2*iv[i];
    }
    RcppResultSet rs;
    rs.add("", iv);
    return(rs.getSEXP());
}
