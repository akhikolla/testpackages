
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
Rcpp::List RcppMatrix_int(RcppMatrix<int> m) {
    RcppResultSet rs;
    rs.add("dim1",  m.getDim1());
    rs.add("dim2",  m.getDim2());
    rs.add("rows",  m.rows());
    rs.add("cols",  m.cols());
    rs.add("p22",   m(1,1));
    std::vector<std::vector<int> > mm = m.stlMatrix();
    rs.add("m",     mm);
    return rs.getReturnList();
}

// [[Rcpp::export]]
Rcpp::List RcppMatrix_double(RcppMatrix<double> m) {
    RcppResultSet rs;
    rs.add("dim1",  m.getDim1());
    rs.add("dim2",  m.getDim2());
    rs.add("rows",  m.rows());
    rs.add("cols",  m.cols());
    rs.add("p22",   m(1,1));
    std::vector<std::vector<double> > mm = m.stlMatrix();
    rs.add("m",     mm);
    return rs.getReturnList();
}

// [[Rcpp::export]]
Rcpp::List RcppMatrix_double_na_nan(RcppMatrix<double> m) {
    RcppResultSet rs;
    rs.add("na_21",  R_IsNA(m(1,0)));
    rs.add("na_22",  R_IsNA(m(1,1)));
    rs.add("nan_31", R_IsNaN(m(2,0)));
    rs.add("nan_32", R_IsNaN(m(2,1)));
    return rs.getReturnList();
}

// [[Rcpp::export]]
Rcpp::List RcppMatrixView_int(RcppMatrixView<int> m) {
    RcppResultSet rs;
    rs.add("dim1",  m.dim1());
    rs.add("dim2",  m.dim2());
    rs.add("rows",  m.rows());
    rs.add("cols",  m.cols());
    rs.add("p22",   m(1,1));
    return rs.getReturnList();
}

// [[Rcpp::export]]
Rcpp::List RcppMatrixView_double(RcppMatrixView<double> m) {
    RcppResultSet rs;
    rs.add("dim1",  m.dim1());
    rs.add("dim2",  m.dim2());
    rs.add("rows",  m.rows());
    rs.add("cols",  m.cols());
    rs.add("p22",   m(1,1));
    return rs.getReturnList();
}

// [[Rcpp::export]]
Rcpp::List RcppVector_int(RcppVector<int> m) {
    RcppResultSet rs;
    rs.add("size",  m.size());
    rs.add("p2",    m(1));
    std::vector<int> v = m.stlVector();
    rs.add("v",     v);
    return rs.getReturnList();
}

// [[Rcpp::export]]
Rcpp::List RcppVector_double(RcppVector<double> m) {
    RcppResultSet rs;
    rs.add("size", m.size());
    rs.add("p2",   m(1));
    std::vector<double> v = m.stlVector();
    rs.add("v",     v);
    return rs.getReturnList();
}

// [[Rcpp::export]]
Rcpp::List RcppVector_double_na_nan(RcppVector<double> m) {
    RcppResultSet rs;
    rs.add("na_2",  R_IsNA(m(1)));
    rs.add("na_3",  R_IsNA(m(2)));
    rs.add("nan_4", R_IsNaN(m(3)));
    rs.add("nan_5", R_IsNaN(m(4)));
    return rs.getReturnList();
}

// [[Rcpp::export]]
Rcpp::List RcppVectorView_int(RcppVectorView<int> m) {
    RcppResultSet rs;
    rs.add("size",  m.size());
    rs.add("p2",    m(1));
    return rs.getReturnList();
}

// [[Rcpp::export]]
Rcpp::List RcppVectorView_double(RcppVectorView<double> m) {
    RcppResultSet rs;
    rs.add("size", m.size());
    rs.add("p2",   m(1));
    return rs.getReturnList();
}

// [[Rcpp::export]]
Rcpp::List RcppStringVector_classic(RcppStringVector s) {
    RcppResultSet rs;
    rs.add("string", s);
    return rs.getReturnList();
}

// [[Rcpp::export]]
RcppStringVector RcppStringVector_wrap(RcppStringVector s) {
    return s;
}

// [[Rcpp::export]]
RcppStringVector RcppStringVector_begin(RcppStringVector s) {
    return wrap(*s.begin());
}

// [[Rcpp::export]]
RcppStringVector RcppStringVector_end(RcppStringVector s) {
    return wrap(s(s.size()-1));
}
