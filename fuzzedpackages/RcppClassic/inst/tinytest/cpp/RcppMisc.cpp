
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
Rcpp::List RcppFrameFunc() {
    std::vector<std::string> names;
    names.push_back("A");
    names.push_back("B");
    names.push_back("C");
    RcppFrame fr(names);

    std::vector<ColDatum> colDatumVector(3);
    colDatumVector[0].setDoubleValue(1.23);
    colDatumVector[1].setIntValue(42);
    colDatumVector[2].setLogicalValue(0);
    fr.addRow(colDatumVector);

    colDatumVector[0].setDoubleValue(4.56);
    colDatumVector[1].setIntValue(21);
    colDatumVector[2].setLogicalValue(1);
    fr.addRow(colDatumVector);

    RcppResultSet rs;
    rs.add("data.frame", fr);
    return rs.getReturnList();
}

// [[Rcpp::export]]
Rcpp::List RcppListFunc() {
    RcppList l;
    l.setSize(3);
    l.append("foo", 1);
    l.append("bar", 2.0);
    l.append("biz", "xyz");
    return l.getList();
}

// [[Rcpp::export]]
double RcppParams_Double(RcppParams p) {
    double y = 2 * p.getDoubleValue("val");
    return y;
}

// [[Rcpp::export]]
int RcppParams_Int(RcppParams p) {
    int y = 2 * p.getIntValue("val");
    return y;
}

// [[Rcpp::export]]
std::string RcppParams_String(RcppParams p) {
    std::string y = p.getStringValue("val");
    y = y + y; // trivial string operation
    return y;
}

// [[Rcpp::export]]
bool RcppParams_Bool(RcppParams p) {
    bool y = p.getBoolValue("val");
    return y;
}

// [[Rcpp::export]]
Rcpp::List RcppParams_Date(RcppParams p) {
    RcppDate y = p.getDateValue("val");
    RcppResultSet rs;
    rs.add("date", y);
    return rs.getReturnList();
}

// [[Rcpp::export]]
Rcpp::List RcppParams_Datetime(RcppParams p) {
    RcppDatetime y = p.getDatetimeValue("val");
    RcppResultSet rs;
    rs.add("datetime", y);
    return rs.getReturnList();
}
