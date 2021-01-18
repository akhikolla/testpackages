// RcppList__backward.h: RcppClassic R/C++ interface class library --
//
// Copyright (C) 2010 - 2018  Dirk Eddelbuettel and Romain Francois
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

#ifndef RcppList__backward_h
#define RcppList__backward_h

template <typename T>
void RcppList::append( const std::string& name, const T& value ) {
    if (currListPosn < 0 || currListPosn >= listSize)
        throw std::range_error("RcppList::append(): list posn out of range");

    SET_VECTOR_ELT(listArg, currListPosn++, Rcpp::wrap(value) );
    names.push_back(name);
}

#endif
