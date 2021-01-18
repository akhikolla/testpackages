// RcppList.h: RcppClassic R/C++ interface class library -- 'list' type support
//
// Copyright (C) 2009 - 2018  Dirk Eddelbuettel
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

#ifndef RcppList_h
#define RcppList_h

#include <RcppCommon.h>

// This class was first used in the RProtoBuf project (currently on r-forge only)
// but is more generally useful and hence moved over here
class RcppList {
public:
    RcppList(void);
    ~RcppList();
    void setSize(int size);

    // defined later because it needs wrap
    template <typename T>
    void append( const std::string& name, const T& value );

    void clearProtectionStack();
    SEXP getList(void) const;

protected:
    friend class RcppResultSet;

private:
    SEXP listArg;
    int listSize, currListPosn, numProtected;
    std::vector<std::string> names;
};

namespace Rcpp{
    template<> inline SEXP wrap<RcppList>( const RcppList& x ){
        return x.getList( ) ;
    }
}

#endif
