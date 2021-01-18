// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// RcppDatetimeeVector.h: RcppClassic R/C++ interface class library -- Datetime vector support
//
// Copyright (C) 2008 - 2009 Dirk Eddelbuettel
// Copyright (C) 2010	     Dirk Eddelbuettel and Romain Francois
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

#ifndef RcppDatetimeVector_h
#define RcppDatetimeVector_h

#include <RcppCommon.h>
#include <classic/RcppDatetime.h>

class RcppDatetimeVector {
public:
    RcppDatetimeVector(SEXP vec);
    RcppDatetimeVector(int n);
    ~RcppDatetimeVector() {};
    const RcppDatetime &operator()(int i) const;
    RcppDatetime &operator()(int i);
    int size() const;
private:
    std::vector<RcppDatetime> v;
};

#endif
