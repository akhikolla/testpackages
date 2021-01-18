// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
//
// RcppClassic.h: R/C++ interface class library
//
// Copyright (C) 2005 - 2006 Dominick Samperi
// Copyright (C) 2008 - 2009 Dirk Eddelbuettel
// Copyright (C) 2009 - 2012 Dirk Eddelbuettel and Romain Francois
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

#ifndef RcppClassic_hpp
#define RcppClassic_hpp

#include <RcppCommon.h>
#include <classic/classic.h>
#include <Rcpp.h>
#include <classic/classic_backward.h>

namespace Rcpp{
    namespace internal{
        SEXP getPosixClasses() ;
        SEXP new_posixt_object( double d) ;
        SEXP new_date_object( double d) ;
    }
    
}

#endif
