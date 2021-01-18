/*
 *  bigmemory: an R package for managing massive matrices using C,
 *  with support for shared memory.
 *
 *  Copyright (C) 2008 John W. Emerson and Michael J. Kane
 *
 *  This file is part of bigmemory.
 *
 *  bigmemory is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published
 *  by the Free Software Foundation; either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 */

#include <limits.h>
#include <math.h>
#include "../inst/synchronicity/util.h"

vector<string> RChar2StringVec( SEXP charVec )
{
  vector<string> ret( Rf_length(charVec) );
  unsigned long i;
  for (i=0; i < ret.size(); ++i)
  {
    ret[i] = string(CHAR(STRING_ELT(charVec, i)));
  }
  return ret;
}

vector<string> RChar2StringVec( SEXP charVec, 
  const vector<unsigned long> &indices )
{
  vector<string> ret( indices.size() );
  unsigned long i;
  for (i=0; i < indices.size(); ++i)
  {
    ret[i] = string(CHAR(STRING_ELT(charVec, indices[i]-1)));
  }
  return ret;
}

std::string RChar2String(SEXP str)
{
  return string(CHAR(STRING_ELT(str, 0)));
}

SEXP StringVec2RChar( const vector<string> &strVec )
{
  if (strVec.empty())
    return R_NilValue;
  SEXP ret = Rf_protect(Rf_allocVector(STRSXP, strVec.size()));
  unsigned long i;
  for (i=0; i < strVec.size(); ++i)
  {
    SET_STRING_ELT(ret, i, Rf_mkChar(strVec[i].c_str()));
  }
  Rf_unprotect(1);
  return ret;
}

SEXP String2RChar(const std::string &str)
{
  SEXP ret = Rf_protect(Rf_allocVector(STRSXP, 1));
  SET_STRING_ELT(ret, 0, Rf_mkChar(str.c_str()));
  Rf_unprotect(1);
  return ret;
}
