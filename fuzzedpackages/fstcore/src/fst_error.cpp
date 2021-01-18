/*
 fstcore - R bindings to the fstlib library

 Copyright (C) 2017-present, Mark AJ Klik

 This file is part of the fstcore R package.

 The fstcore R package is free software: you can redistribute it and/or modify it
 under the terms of the GNU Affero General Public License version 3 as
 published by the Free Software Foundation.

 The fstcore R package is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License
 for more details.

 You should have received a copy of the GNU Affero General Public License along
 with the fstcore R package. If not, see <http://www.gnu.org/licenses/>.

 You can contact the author at:
 - fstcore R package source repository : https://github.com/fstpackage/fstcore
*/


#include <Rcpp.h>


SEXP fst_error(const char* error_message)
{
  SEXP fst_error;

  PROTECT(fst_error = Rf_mkString(error_message));

  Rf_classgets(fst_error, Rf_mkString("fst_error"));

  UNPROTECT(1);

  return fst_error;
}
