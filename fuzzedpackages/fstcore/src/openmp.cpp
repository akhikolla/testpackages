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

#include <interface/openmphelper.h>

#ifdef _OPENMP
#include <pthread.h>
#endif


static int FstThreadsAtFork = 1;
static bool RestoreAfterFork = true;


SEXP getnrofthreads()
{
  return Rf_ScalarInteger(GetFstThreads());
}


int setnrofthreads(SEXP nrOfThreads)
{
  SEXP intVal = Rf_coerceVector(nrOfThreads, INTSXP);

  if (!Rf_isInteger(intVal) || Rf_length(intVal) != 1 || INTEGER(intVal)[0] < 0)
  {
    // catches NA too since NA is -ve
    Rcpp::stop("Argument to threads_fst must be a single integer >= 0, default 0 is recommended to use all cores");
  }

  return ThreadsFst(*INTEGER(intVal));
}


void restore_after_fork(bool restore)
{
  RestoreAfterFork = restore;
}

// Folowing code was adopted from package data.table:
// TODO: implement in fstcore to guard against forking issues with OpenMP

void when_fork()
{
  // use call to fstlib that does not call any OpenMP library functions internally
  FstThreadsAtFork = GetThreads();

  SetThreads(1);
}


void when_unfork()
{
  // use call to fstlib that does not call any OpenMP library functions internally
  if (RestoreAfterFork)
  {
    SetThreads(FstThreadsAtFork);
  }
}


int avoid_openmp_hang_within_fork(DllInfo *dll)
{
  // Called once on loading package
#ifdef _OPENMP
  return pthread_atfork(&when_fork, &when_unfork, NULL);
#endif

  return 0;
}


SEXP hasopenmp()
{
  return Rf_ScalarLogical(HasOpenMP());
}

