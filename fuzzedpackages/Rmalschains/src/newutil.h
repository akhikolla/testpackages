/**
 * Copyright 2008, Daniel Molina Cabrera <danimolina@gmail.com>
 * 
 * This file is part of software Realea
 * 
 * Realea is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Realea is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _NEWUTIL_H

#define _NEWUTIL_H 1

#include <queue>

using namespace std;

#define WANT_STREAM                  // include.h will get stream fns
#define WANT_MATH                    // include.h will get math fns
#define SETUP_C_SUBSCRIPTS

#define USING_DOUBLE

#include "newmat.h"
#include "newmatap.h"                // need matrix applications
#include "newmatio.h"                // need matrix output routines
#include "real.h"

#ifdef use_namespace
using namespace NEWMAT;              // access NEWMAT namespace
#endif

typedef Matrix MyMatrix;
typedef ReturnMatrix MyReturnMatrix;

MyReturnMatrix sqrt(const MyMatrix &mat); 

double sqrt(int x); 

double norm(const MyMatrix &mat); 

MyReturnMatrix eye(int N); 

Real min(ColumnVector &mat);
Real max(ColumnVector &mat);

double pow2(double x);
double sqrt_double(double x);
double log_double(double x); 

MyReturnMatrix log(const MyMatrix &mat); 
double log(int x); 
/**
 * Devuelve la media del vector columna
 */
Real mean(const ColumnVector &mat);

MyReturnMatrix pow(const ColumnVector &mat, double exp); 

double pow2_double(double x);

MyReturnMatrix pow2(ColumnVector &weights); 
MyReturnMatrix pow2(RowVector &weights); 

MyReturnMatrix log(const MyMatrix &mat); 

MyReturnMatrix pow2_m(const MyMatrix &mat); 

void set_sort_matrix(Matrix *mat); 

bool sort_index_matrix(int i, int j); 
/**
 * Devuelve la media de la diagonal
 */
Real mean_diag(const DiagonalMatrix &mat);

/**
 * Devuelve la mediana del vector columna
 */

Real median(const ColumnVector &mat); 
Real median(const RowVector &mat); 

/**
 * Permite aplicar una división elemento a elemento de do vectores columnas
 */
MyReturnMatrix DivVectors(ColumnVector &A, const ColumnVector &B);

/**
 * Computes the percentiles in vector perc from vector inar
 * returns vector with length(res)==length(perc)
 * idx: optional index-array indicating sorted order
 */
MyReturnMatrix myprctile(RowVector &inar, int *perc, int nperc);

void copyRow(queue<Real> num, RowVector &row); 

void copyColumn(DiagonalMatrix diag, ColumnVector &row); 
void copyFromColumn(ColumnVector &row, DiagonalMatrix &diag);

MyReturnMatrix DotVectors(ColumnVector &A, const ColumnVector &B);

MyReturnMatrix DivVectors(ColumnVector &A, const ColumnVector &B);

Real min_positive(queue<Real> num); 
Real min_positive(ColumnVector &num); 


void range(int min, int max, int *Rang); 

/**
 * Permite obtener una lista de columnas
 */
void getColumns(Matrix &m, int *cols, int num, Matrix &result);

void checkDiag(Matrix &C, DiagonalMatrix &Diag); 

void checkAxis(ColumnVector &xmean, Real ccov, Real cs, Real damps, int countiter, ColumnVector &sigma, Matrix &C, Matrix &BD);

/**
 * Allow to create a ColumnVector from a vector
 *
 * @param vector vector
 * @param N size
 * @param column output
 */
void copyToColumn(Real vector, unsigned N, ColumnVector *column); 
void copyToColumn(vector<tReal> array, ColumnVector *column); 

MyReturnMatrix DotDivide(const MyMatrix &elem1, const MyMatrix &elem2); 
#endif
