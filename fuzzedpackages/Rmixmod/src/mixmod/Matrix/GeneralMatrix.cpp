/***************************************************************************
                             SRC/mixmod/Matrix/GeneralMatrix.cpp  description
    copyright            : (C) MIXMOD Team - 2001-2016
    email                : contact@mixmod.org
 ***************************************************************************/

/***************************************************************************
    This file is part of MIXMOD
    
    MIXMOD is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MIXMOD is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MIXMOD.  If not, see <http://www.gnu.org/licenses/>.

    All informations available on : http://www.mixmod.org                                                                                               
***************************************************************************/

#include "mixmod/Matrix/GeneralMatrix.h"
#include "mixmod/Matrix/DiagMatrix.h"

namespace XEM {

//------------
// Constructor
//------------
GeneralMatrix::GeneralMatrix() {
	_value = NULL;
	_store = NULL;
	THROW(OtherException, wrongConstructorType);
}

GeneralMatrix::GeneralMatrix(int64_t pbDimension, double d) : Matrix(pbDimension) {
	_value = new MATH::Matrix(pbDimension, pbDimension);
	_store = _value->Store();

	_s_storeDim = pbDimension * pbDimension;
	(*this) = 1.0;
}

// copy constructor
GeneralMatrix::GeneralMatrix(GeneralMatrix * A) : Matrix(A) {
	_value = new MATH::Matrix(_s_pbDimension, _s_pbDimension);
	_store = _value->Store();

	_s_storeDim = _s_pbDimension * _s_pbDimension;

	recopyTab(A->getStore(), _store, _s_storeDim);
}

//----------
//Destructor
//----------
GeneralMatrix::~GeneralMatrix() {
	if (_value) {
		delete _value;
	}
	_value = NULL;
}

double GeneralMatrix::determinant(Exception& errorType) {
	THROW(OtherException, nonImplementedMethod);
}

void GeneralMatrix::equalToMatrixMultiplyByDouble(Matrix* D, double d) {
	THROW(OtherException, nonImplementedMethod);
}

void GeneralMatrix::multiply(double * V, int64_t nk, GeneralMatrix * Q) {
	int64_t indiceV, indiceQ, indice = 0;
	double * storeQ = Q->getStore();

	for (indiceV = 0; indiceV < _s_pbDimension; indiceV++) {
		for (indiceQ = 0; indiceQ < nk; indiceQ++) {
			_store[indice] = V[indiceV] * storeQ[indiceQ] 
					+ V[indiceV + _s_pbDimension] * storeQ[indiceQ + _s_pbDimension];
			indice++;
		}
	}
}

double* GeneralMatrix::getDiagonalStore() {
	THROW(OtherException, wrongMatrixType);
}

double* GeneralMatrix::getSymmetricStore() {
	THROW(OtherException, wrongMatrixType);
}

double* GeneralMatrix::getGeneralStore() {
	return (_store);
}

double GeneralMatrix::getSphericalStore() {
	THROW(OtherException, wrongMatrixType);
}

double GeneralMatrix::compute_trace_W_C(Matrix * C) {
	THROW(OtherException, nonImplementedMethod);
}

void GeneralMatrix::computeShape_as__diag_Ot_this_O(
		DiagMatrix* & Shape, GeneralMatrix* & Ori, double diviseur) 
{
	THROW(OtherException, nonImplementedMethod);
}

void GeneralMatrix::inverse(Matrix * & A) {
	THROW(OtherException, nonImplementedMethod);
}

void GeneralMatrix::compute_product_Lk_Wk(Matrix* Wk, double L) {
	THROW(OtherException, nonImplementedMethod);
}

double GeneralMatrix::norme(double * xMoinsMean) {
	THROW(OtherException, nonImplementedMethod);
}

// (this) will be A / d
void GeneralMatrix::equalToMatrixDividedByDouble(Matrix * A, double d) {
	THROW(OtherException, nonImplementedMethod);
}

// add :  cik * xMoinsMean * xMoinsMean'  to this
void GeneralMatrix::add(double * xMoinsMean, double cik) {
	THROW(OtherException, nonImplementedMethod);
}

// add : diag( cik * xMoinsMean * xMoinsMean' )  to this
/*void GeneralMatrix::addDiag(double * xMoinsMean, double cik){
  THROW(InputException,nonImplementedMethod);
}*/

double GeneralMatrix::putSphericalValueInStore(double & store) {
	THROW(OtherException, wrongMatrixType);
}

double GeneralMatrix::addSphericalValueInStore(double & store) {
	THROW(OtherException, wrongMatrixType);
}

double * GeneralMatrix::putDiagonalValueInStore(double * store) {
	THROW(OtherException, wrongMatrixType);
}

double * GeneralMatrix::addDiagonalValueInStore(double * store) {
	THROW(OtherException, wrongMatrixType);
}

double * GeneralMatrix::putSymmetricValueInStore(double * store) {
	THROW(OtherException, wrongMatrixType);
}

double * GeneralMatrix::addSymmetricValueInStore(double * store) {
	THROW(OtherException, wrongMatrixType);
}

double* GeneralMatrix::putGeneralValueInStore(double * store) {
	for (int64_t p = 0; p < _s_storeDim; p++) {
		store[p] = _store[p];
	}
	return (store);
}

double* GeneralMatrix::addGeneralValueInStore(double * store) {
	for (int64_t p = 0; p < _s_storeDim; p++) {
		store[p] += _store[p];
	}
	return (store);
}

// set the value of (d x Identity) to this  
void GeneralMatrix::operator=(const double& d) {
	//THROW(InputException,nonImplementedMethod);
	int64_t indice = 0;
	while (indice < _s_storeDim) {
		for (int64_t i = 0; i < _s_pbDimension; i++) {
			for (int64_t j = 0; j < _s_pbDimension; j++) {
				if (i == j) {
					_store[indice] = d;
				}
				else {
					_store[indice] = 0.0;
				}
				indice++;
			}
		}
	}
}

// divide each element by d
void GeneralMatrix::operator/=(const double& d) {
	THROW(OtherException, nonImplementedMethod);
}

// multiply each element by d
void GeneralMatrix::operator*=(const double& d) {
	THROW(OtherException, nonImplementedMethod);
}

void GeneralMatrix::operator=(Matrix * M) {
	M->putGeneralValueInStore(_store);
}

//add M to this
void GeneralMatrix::operator+=(Matrix* M) {
	M->addGeneralValueInStore(_store);
}

void GeneralMatrix::edit(std::ostream& flux, std::string before, std::string sep, int64_t dim) {
	if (dim <= 0) {
		//ugly fix [bauder]: TODO = understand why dim can be 0 
		//(it shouldn't... but it does on example da2), and fix in upper classes.
		return;
	}
	for (int64_t p = 0; p < _s_pbDimension; p++) {
		flux << before;
		// TODO: << operator for double*
		double* row = _value->GetRow(p);
		for (int64_t i=0; i<dim-1; i++) 
			flux << row[i] << ",";
		flux << row[dim-1];
		flux << sep;
	}
}

/// compute this as : multi * (O * S * O' )
void GeneralMatrix::compute_as__multi_O_S_O(double multi, GeneralMatrix* & O, DiagMatrix* & S) {
	THROW(OtherException, nonImplementedMethod);
}

/// compute this as : (O * S * O' )
//void GeneralMatrix::compute_as_O_S_O(Matrix & O, double* & S_store){
void GeneralMatrix::compute_as_O_S_O(GeneralMatrix* & O, double* & S_store) {
	THROW(OtherException, nonImplementedMethod);
}

//compute trace of a symmetric matrix
double GeneralMatrix::computeTrace() {
	THROW(OtherException, nonImplementedMethod);
}

void GeneralMatrix::computeSVD(DiagMatrix* & S, GeneralMatrix* & O) {
	THROW(OtherException, nonImplementedMethod);
}

void GeneralMatrix::compute_as_M_tM(GeneralMatrix* M, int64_t d) {
	THROW(OtherException, nonImplementedMethod);
}

void GeneralMatrix::compute_as_M_V(GeneralMatrix* M, double * V) {
	THROW(OtherException, nonImplementedMethod);
}

// compute M as : M = ( O * S^{-1} * O' ) * this
void GeneralMatrix::compute_M_as__O_Sinverse_Ot_this(
		GeneralMatrix & M, GeneralMatrix* & O, DiagMatrix* & S) 
{
	THROW(OtherException, nonImplementedMethod);
}

void GeneralMatrix::input(std::ifstream & fi) {
	int64_t i, j, r = 0;

	for (i = 0; i < _s_pbDimension; i++) {
		for (j = 0; j < _s_pbDimension; j++, r++)
			_store[r] = getDoubleFromStream(fi);
	}
}

void GeneralMatrix::input(double ** variances) {
	int64_t i, j, r = 0;

	for (i = 0; i < _s_pbDimension; i++) {
		for (j = 0; j < _s_pbDimension; j++, r++) {
			_store[r] = variances[i][j];
		}
	}
}

void GeneralMatrix::input(std::ifstream & fi, int64_t dim) {
	int64_t i, j, r = 0;

	for (i = 0; i < _s_pbDimension; i++) {
		for (j = 0; j < dim; j++, r++)
			_store[r] = getDoubleFromStream(fi);
		for (j = dim; j < _s_pbDimension; j++, r++) {
			_store[r] = 0.0;
		}
	}
}

// gives : det(diag(this))
double GeneralMatrix::detDiag(Exception& errorType) {
	THROW(OtherException, nonImplementedMethod);
}

// trace( this * O * S^{-1} * O' )
double GeneralMatrix::trace_this_O_Sm1_O(GeneralMatrix* & O, DiagMatrix* & S) {
	THROW(OtherException, nonImplementedMethod);
}

double** GeneralMatrix::storeToArray() const {

	int64_t i, j, k = 0;
	double** newStore = new double*[_s_pbDimension];
	for (i = 0; i < _s_pbDimension; ++i) {
		newStore[i] = new double[_s_pbDimension];
		for (j = 0; j < _s_pbDimension; ++j) {
			newStore[i][j] = _store[k];
			k++;
		}
	}

	return newStore;
}

}
