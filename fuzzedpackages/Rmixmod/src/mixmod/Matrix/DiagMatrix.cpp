/***************************************************************************
                             SRC/mixmod/Matrix/DiagMatrix.cpp  description
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
#include "mixmod/Matrix/DiagMatrix.h"
#include "mixmod/Matrix/GeneralMatrix.h"
#include "mixmod/Matrix/SymmetricMatrix.h"

namespace XEM {

//------------
// Constructor
//------------
DiagMatrix::DiagMatrix() {
	_store = NULL;
	THROW(OtherException, wrongConstructorType);
}

DiagMatrix::DiagMatrix(int64_t pbDimension, double d) : Matrix(pbDimension) {
	_store = new double[_s_pbDimension];
	for (int64_t i = 0; i < _s_pbDimension; i++) {
		_store[i] = d;
	}
}

DiagMatrix::DiagMatrix(DiagMatrix * A) : Matrix(A) {
	_store = copyTab(A->getStore(), _s_pbDimension);
}

//----------
//Destructor
//----------
DiagMatrix::~DiagMatrix() {
	if (_store) {
		delete[] _store;
	}
	_store = NULL;
}

double DiagMatrix::determinant(Exception& errorType) {
	int64_t p;
	double det = _store[0];
	for (p = 1; p < _s_pbDimension; p++) {
		det *= _store[p];
	}

	if (det < minDeterminantValue)
		throw NumericException(dynamic_cast<NumericException&> (errorType));

	return det;
}

void DiagMatrix::compute_product_Lk_Wk(Matrix* Wk, double L) {
	THROW(OtherException, nonImplementedMethod);
}

void DiagMatrix::inverse(Matrix * & Inv) {
	//cout<<"Inv diag :  "<<Inv<<endl;
	if (Inv == NULL) {
		Inv = new DiagMatrix(_s_pbDimension);
	}
	double * Inv_store = new double[_s_pbDimension];
	int64_t p;
	for (p = 0; p < _s_pbDimension; p++) {
		Inv_store[p] = 1.0 / _store[p];
	}

	Inv->setDiagonalStore(Inv_store);

	delete [] Inv_store;
}

double* DiagMatrix::getDiagonalStore() {
	return (_store);
}

double* DiagMatrix::getSymmetricStore() {
	THROW(OtherException, wrongMatrixType);
}

double* DiagMatrix::getGeneralStore() {
	THROW(OtherException, wrongMatrixType);
}

double DiagMatrix::getSphericalStore() {
	THROW(OtherException, wrongMatrixType);
}

double DiagMatrix::norme(double * xMoinsMean) {
	int64_t p;
	double termesDiag = 0.0;
	double xMoinsMean_p;

	for (p = 0; p < _s_pbDimension; p++) {
		xMoinsMean_p = xMoinsMean[p];
		termesDiag += xMoinsMean_p * xMoinsMean_p * _store[p];
	}
	return termesDiag;
}

void DiagMatrix::compute_as__multi_O_S_O(double multi, GeneralMatrix* & O, DiagMatrix *& S) {
	THROW(OtherException, nonImplementedMethod);
}

double DiagMatrix::trace_this_O_Sm1_O(GeneralMatrix* & O, DiagMatrix* & S) {
	THROW(OtherException, nonImplementedMethod);
}

double DiagMatrix::compute_trace_W_C(Matrix * C) {
	THROW(OtherException, nonImplementedMethod);
}

void DiagMatrix::computeShape_as__diag_Ot_this_O(
		DiagMatrix* & Shape, GeneralMatrix* & Ori, double diviseur) 
{
	THROW(OtherException, nonImplementedMethod);
}

double DiagMatrix::putSphericalValueInStore(double & store) {
	store = 0.0;
	int64_t p;
	for (p = 0; p < _s_pbDimension; p++) {
		store += _store[p];
	}
	store /= _s_pbDimension;
	return (store);
}

double DiagMatrix::addSphericalValueInStore(double & store) {
	int64_t p;
	for (p = 0; p < _s_pbDimension; p++) {
		store += _store[p];
	}
	store /= _s_pbDimension;
	return (store);
}

double* DiagMatrix::putDiagonalValueInStore(double * store) {
	for (int64_t p = 0; p < _s_pbDimension; p++) {
		store[p] = _store[p];
	}
	return (store);
}

double* DiagMatrix::addDiagonalValueInStore(double * store) {
	for (int64_t p = 0; p < _s_pbDimension; p++) {
		store[p] += _store[p];
	}
	return (store);
}

double* DiagMatrix::addSymmetricValueInStore(double * store) {
	// return the store of of a symmetric matrix with this on the diag
	//int64_t dimStore = _s_pbDimension*(_s_pbDimension+1)/2;
	// double * store = new double[dimStore];

	int64_t p, q, r;
	for (p = 0, r = 0; p < _s_pbDimension; p++, r++) {
		for (q = 0; q < p; q++, r++) {
			store[r] = 0.0;
		}
		store[r] += _store[p];
	}
	return (store);
}

double* DiagMatrix::putSymmetricValueInStore(double * store) {
	// return the store of of a symmetric matrix with this on the diag
	// int64_t dimStore = _s_pbDimension*(_s_pbDimension+1)/2;
	// double * store = new double[dimStore];

	int64_t p, q, r;
	for (p = 0, r = 0; p < _s_pbDimension; p++, r++) {
		for (q = 0; q < p; q++, r++) {
			store[r] = 0.0;
		}
		store[r] = _store[p];
	}
	return (store);
}

double* DiagMatrix::putGeneralValueInStore(double * store) {
	THROW(OtherException, wrongMatrixType);
}

double* DiagMatrix::addGeneralValueInStore(double * store) {
	THROW(OtherException, wrongMatrixType);
}

void DiagMatrix::computeSVD(DiagMatrix* & S, GeneralMatrix* & O) {
	THROW(OtherException, nonImplementedMethod);
}

double DiagMatrix::computeTrace() {

	double trace = 0.0;
	for (int64_t i = 0; i < _s_pbDimension; i++) {
		trace += _store[i];
	}
	return trace;
}

// (this) will be A / d
void DiagMatrix::equalToMatrixDividedByDouble(Matrix * A, double d) {
	A->putDiagonalValueInStore(_store);

	int64_t p;
	for (p = 0; p < _s_pbDimension; p++) {
		_store[p] /= d;
	}
}

void DiagMatrix::equalToMatrixMultiplyByDouble(Matrix* D, double d) {
	D->putDiagonalValueInStore(_store);
	int64_t p;
	for (p = 0; p < _s_pbDimension; p++) {
		_store[p] *= d;
	}
}

// add :  cik * xMoinsMean * xMoinsMean'  to this
void DiagMatrix::add(double * xMoinsMean, double cik) {

	int64_t p;
	double xMoinsMean_p;

	for (p = 0; p < _s_pbDimension; p++) {
		xMoinsMean_p = xMoinsMean[p];
		_store[p] += cik * xMoinsMean_p * xMoinsMean_p;
	}//end for p
}

// add : diag( cik * xMoinsMean * xMoinsMean' )  to this
/*void DiagMatrix::addDiag(double * xMoinsMean, double cik){
  
  int64_t p;
  double xMoinsMean_p;

  for(p=0; p<_s_pbDimension ; p++){
	xMoinsMean_p = xMoinsMean[p];
	_store[p]   +=  cik * xMoinsMean_p * xMoinsMean_p;
  }//end for p
  
}*/

// set the value of (d x Identity) to this  
void DiagMatrix::operator=(const double& d) {
	int64_t p;

	for (p = 0; p < _s_pbDimension; p++) {
		_store[p] = d;
	}
}

// divide each element by d
void DiagMatrix::operator/=(const double& d) {
	int64_t p;
	for (p = 0; p < _s_pbDimension; p++) {
		_store[p] /= d;
	}
}

// multiply each element by d
void DiagMatrix::operator*=(const double& d) {
	int64_t p;
	for (p = 0; p < _s_pbDimension; p++) {
		_store[p] *= d;
	}
}

//add M to this
void DiagMatrix::operator+=(Matrix* M) {
	M -> addDiagonalValueInStore(_store);
}

void DiagMatrix::operator=(Matrix* M) {
	M -> putDiagonalValueInStore(_store);
}

void DiagMatrix::input(std::ifstream & fi) {
	int64_t p, q;
	double garbage;

	for (p = 0; p < _s_pbDimension; p++) {
		// useless because all are 0
    for (q = 0; q < p; q++)
			getDoubleFromStream(fi);

    // here i==j so we are in the diagonal
		_store[p] = getDoubleFromStream(fi);

    // useless because all are 0
    for (q = p + 1; q < _s_pbDimension; q++)
			getDoubleFromStream(fi);
  }
}

void DiagMatrix::input(double ** variances) {
	int64_t p, q;
	for (p = 0; p < _s_pbDimension; p++) {
		// useless because all are 0
		for (q = 0; q < p; q++) {
		}

		// here i==j so we are in the diagonal
		_store[p] = variances[p][q];

		// useless because all are 0
		for (q = p + 1; q < _s_pbDimension; q++) {
		}
	}
}

double DiagMatrix::detDiag(Exception& errorType) {
	return determinant(errorType);
}

void DiagMatrix::sortDiagMatrix() {
	int64_t max;
	for (int64_t i = 0; i < _s_pbDimension; i++) {
		max = i;
		for (int64_t j = i + 1; j < _s_pbDimension; j++) {
			if (_store[j] > _store[max]) {
				max = j;
			}
		}
		if (max != i) { // swich
			double tmp = _store[i];
			_store[i] = _store[max];
			_store[max] = tmp;
		}
	}

	/*for (int64_t i=1; i<= _s_pbDimension;i++){
			//search the max eigenvalue
			max = i;
	 for (int64_t j=i; j<_s_pbDimension; j++){
				if (_store[j-1] > _store[max-1])
						max = j;
				if (max != i){
					// switch
					double tmp = _store[max-1];
					_store[max-1] = _store[i-1];
					_store[i-1] = tmp;
				}
			}
		}*/
}

double** DiagMatrix::storeToArray() const {

	int64_t i, j;
	double** newStore = new double*[_s_pbDimension];
	for (i = 0; i < _s_pbDimension; ++i) {
		newStore[i] = new double[_s_pbDimension];
	}
	for (i = 0; i < _s_pbDimension; ++i) {

		for (j = 0; j < _s_pbDimension; ++j) {
			if (i == j) {
				newStore[i][j] = _store[i];
			}
			else {
				newStore[i][j] = 0;
			}
		}
	}

	return newStore;
}

}
