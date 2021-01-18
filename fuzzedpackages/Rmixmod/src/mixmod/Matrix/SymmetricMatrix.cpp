/***************************************************************************
                             SRC/mixmod/Matrix/SymmetricMatrix.cpp  description
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

#include "mixmod/Matrix/SymmetricMatrix.h"
#include "mixmod/Matrix/DiagMatrix.h"
#include "mixmod/Matrix/GeneralMatrix.h"

namespace XEM {

//------------
// Constructor
//------------
SymmetricMatrix::SymmetricMatrix() {
	_value = NULL;
	_store = NULL;
	THROW(OtherException, wrongConstructorType);
}

SymmetricMatrix::SymmetricMatrix(int64_t pbDimension, double d) : Matrix(pbDimension) {
	_value = new MATH::SymmetricMatrix(_s_pbDimension);
	_store = _value->Store();

	_s_storeDim = _s_pbDimension * (_s_pbDimension + 1) / 2;
	(*this) = d;
}

// copy constructor
SymmetricMatrix::SymmetricMatrix(SymmetricMatrix * A) : Matrix(A) {
	_value = new MATH::SymmetricMatrix(_s_pbDimension);
	_store = _value->Store();

	_s_storeDim = _s_pbDimension * (_s_pbDimension + 1) / 2;

	recopyTab(A->getStore(), _store, _s_storeDim);
}

//----------
//Destructor
//----------
SymmetricMatrix::~SymmetricMatrix() {
	if (_value) {
		delete _value;
	}
	_value = NULL;
}

double SymmetricMatrix::determinant(Exception& errorType) {
	double det = 0;
	try {
		det = _value->determinant(_store);
	}
	catch (...) {
		throw errorType;
	}
	if (fabs(det) < minDeterminantValue) {
		throw NumericException(dynamic_cast<NumericException&> (errorType));
	}
	return det;
}

void SymmetricMatrix::equalToMatrixMultiplyByDouble(Matrix* D, double d) {
	THROW(OtherException, nonImplementedMethod);
}

double* SymmetricMatrix::getDiagonalStore() {
	THROW(OtherException, wrongMatrixType);
}

double* SymmetricMatrix::getSymmetricStore() {
	return (_store);
}

double* SymmetricMatrix::getGeneralStore() {
	return (_store);
}

double SymmetricMatrix::getSphericalStore() {
	THROW(OtherException, wrongMatrixType);
}

void SymmetricMatrix::compute_M_tM(double* V, int64_t l) {
	int64_t indice1 = l - 1, indice2;
	int64_t indiceStoreGammak = _s_storeDim - 1;
	int64_t dim = l / _s_pbDimension;
	while (indice1 > 0) {
		for (int64_t j = 0; j < dim; j++) {
			_store[indiceStoreGammak] += V[(indice1 - j)] * V[(indice1 - j)];
		}
		indiceStoreGammak -= 1;
		indice2 = indice1 - dim;
		while (indice2 > 0) {
			for (int64_t j = 0; j < dim; j++) {
				_store[indiceStoreGammak] += V[(indice1 - j)] * V[(indice2 - j)];
			}
			indice2 -= dim;
			indiceStoreGammak -= 1;
		}
		indice1 -= dim;
	}
}

void SymmetricMatrix::compute_product_Lk_Wk(Matrix* Wk, double L) {
	double * Wk_store;
	Wk_store = Wk->getSymmetricStore();
	for (int64_t p = 0; p < _s_storeDim; p++) {
		_store[p] += Wk_store[p] / L;
	}
}

double SymmetricMatrix::compute_trace_W_C(Matrix * C) {
	double tabLambdak_k = 0.0;
	double termesHorsDiag;
	int64_t p, q, r;
	double * C_store = C->getSymmetricStore();
	termesHorsDiag = 0.0;
	for (p = 0, r = 0; p < _s_pbDimension; p++, r++) {
		for (q = 0; q < p; q++, r++) {
			termesHorsDiag += _store[r] * C_store[r];
		}
		tabLambdak_k += _store[r] * C_store[r];
	}
	tabLambdak_k += 2.0 * termesHorsDiag;
	return tabLambdak_k;
}

void SymmetricMatrix::inverse(Matrix * & Inv) {
	//cout<<"Inv Symm :  "<<Inv<<endl;
	if (Inv == NULL) {
		Inv = new SymmetricMatrix(_s_pbDimension);
	}

	MATH::SymmetricMatrix* value_Inv = _value->Inverse(_store);

	Inv->setSymmetricStore(value_Inv->Store());
	//cout<<"Inv Symm :  "<<Inv<<endl;
	delete value_Inv;
}

double SymmetricMatrix::norme(double * xMoinsMean) {
	int64_t p, q, r;
	double termesHorsDiag = 0.0;
	double termesDiag = 0.0;
	double xMoinsMean_p;

	for (p = 0, r = 0; p < _s_pbDimension; p++, r++) {
		xMoinsMean_p = xMoinsMean[p];
		for (q = 0; q < p; q++, r++) {
			termesHorsDiag += xMoinsMean_p * xMoinsMean[q] * _store[r];
		}
		termesDiag += xMoinsMean_p * xMoinsMean_p * _store[r];
	}
	termesDiag += 2.0 * termesHorsDiag;
	return termesDiag;
}

double SymmetricMatrix::putSphericalValueInStore(double & store) {
	int64_t p, r;
	int64_t increment = 2;
	store = 0.0;

	for (p = 0, r = 0; p < _s_pbDimension; p++) {
		store += _store[r];
		r += increment;
		increment++;
	}
	store /= _s_pbDimension;
	return (store);

}

double SymmetricMatrix::addSphericalValueInStore(double & store) {
	int64_t p, r;
	int64_t increment = 2;
	for (p = 0, r = 0; p < _s_pbDimension; p++) {
		store += _store[r];
		r += increment;
		increment++;
	}
	store /= _s_pbDimension;
	return (store);
}

double* SymmetricMatrix::putDiagonalValueInStore(double * store) {
	int64_t p, r;
	int64_t increment = 2;
	for (p = 0, r = 0; p < _s_pbDimension; p++) {
		store[p] = _store[r];
		r += increment;
		increment++;
	}
	return (store);
}

double* SymmetricMatrix::addDiagonalValueInStore(double * store) {
	int64_t p, r;
	int64_t increment = 2;
	for (p = 0, r = 0; p < _s_pbDimension; p++) {
		store[p] += _store[r];
		r += increment;
		increment++;
	}
	return (store);
}

double* SymmetricMatrix::putSymmetricValueInStore(double * store) {
	for (int64_t p = 0; p < _s_storeDim; p++) {
		store[p] = _store[p];
	}
	return (store);
}

double* SymmetricMatrix::addSymmetricValueInStore(double * store) {
	for (int64_t p = 0; p < _s_storeDim; p++) {
		store[p] += _store[p];
	}
	return (store);
}

double* SymmetricMatrix::putGeneralValueInStore(double * store) {
	THROW(OtherException, wrongMatrixType);
}

double* SymmetricMatrix::addGeneralValueInStore(double * store) {
	THROW(OtherException, wrongMatrixType);
}

// (this) will be A / d
void SymmetricMatrix::equalToMatrixDividedByDouble(Matrix * A, double d) {
	A->putSymmetricValueInStore(_store);

	int64_t p;
	for (p = 0; p < _s_storeDim; p++) {
		_store[p] /= d;
	}
}

// add :  cik * xMoinsMean * xMoinsMean'  to this
void SymmetricMatrix::add(double * xMoinsMean, double cik) {

	int64_t p, q, r;
	double xMoinsMean_p;

	for (p = 0, r = 0; p < _s_pbDimension; p++, r++) {
		xMoinsMean_p = xMoinsMean[p];
		for (q = 0; q < p; q++, r++) {
			_store[r] += cik * xMoinsMean_p * xMoinsMean[q];
		} // end for q
		_store[r] += cik * xMoinsMean_p * xMoinsMean_p;
	}//end for p
}

// add : diag( cik * xMoinsMean * xMoinsMean' )  to this
/*void SymmetricMatrix::addDiag(double * xMoinsMean, double cik){
  int64_t p,q,r;
  double xMoinsMean_p;

  for(p=0,r=0 ; p<_s_pbDimension ; p++,r++){
	xMoinsMean_p = xMoinsMean[p];
	for(q=0 ; q<p ; q++,r++) ;
	_store[r]  +=  cik * xMoinsMean_p * xMoinsMean_p;
  }  
}*/

// set the value of (d x Identity) to this  
void SymmetricMatrix::operator=(const double& d) {

	int64_t p, q, r;
	for (p = 0, r = 0; p < _s_pbDimension; p++, r++) {
		for (q = 0; q < p; q++, r++) {
			_store[r] = 0.0;
		}
		_store[r] = d;
	}
}

// divide each element by d
void SymmetricMatrix::operator/=(const double& d) {
	int64_t p;
	for (p = 0; p < _s_storeDim; p++) {
		_store[p] /= d;
	}
}

// multiply each element by d
void SymmetricMatrix::operator*=(const double& d) {
	int64_t p;
	for (p = 0; p < _s_storeDim; p++) {
		_store[p] *= d;
	}
}

void SymmetricMatrix::operator=(Matrix* M) {
	M->putSymmetricValueInStore(_store);
}

//add M to this
void SymmetricMatrix::operator+=(Matrix* M) {
	M->addSymmetricValueInStore(_store);
}

// compute Shape as diag(Ot . this . O ) / diviseur
void SymmetricMatrix::computeShape_as__diag_Ot_this_O(
		DiagMatrix* & Shape, GeneralMatrix* & Ori, double diviseur) 
{
	int64_t i_index, j_index;

	int64_t p, q, r, j;
	double * O_store = Ori->getStore();
	double * Shape_store = Shape->getStore();

	double termesDiag, termesHorsDiag;
	double tmp;

	for (j = 0; j < _s_pbDimension; j++) {
		// computation of the [j,j] term of the diagonal

		//-- reset
		termesDiag = 0.0;
		termesHorsDiag = 0.0;

		i_index = j;
		for (p = 0, r = 0; p < _s_pbDimension; p++, r++, i_index += _s_pbDimension) {
			j_index = j;
			for (q = 0; q < p; q++, r++, j_index += _s_pbDimension) {
				termesHorsDiag += O_store[i_index] * O_store[j_index] * _store[r];
			}
			tmp = O_store[i_index];
			termesDiag += tmp * tmp * _store[r];
		}

		termesDiag += 2.0 * termesHorsDiag;
		termesDiag /= diviseur;
		Shape_store[j] = termesDiag;
	}
}

// compute this as : multi * (O * S * O' )
void SymmetricMatrix::compute_as__multi_O_S_O(double multi, GeneralMatrix* & O, DiagMatrix* & S) {

	int64_t i_index = 0;
	int64_t j_index;
	int64_t p, q, r, l;

	double * O_store = O->getStore();
	double * S_store = S->getStore();
	double tmp;

	for (p = 0, r = 0; p < _s_pbDimension; p++, i_index += _s_pbDimension) {
		j_index = 0;
		for (q = 0; q <= p; q++, r++, j_index += _s_pbDimension) {
			// compute this[i,j] = \multi * sum_{l} ( O[i,l] * 0[j,l] * S|l] ) 
			tmp = 0.0;
			for (l = 0; l < _s_pbDimension; l++) {
				tmp += O_store[i_index + l] * O_store[j_index + l] * S_store[l];
			}
			tmp *= multi;
			_store[r] = tmp;
		}
	}
}

// compute this as : (O * S * O' )
void SymmetricMatrix::compute_as_O_S_O(GeneralMatrix* & O, double* & S_store) {

	int64_t i_index = 0;
	int64_t j_index;
	int64_t p, q, r, l;

	for (int64_t i = 0; i < _s_storeDim; i++) {
		_store[i] = 0;
	}

	double * O_store = O->getStore();
	double tmp;
	for (p = 0, r = 0; p < _s_pbDimension; p++, i_index += _s_pbDimension) {
		j_index = 0;
		for (q = 0; q <= p; q++, r++, j_index += _s_pbDimension) {
			tmp = 0.0;
			for (l = 0; l < _s_pbDimension; l++) {
				tmp += O_store[i_index + l] * O_store[j_index + l] * S_store[l];
			}
			_store[r] = tmp;
		}
	}
}

//compute trace of a symmetric matrix
double SymmetricMatrix::computeTrace() {
	int64_t i;
	int64_t indice = 0;
	double trace = 0.0;


	i = 0;
	while (indice < _s_storeDim) {
		trace += _store[indice];
		i++;
		indice += i + 1;
	}
	return trace;
}

void SymmetricMatrix::computeSVD(DiagMatrix* & S, GeneralMatrix* & O) {
	int64_t dim = O->getPbDimension();
	MATH::DiagonalMatrix * tabShape_k = new MATH::DiagonalMatrix(dim);
	MATH::Matrix * tabOrientation_k = new MATH::Matrix(dim, dim);
	_value->computeSVD(tabShape_k, tabOrientation_k, _store);

	double * storeS = S->getStore();
	double * storeO = O->getStore();

	double * storeTabShape_k = tabShape_k->Store();
	double * storeTabOrientation_k = (*tabOrientation_k).Store();

	recopyTab(storeTabShape_k, storeS, dim);
	recopyTab(storeTabOrientation_k, storeO, dim * dim);

	delete tabShape_k;
	delete tabOrientation_k;
}

void SymmetricMatrix::compute_as_M_tM(GeneralMatrix* M, int64_t d) {

	int64_t indiceStoreM1 = 0, indiceStoreM2;
	int64_t indice = 0;
	int64_t k1 = 0, k2;
	int64_t DimStoreM = _s_pbDimension*_s_pbDimension;
	double * storeM = M->getStore();

	for (int64_t i = 0; i < _s_storeDim; i++) {
		_store[i] = 0;
	}

	while (indiceStoreM1 < DimStoreM) {
		k2 = k1;
		indiceStoreM2 = indiceStoreM1;
		while (indiceStoreM2 < DimStoreM) {

			for (int64_t j = 0; j < d; j++) {
				// attention vecteur contenant la matrice triangulaire supérieure
				_store[indice] += storeM[(indiceStoreM1 + j)] * storeM[(indiceStoreM2 + j)];
			}
			indiceStoreM2 = (k2 + 1) * _s_pbDimension;
			k2 += 1;
			indice += 1;
		}
		indiceStoreM1 = (k1 + 1) * _s_pbDimension;
		k1 += 1;
	}
}

void SymmetricMatrix::compute_as_M_V(SymmetricMatrix* M, double * V) {

	for (int64_t i = 0; i < _s_pbDimension; i++) {
		_store[i] = 0;
	}
	int64_t indiceV = 0, k = 0, indice = 0;
	int64_t indiceM = 0;
	double* storeM = M->getStore();

	while (indice < _s_pbDimension) {
		for (int64_t i = 0; i < (_s_pbDimension - k); i++) {
			_store[indice] += V[indiceV + i] * storeM[indiceM + i];
		}

		for (int64_t j = 1; j < (_s_pbDimension - k); j++) {
			_store[indice + j] += V[indiceV] * storeM[indiceM + j];
		}
		indiceM += (_s_pbDimension - k);
		k += 1;
		indiceV += 1;
		indice += 1;
	}
}

// compute M as : M = ( O * S^{-1} * O' ) * this
void SymmetricMatrix::compute_M_as__O_Sinverse_Ot_this(
		GeneralMatrix & M, GeneralMatrix* & O, DiagMatrix* & S) 
{
	double * M_store = M.getStore();
	double * O_store = O->getStore();
	double * S_store = S->getStore();

	int64_t i, j, l, p, r;
	int64_t O1_index = 0;
	int64_t O2_index;
	int64_t r_decalage;
	double tmp, omega;
	int64_t fillindex = 0;

	for (i = 0; i < _s_pbDimension; i++) {
		for (j = 0; j < _s_pbDimension; j++, fillindex++) {
			// filling tabMtmpk_store[i,j]

			tmp = 0.0;
			r_decalage = _s_pbDimension - j + 1;

			r = j;
			p = 0;
			O2_index = 0;


			while (p < j) {
				omega = 0.0;
				for (l = 0; l < _s_pbDimension; l++) {

					omega += O_store[O1_index + l] * O_store[O2_index + l] / S_store[l];
				}
				tmp += omega * _store[r];
				r += r_decalage;
				r_decalage--;
				p++;
				O2_index += _s_pbDimension;
			}

			while (p < _s_pbDimension) {
				omega = 0.0;

				for (l = 0; l < _s_pbDimension; l++) {
					omega += O_store[O1_index + l] * O_store[O2_index + l] / S_store[l];
				}
				tmp += omega * _store[r];
				r++;
				O2_index += _s_pbDimension;
				p++;
			}


			M_store[fillindex] = tmp;
		}
		O1_index += _s_pbDimension;
	}
}

void SymmetricMatrix::input(std::ifstream & fi) {
	int64_t i, j, r = 0;
	double garbage;

	for (i = 0; i < _s_pbDimension; i++) {
    for (j = 0; j < i + 1; j++) {
			_store[r] = getDoubleFromStream(fi);
      r++;
    }
    for (j = i + 1; j < _s_pbDimension; j++) {
			// we don't need values (all are 0 ?)
			getDoubleFromStream(fi);
		}
	}
}

void SymmetricMatrix::input(double ** variances) {
	int64_t i, j, r = 0;

	for (i = 0; i < _s_pbDimension; i++) {
		for (j = 0; j < i + 1; j++) {
			_store[r] = variances[i][j];
			r++;
		}
		for (j = i + 1; j < _s_pbDimension; j++) {
		}
	}
}

// gives : det(diag(this))
double SymmetricMatrix::detDiag(Exception& errorType) {
	int64_t p, q, r;
	double det = 1.0;

	for (p = 0, r = 0; p < _s_pbDimension; p++, r++) {
		for (q = 0; q < p; q++, r++);
		det *= _store[r];
	}
	if (det < minDeterminantValue)
		throw errorType;
	return det;
}

// trace( this * O * S^{-1} * O' )
double SymmetricMatrix::trace_this_O_Sm1_O(GeneralMatrix* & O, DiagMatrix* & S) {
	double * O_store = O->getStore();
	double * S_store = S->getStore();

	double trace = 0.0;
	double termesHorsDiag = 0.0;
	double tmp, tmp2;

	int64_t i_index = 0;
	int64_t j_index;
	int64_t p, q, r, l;

	for (p = 0, r = 0; p < _s_pbDimension; p++, r++, i_index += _s_pbDimension) {
		j_index = 0;
		for (q = 0; q < p; q++, r++, j_index += _s_pbDimension) {
			tmp = 0.0;

			for (l = 0; l < _s_pbDimension; l++) {
				tmp += O_store[i_index + l] * O_store[j_index + l] / S_store[l];
			}
			tmp *= _store[r];
			termesHorsDiag += tmp;
		}
		tmp = 0.0;
		for (l = 0; l < _s_pbDimension; l++) {
			tmp2 = O_store[i_index + l];
			tmp += tmp2 * tmp2 / S_store[l];
		}
		tmp *= _store[r];
		trace += tmp;
	}
	trace += 2.0 * termesHorsDiag;

	return trace;
}

double** SymmetricMatrix::storeToArray() const {

	int64_t i, j, k = (_s_storeDim - 1);
	double** newStore = new double*[_s_pbDimension];

	for (i = 0; i < _s_pbDimension; ++i) {
		newStore[i] = new double[_s_pbDimension];
	}
	for (i = (_s_pbDimension - 1); i>-1; --i) {
		newStore[i][i] = _store[k];
		k--;
		for (j = (i - 1); j>-1; --j) {
			newStore[i][j] = _store[k];
			newStore[j][i] = _store[k];
			k--;
		}

	}

	return newStore;
}

}
