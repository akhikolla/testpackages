/***************************************************************************
                             SRC/mixmod/Matrix/SymmetricMatrix.h  description
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
#ifndef XEMSYMMETRICMATRIX_H
#define XEMSYMMETRICMATRIX_H

#include "mixmod/Matrix/Matrix.h"
#include "mixmod/Utilities/maths/SelectLibrary.h"

namespace XEM {

/**
  @brief class GeneralMatrix
  @author F Langrognet & A Echenim
 */

class DiagMatrix;

class SymmetricMatrix : public Matrix {

public:

	/// Default constructor
	SymmetricMatrix();

	/// constructor : d*Id
	/// default value : Id
	SymmetricMatrix(int64_t pbDimension, double d = 1.0);

	SymmetricMatrix(SymmetricMatrix * A);

	/// Destructor
	virtual ~SymmetricMatrix();

	/// Compute determinant of symmetric matrix
	double determinant(Exception& errorType);

	/// Return store of symmetric matrix
	double * getStore();

	/// Return newmat symmetric matrix
	MATH::SymmetricMatrix * getValue();

	/// Return dimension of store
	int64_t getStoreDim();

	/// Inverse symmetric matrix
	void inverse(Matrix * & A);

	/// compute (x - mean)' this (x - mean) 
	double norme(double * xMoinsMean);

	/// compute : this =  A / d
	void equalToMatrixDividedByDouble(Matrix * A, double d);

	/// compute : this =  A * d
	void equalToMatrixMultiplyByDouble(Matrix*D, double d);

	/// add :  cik * xMoinsMean * xMoinsMean'  to this
	void add(double * xMoinsMean, double cik);

	// add : diag( cik * xMoinsMean * xMoinsMean' )  to this
	//void addDiag(double * xMoinsMean, double cik);

	/// Return store of a spherical matrix in a symmetric one
	double putSphericalValueInStore(double & store);
	/// Add store of a spherical matrix in a symmetric one
	double addSphericalValueInStore(double & store);

	double getSphericalStore();

	/// Return store of a diagonal matrix
	double* putDiagonalValueInStore(double * store);
	/// Add store of a diagonal matrix in a diagonal one
	double* addDiagonalValueInStore(double * store);

	double* getDiagonalStore();

	/// Return store of a diagonal matrix
	double* putSymmetricValueInStore(double * store);
	/// Add store of a diagonal matrix in a diagonal one
	double* addSymmetricValueInStore(double * store);

	double* getSymmetricStore();

	/// Return store of a diagonal matrix
	double* putGeneralValueInStore(double * store);
	/// Add store of a diagonal matrix in a diagonal one
	double* addGeneralValueInStore(double * store);

	double* getGeneralStore();

	/// compute general matrix SVD decomposition
	void computeSVD(DiagMatrix* & S, GeneralMatrix* & O);


	/// this =  (d * Identity)
	void operator=(const double& d);
	/// this =  this /  (Identity * d)
	void operator/=(const double& d);
	/// this =  this * ( Identity * d)
	void operator*=(const double& d);
	/// this =  this + matrix
	void operator+=(Matrix* M);
	/// this = matrix
	void operator=(Matrix* M);

	
	/// read symmetric matrix store in file
	void input(std::ifstream & fi);
	virtual void input(double ** variances);

	/* ///compute SVD decomposition for a symmetric matrix
	 void computeSVD(XEMDiagMatrix* & S, XEMGeneralMatrix* & O);*/

	/// compute Shape as diag(Ot . this . O ) / diviseur
	void computeShape_as__diag_Ot_this_O(DiagMatrix* & Shape, GeneralMatrix* & Ori, double diviseur = 1.0);
	//  double trace_this_O_Sm1_O(XEMGeneralMatrix* & O, XEMDiagMatrix* & S);

	/// compute this as : multi * (O * S * O' )
	void compute_as__multi_O_S_O(double multi, GeneralMatrix* & O, DiagMatrix* & S);

	/// compute this as O*S*O'
	void compute_as_O_S_O(GeneralMatrix* & O, double* & S_store);
	/// compute trace of this
	double computeTrace();

	/// compute this as M * M'
	void compute_as_M_tM(GeneralMatrix* M, int64_t d);

	/// compute this as matrix * vector
	void compute_as_M_V(SymmetricMatrix* M, double * V);

	/// compute this as double * matrix
	void compute_product_Lk_Wk(Matrix* Wk, double L);

	/// copute trace of W * C
	double compute_trace_W_C(Matrix * C);

	/// compute M as : M = ( O * S^{-1} * O' ) * this
	void compute_M_as__O_Sinverse_Ot_this(GeneralMatrix & M, GeneralMatrix* & O, DiagMatrix* & S);

	/// compute this as vector * vector'
	void compute_M_tM(double* V, int64_t l);

	/// gives : det(diag(this))
	double detDiag(Exception& errorType);

	/// trace( this * O * S^{-1} * O' )
	double trace_this_O_Sm1_O(GeneralMatrix* & O, DiagMatrix* & S);

	///set store

	void setSymmetricStore(double * store);
	void setGeneralStore(double * store);
	void setDiagonalStore(double * store);
	void setSphericalStore(double store);
	double** storeToArray() const;

protected:

	/// Symmetric matrix as in mathematical library (if specified)
	MATH::SymmetricMatrix * _value;

	/// store of matrix
	double * _store;
	
	/// dimension of store
	int64_t _s_storeDim;
};

// TODO static :
// int64_t XEMGeneralMatrix::_s_storeDim = 0;

inline double * SymmetricMatrix::getStore() {
	return _store;
}

inline MATH::SymmetricMatrix * SymmetricMatrix::getValue() {
	return _value;
}

inline int64_t SymmetricMatrix::getStoreDim() {
	return _s_storeDim;
}

inline void SymmetricMatrix::setSymmetricStore(double * store) {
	// _store = store;
	recopyTab(store, _store, _s_storeDim);
}

inline void SymmetricMatrix::setSphericalStore(double store) {
	THROW(OtherException, wrongMatrixType);
}

inline void SymmetricMatrix::setGeneralStore(double * store) {
	THROW(OtherException, wrongMatrixType);
}

inline void SymmetricMatrix::setDiagonalStore(double * store) {
	THROW(OtherException, wrongMatrixType);
}

/* TODO static
inline void XEMGeneralMatrix::initiate(){
  _s_storeDim = _s_pbDimension * _s_pbDimension / 2; 
} 
 */

}

#endif
