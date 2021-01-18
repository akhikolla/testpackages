/***************************************************************************
                             SRC/mixmod/Matrix/GeneralMatrix.h  description
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
#ifndef XEMGENERALMATRIX_H
#define XEMGENERALMATRIX_H

#include "mixmod/Matrix/Matrix.h"
#include "mixmod/Utilities/maths/SelectLibrary.h"

namespace XEM {

// pre-declaration
class DiagMatrix;

/**
  @brief class GeneralMatrix
  @author F Langrognet & A Echenim
 */

class GeneralMatrix : public Matrix {

public:

	/// Default constructor
	GeneralMatrix();

	/// constructor : d*Id
	///default value = Id
	GeneralMatrix(int64_t pbDimension, double d = 1.0);

	GeneralMatrix(GeneralMatrix * A);

	/// Destructor
	virtual ~GeneralMatrix();

	/// compute determinant of general matrix
	double determinant(Exception& errorType);

	/// return store of general matrix
	double * getStore();

	/// return newmat general matrix
	MATH::Matrix * getValue();

	/// return dimension of store
	int64_t getStoreDim();

	/// inverse general matrix
	void inverse(Matrix * & A);

	void compute_product_Lk_Wk(Matrix* Wk, double L);

	/// compute (x - mean)' this (x - mean) 
	double norme(double * xMoinsMean);

	/// this =  A / d
	void equalToMatrixDividedByDouble(Matrix * A, double d);
	/// this =   A * d
	void equalToMatrixMultiplyByDouble(Matrix*D, double d);


	/// add :  cik * xMoinsMean * xMoinsMean'  to this
	void add(double * xMoinsMean, double cik);

	// add : diag( cik * xMoinsMean * xMoinsMean' )  to this
	//void addDiag(double * xMoinsMean, double cik);

	/// return store of a spherical matrix in a general one
	double putSphericalValueInStore(double & store);
	/// add store of a spherical matrix in a general one
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

	/// this =  (d x Identity)
	void operator=(const double& d);

	/// this =  this / (d x Identity)
	void operator/=(const double& d);

	/// this =  this * (d x Identity)
	void operator*=(const double& d);

	/// this =  this + matrix
	void operator+=(Matrix* M);

	/// this =  matrix
	void operator=(Matrix* M);


	/// edit general matrix
	void edit(std::ostream& flux, std::string before, std::string sep, int64_t dim);

	/// read general matrix from input file
	void input(std::ifstream & fi);
	/// read general matrix from input file
	void input(std::ifstream & fi, int64_t dim);
	virtual void input(double ** variances);

	
	/// compute general matrix SVD decomposition
	void computeSVD(DiagMatrix* & S, GeneralMatrix* & O);

	/// compute Shape as diag(Ot . this . O ) / diviseur
	void computeShape_as__diag_Ot_this_O(DiagMatrix* & Shape, GeneralMatrix* & Ori, double diviseur = 1.0);

	/// compute this as : multi * (O * S * O' )
	void compute_as__multi_O_S_O(double multi, GeneralMatrix* & O, DiagMatrix *& S);

	/// compute this as O * S *O'
	void compute_as_O_S_O(GeneralMatrix* & O, double* & S_store);

	/// compute trace of this
	double computeTrace();

	/// compute this as matrix * matrix'
	void compute_as_M_tM(GeneralMatrix* M, int64_t d);

	/// compute this as matrix * vector
	void compute_as_M_V(GeneralMatrix* M, double * V);
	/// compute this as vector multiplied by matrix
	void multiply(double * V, int64_t nk, GeneralMatrix * Q);

	/// compute M as : M = ( O * S^{-1} * O' ) * this
	void compute_M_as__O_Sinverse_Ot_this(GeneralMatrix & M, GeneralMatrix* & O, DiagMatrix* & S);
	double compute_trace_W_C(Matrix * C);
	//  void computeShape_as__diag_Ot_this_O(XEMDiagMatrix* & Shape, XEMGeneralMatrix* & Ori, double diviseur = 1.0);
	/// gives : det(diag(this))
	double detDiag(Exception& errorType);

	/// trace( this * O * S^{-1} * O' )
	double trace_this_O_Sm1_O(GeneralMatrix* & O, DiagMatrix* & S);

	//void refreshStore();  

	void setSymmetricStore(double * store);
	void setGeneralStore(double * store);
	void setDiagonalStore(double * store);
	void setSphericalStore(double store);

	double** storeToArray() const;

protected:

	// General matrix as in mathematical library
	MATH::Matrix * _value;

	double * _store;

	int64_t _s_storeDim;
};

// TODO static :
// int64_t XEMGeneralMatrix::_s_storeDim = 0;

inline double * GeneralMatrix::getStore() {
	return _store;
}

inline MATH::Matrix * GeneralMatrix::getValue() {
	return _value;
}

inline int64_t GeneralMatrix::getStoreDim() {
	return _s_storeDim;
}

inline void GeneralMatrix::setSymmetricStore(double * store) {
	THROW(OtherException, wrongMatrixType);
}

inline void GeneralMatrix::setSphericalStore(double store) {
	THROW(OtherException, wrongMatrixType);
}

inline void GeneralMatrix::setGeneralStore(double * store) {
	//_store = store;
	recopyTab(store, _store, _s_storeDim);
}

inline void GeneralMatrix::setDiagonalStore(double * store) {
	THROW(OtherException, wrongMatrixType);
}

/* TODO static
inline void GeneralMatrix::initiate(){
  _s_storeDim = _s_pbDimension * _s_pbDimension / 2; 
} 
 */

}

#endif
