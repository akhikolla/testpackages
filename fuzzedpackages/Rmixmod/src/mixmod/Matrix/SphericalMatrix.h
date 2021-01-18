/***************************************************************************
                             SRC/mixmod/Matrix/SphericalMatrix.h  description
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
#ifndef XEMSPHERICALMATRIX_H
#define XEMSPHERICALMATRIX_H

#include "mixmod/Matrix/Matrix.h"

namespace XEM {

// pre-declaration
class GeneralMatrix;
class DiagMatrix;

/**
  @brief class XEMSphericalMatrix
  @author F Langrognet & A Echenim
 */

class SphericalMatrix : public Matrix {

public:

	/// Default constructor
	SphericalMatrix();

	/// contructor
	/// default initialisation : Id
	SphericalMatrix(int64_t pbDimension, double initValue = 1.0);

	SphericalMatrix(SphericalMatrix * A);

	/// Desctructor
	virtual ~SphericalMatrix();

	/// compute determinant of spherical matrix
	double determinant(Exception& errorType);

	/// return store of spherical matrix
	double getStore();


	/// inverse spherical matrix
	void inverse(Matrix * & A);

	void compute_product_Lk_Wk(Matrix* Wk, double L);

	/// add a to the value of this
	void addToValue(double a);

	/// compute (x - mean)' this (x - mean) 
	double norme(double * xMoinsMean);

	/// (this) will be A / d
	void equalToMatrixDividedByDouble(Matrix * A, double d);

	/// (this) will be A * d
	void equalToMatrixMultiplyByDouble(Matrix*D, double d);

	/// compute trace of spherical matrix
	double computeTrace();


	/// add :  cik * xMoinsMean * xMoinsMean'  to this
	void add(double * xMoinsMean, double cik);

	/// add : diag( cik * xMoinsMean * xMoinsMean' )  to this
	void addDiag(double * xMoinsMean, double cik);

	/// this  = d * Identity
	void operator=(const double& d);

	/// this = this / (d * Identity)
	void operator/=(const double& d);

	/// this = this *  (d * Identity)
	void operator*=(const double& d);
	/// this = this + matrix
	void operator+=(Matrix* M);
	/// this = matrix
	void operator=(Matrix* M);


	/// read spherical matrix store in file
	void input(std::ifstream & fi);
	virtual void input(double ** variances);

	/// return store of a spherical matrix
	double putSphericalValueInStore(double & store);
	/// add store of a spherical matrix
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

	void compute_as__multi_O_S_O(double multi, GeneralMatrix* & O, DiagMatrix *& S);
	double trace_this_O_Sm1_O(GeneralMatrix* & O, DiagMatrix* & S);
	double compute_trace_W_C(Matrix * C);
	void computeShape_as__diag_Ot_this_O(DiagMatrix* & Shape, GeneralMatrix* & Ori, double diviseur = 1.0);
	/// gives : det(diag(this))
	double detDiag(Exception& errorType);

	void setSymmetricStore(double * store);
	void setGeneralStore(double * store);
	void setDiagonalStore(double * store);
	void setSphericalStore(double store);
	double** storeToArray() const;

protected:

	double _store;
};

inline double SphericalMatrix::getStore() {
	return _store;
}

inline void SphericalMatrix::setSymmetricStore(double * store) {
	THROW(OtherException, wrongMatrixType);
}

inline void SphericalMatrix::setSphericalStore(double store) {
	_store = store;
}

inline void SphericalMatrix::setGeneralStore(double * store) {
	THROW(OtherException, wrongMatrixType);
}

inline void SphericalMatrix::setDiagonalStore(double * store) {
	THROW(OtherException, wrongMatrixType);
}

}

#endif
