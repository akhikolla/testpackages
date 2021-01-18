/***************************************************************************
                             SRC/mixmod/Matrix/Matrix.h  description
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
#ifndef XEMMATRIX_H
#define XEMMATRIX_H

#include "mixmod/Utilities/Util.h"

namespace XEM {

// pre-declaration
class GeneralMatrix;
class DiagMatrix;

/**
  @brief Base class for Matrix
  @author F Langrognet & A Echenim
 */

class Matrix {

public:

	/// Default constructor
	Matrix();

	Matrix(int64_t pbDimension);

	Matrix(Matrix * A);

	/// Desctructor
	virtual ~Matrix();

	int64_t getPbDimension();

	// TODO static
	int64_t _s_pbDimension;


	/// fill this with the inverse matrix of A
	virtual void inverse(Matrix * & A) = 0;

	virtual void compute_product_Lk_Wk(Matrix* Wk, double L) = 0;

	/// compute (x - mean)' this (x - mean) 
	virtual double norme(double * xMoinsMean) = 0;

	/// compute determinant
	virtual double determinant(Exception& errorType) = 0;

	/// compute trace
	virtual double computeTrace() = 0;

	/// this will be A / d
	virtual void equalToMatrixDividedByDouble(Matrix * A, double d) = 0;

	/// this will be A * d
	virtual void equalToMatrixMultiplyByDouble(Matrix*D, double d) = 0;


	/// add :  cik * xMoinsMean * xMoinsMean'  to this
	virtual void add(double * xMoinsMean, double cik) = 0;

	// add : diag( cik * xMoinsMean * xMoinsMean' )  to this
	//virtual void addDiag(double * xMoinsMean, double cik) = 0;

	/// this = d * Identity
	virtual void operator=(const double& d) = 0;
	/// this = this / (d * Identity)
	virtual void operator/=(const double& d) = 0;
	/// this = this * (d * Identity)
	virtual void operator*=(const double& d) = 0;
	/// this = this + matrix
	virtual void operator+=(Matrix* M) = 0;
	/// this = matrix
	virtual void operator=(Matrix* M) = 0;


	void edit(std::ostream& flux, std::string before);
	/// read matrix from an input file
	virtual void input(std::ifstream & fi) = 0;
	virtual void input(double ** variances) = 0;

	// pour ne plus faire de transtypages

	/// return store of a spherical matrix
	virtual double putSphericalValueInStore(double & store) = 0;
	/// add store of a spherical matrix
	virtual double addSphericalValueInStore(double & store) = 0;

	virtual double getSphericalStore() = 0;

	/// Return store of a diagonal matrix
	virtual double* putDiagonalValueInStore(double * store) = 0;
	/// Add store of a diagonal matrix 
	virtual double* addDiagonalValueInStore(double * store) = 0;

	virtual double* getDiagonalStore() = 0;

	/// Return store of a diagonal matrix
	virtual double* putSymmetricValueInStore(double * store) = 0;
	/// Add store of a diagonal matrix in a diagonal one
	virtual double* addSymmetricValueInStore(double * store) = 0;

	virtual double* getSymmetricStore() = 0;

	/// Return store of a diagonal matrix
	virtual double* putGeneralValueInStore(double * store) = 0;
	/// Add store of a diagonal matrix in a diagonal one
	virtual double* addGeneralValueInStore(double * store) = 0;

	virtual double* getGeneralStore() = 0;

	virtual double** storeToArray() const = 0;

	/// gives : det(diag(this))
	virtual double detDiag(Exception& errorType) = 0;

	///compute singular vector decomposition
	virtual void computeSVD(DiagMatrix* & S, GeneralMatrix* & O) = 0;

	virtual void compute_as__multi_O_S_O(double multi, GeneralMatrix* & O, DiagMatrix *& S) = 0;
	virtual double trace_this_O_Sm1_O(GeneralMatrix* & O, DiagMatrix* & S) = 0;
	virtual double compute_trace_W_C(Matrix * C) = 0;
	virtual void computeShape_as__diag_Ot_this_O(DiagMatrix* & Shape, GeneralMatrix* & Ori, double diviseur = 1.0) = 0;
	// temporary table of Size pbDimension, comes from XEMGaussianData
	// used in norme in the General case
	double * _tmpTab;

	///set store
	virtual void setSymmetricStore(double * store) = 0;
	virtual void setGeneralStore(double * store) = 0;
	virtual void setDiagonalStore(double * store) = 0;
	virtual void setSphericalStore(double store) = 0;
};

inline int64_t Matrix::getPbDimension() {
	return _s_pbDimension;
}

}

#endif
