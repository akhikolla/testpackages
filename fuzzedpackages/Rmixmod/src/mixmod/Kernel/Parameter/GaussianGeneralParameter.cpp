/***************************************************************************
                             SRC/mixmod/Kernel/Parameter/GaussianGeneralParameter.cpp  description
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
#include "mixmod/Kernel/Parameter/GaussianGeneralParameter.h"
#include "mixmod/Kernel/Parameter/GaussianEDDAParameter.h"
#include "mixmod/Kernel/IO/GaussianData.h"
#include "mixmod/Kernel/Model/Model.h"
#include "mixmod/Kernel/Model/ModelType.h"
#include "mixmod/Matrix/DiagMatrix.h"
#include "mixmod/Matrix/SymmetricMatrix.h"
#include "mixmod/Matrix/GeneralMatrix.h"

namespace XEM {

/****************/
/* Constructors */
/****************/

GaussianGeneralParameter::GaussianGeneralParameter() {
	THROW(OtherException, wrongConstructorType);
}

//-------------------------------------------------------------------------------------
// constructor called by XEMModel
//-------------------------------------------------------------------------------------
GaussianGeneralParameter::GaussianGeneralParameter(Model * iModel, ModelType * iModelType) 
: GaussianEDDAParameter(iModel, iModelType) 
{
	int64_t k;
	_tabShape = new DiagMatrix*[_nbCluster];
	_tabOrientation = new GeneralMatrix*[_nbCluster];
	_tabLambda = new double [_nbCluster];
	_W = new SymmetricMatrix(_pbDimension); //Id

	for (k = 0; k < _nbCluster; k++) {
		_tabShape[k] = new DiagMatrix(_pbDimension); // Id
		_tabOrientation[k] = new GeneralMatrix(_pbDimension); //Id
		_tabLambda[k] = 1.0;
		_tabSigma[k] = new SymmetricMatrix(_pbDimension); //Id
		_tabInvSigma[k] = new SymmetricMatrix(_pbDimension); //Id
		_tabWk[k] = new SymmetricMatrix(_pbDimension); //Id
	}

	__storeDim = _pbDimension * (_pbDimension + 1) / 2;
}

//-------------------------------------------------------------------------------------
//constructeur avec une initialisation USER from file
//-------------------------------------------------------------------------------------
GaussianGeneralParameter::GaussianGeneralParameter(
		int64_t iNbCluster, 
		int64_t iPbDimension, 
		ModelType * iModelType, 
		std::string & iFileName)
: GaussianEDDAParameter(iNbCluster, iPbDimension, iModelType) 
{
	int64_t k;
	__storeDim = _pbDimension * (_pbDimension + 1) / 2;
	_tabShape = new DiagMatrix*[_nbCluster];
	_tabOrientation = new GeneralMatrix*[_nbCluster];
	_tabLambda = new double [_nbCluster];

	for (k = 0; k < _nbCluster; k++) {
		_tabShape[k] = new DiagMatrix(_pbDimension); //Id
		_tabOrientation[k] = new GeneralMatrix(_pbDimension); //Id
		_tabLambda[k] = 1.0;

		// _tabSigma, _tabInvSigma, _tabWk will be initialized in XEMGaussianEDDAparameter
		_tabInvSigma[k] = new SymmetricMatrix(_pbDimension); //Id
		_tabSigma[k] = new SymmetricMatrix(_pbDimension); // Id
		_tabWk[k] = new SymmetricMatrix(_pbDimension); //Id
		*_tabWk[k] = 1.0;
	}
	_W = new SymmetricMatrix(_pbDimension); //Id

	// read parameters in file iFileName//
	if (iFileName.compare("") != 0) {
		std::ifstream paramFile(iFileName.c_str(), ios::in);
		if (!paramFile.is_open()) {
			THROW(InputException, wrongParamFileName);
		}
    input(paramFile);
		paramFile.close();
	}

	updateTabInvSigmaAndDet(); // method of XEMGaussianParameter
}

//-------------------------------------------------------------------------------------
//constructeur avec une initialisation USER from file for heterogeneous models
//-------------------------------------------------------------------------------------
GaussianGeneralParameter::GaussianGeneralParameter(
		int64_t iNbCluster, 
		int64_t iPbDimension, 
		ModelType * iModelType, 
		std::string & iFileName,
    int64_t iNbVariable_binary,
    std::vector< int64_t > inbFactor)
: GaussianEDDAParameter(iNbCluster, iPbDimension, iModelType) 
{
	int64_t k;
	__storeDim = _pbDimension * (_pbDimension + 1) / 2;
	_tabShape = new DiagMatrix*[_nbCluster];
	_tabOrientation = new GeneralMatrix*[_nbCluster];
	_tabLambda = new double [_nbCluster];

	for (k = 0; k < _nbCluster; k++) {
		_tabShape[k] = new DiagMatrix(_pbDimension); //Id
		_tabOrientation[k] = new GeneralMatrix(_pbDimension); //Id
		_tabLambda[k] = 1.0;

		// _tabSigma, _tabInvSigma, _tabWk will be initialized in XEMGaussianEDDAparameter
		_tabInvSigma[k] = new SymmetricMatrix(_pbDimension); //Id
		_tabSigma[k] = new SymmetricMatrix(_pbDimension); // Id
		_tabWk[k] = new SymmetricMatrix(_pbDimension); //Id
		*_tabWk[k] = 1.0;
	}
	_W = new SymmetricMatrix(_pbDimension); //Id

	// read parameters in file iFileName//
	if (iFileName.compare("") != 0) {
		std::ifstream paramFile(iFileName.c_str(), ios::in);
		if (!paramFile.is_open()) {
			THROW(InputException, wrongParamFileName);
		}
    input(paramFile, iNbVariable_binary, inbFactor);
		paramFile.close();
	}

	updateTabInvSigmaAndDet(); // method of XEMGaussianParameter
}
//-------------------------------------------------------------------------------------
//constructeur avec une initialisation USER from containers
//-------------------------------------------------------------------------------------
GaussianGeneralParameter::GaussianGeneralParameter(
		int64_t iNbCluster, 
		int64_t iPbDimension, 
		ModelType * iModelType, 
		double * proportions, 
		double ** means, 
		double *** variances)
: GaussianEDDAParameter(iNbCluster, iPbDimension, iModelType) 
{
	int64_t k;
	__storeDim = _pbDimension * (_pbDimension + 1) / 2;
	_tabShape = new DiagMatrix*[_nbCluster];
	_tabOrientation = new GeneralMatrix*[_nbCluster];
	_tabLambda = new double [_nbCluster];

	for (k = 0; k < _nbCluster; k++) {
		_tabShape[k] = new DiagMatrix(_pbDimension); //Id
		_tabOrientation[k] = new GeneralMatrix(_pbDimension); //Id
		_tabLambda[k] = 1.0;

		// _tabSigma, _tabInvSigma, _tabWk will be initialized in XEMGaussianEDDAparameter
		_tabInvSigma[k] = new SymmetricMatrix(_pbDimension); //Id
		_tabSigma[k] = new SymmetricMatrix(_pbDimension); // Id
		_tabWk[k] = new SymmetricMatrix(_pbDimension); //Id
		*_tabWk[k] = 1.0;
	}
	_W = new SymmetricMatrix(_pbDimension); //Id

	input(proportions, means, variances);

	updateTabInvSigmaAndDet(); // method of XEMGaussianParameter
}

//-------------------------------------------------------------------------------------
//constructeur par copie
//-------------------------------------------------------------------------------------
GaussianGeneralParameter::GaussianGeneralParameter(const GaussianGeneralParameter * iParameter) 
: GaussianEDDAParameter(iParameter) 
{
	int64_t k;
	__storeDim = _pbDimension * (_pbDimension + 1) / 2;
	_tabShape = new DiagMatrix*[_nbCluster];
	_tabOrientation = new GeneralMatrix*[_nbCluster];
	_tabLambda = new double[_nbCluster];
	DiagMatrix ** iTabShape = iParameter->getTabShape();
	GeneralMatrix ** iTabOrientation = iParameter->getTabOrientation();
	double * iTabLambda = iParameter->getTabLambda();
	Matrix ** iTabSigma = iParameter->getTabSigma();
	Matrix ** iTabInvSigma = iParameter->getTabInvSigma();
	Matrix ** iTabWk = iParameter->getTabWk();
	_W = new SymmetricMatrix((SymmetricMatrix *) (iParameter->getW()));

	for (k = 0; k < _nbCluster; k++) {
		//_tabInvSigma[k] = NULL;
		_tabShape[k] = new DiagMatrix(iTabShape[k]); // copy constructor
		_tabOrientation[k] = new GeneralMatrix(iTabOrientation[k]); // copy constructor   
		_tabLambda[k] = iTabLambda[k];

		_tabSigma[k] = new SymmetricMatrix(_pbDimension);
		(* _tabSigma[k]) = iTabSigma[k];
		_tabWk[k] = new SymmetricMatrix(_pbDimension);
		(* _tabWk[k]) = iTabWk[k];
		_tabInvSigma[k] = new SymmetricMatrix(_pbDimension);
		(* _tabInvSigma[k]) = iTabInvSigma[k];
	}
}

/**************/
/* Destructor */
/**************/
GaussianGeneralParameter::~GaussianGeneralParameter() {
	int64_t k;

	if (_tabShape) {
		for (k = 0; k < _nbCluster; k++) {
			delete _tabShape[k];
			_tabShape[k] = NULL;
		}
		delete[] _tabShape;
		_tabShape = NULL;
	}

	if (_tabOrientation) {
		for (k = 0; k < _nbCluster; k++) {
			delete _tabOrientation[k];
			_tabOrientation[k] = NULL;
		}
		delete[] _tabOrientation;
		_tabOrientation = NULL;
	}

	if (_tabLambda) {
		delete[] _tabLambda;
		_tabLambda = NULL;
	}

	if (_tabInvSigma) {
		for (k = 0; k < _nbCluster; k++) {
			delete _tabInvSigma[k];
			_tabInvSigma[k] = NULL;
		}
	}

	if (_tabSigma) {
		//cout<<"destructeur gaussgene"<<endl;
		for (k = 0; k < _nbCluster; k++) {
			delete _tabSigma[k];
			_tabSigma[k] = NULL;
		}
	}
}

//------------------------
// reset to default values
//------------------------
void GaussianGeneralParameter::reset() {
	int64_t k;
	for (k = 0; k < _nbCluster; k++) {
		*(_tabShape[k]) = 1.0;
		*(_tabOrientation[k]) = 1.0;
		_tabLambda[k] = 1.0;
	}
	GaussianEDDAParameter::reset();
}

/*********/
/* clone */
/*********/
Parameter * GaussianGeneralParameter::clone() const {
	GaussianGeneralParameter * newParam = new GaussianGeneralParameter(this);
	return (newParam);
}

/************/
/* initUSER */
/************/
void GaussianGeneralParameter::initUSER(Parameter * iParam) {
	GaussianEDDAParameter::initUSER(iParam);
	updateTabInvSigmaAndDet();

	GaussianGeneralParameter * param = (GaussianGeneralParameter *) iParam;
	DiagMatrix ** iTabShape = param->getTabShape();
	GeneralMatrix ** iTabOrientation = param->getTabOrientation();
	double * iTabLambda = param->getTabLambda();
	int64_t k;
	for (k = 0; k < _nbCluster; k++) {
		(*_tabShape[k]) = (iTabShape[k]); // copy constructor
		(* _tabOrientation[k]) = iTabOrientation[k]; // copy constructor
		_tabLambda[k] = iTabLambda[k];
	}
}

/***********************/
/* computeTabSigma_L_C */
/***********************/
void GaussianGeneralParameter::computeTabSigma_L_C() {
	int64_t k;
	double totalWeight = (_model->getGaussianData())->_weightTotal;
	for (k = 0; k < _nbCluster; k++) {
		//Sigma_k = W / totalWeight
		_tabSigma[k]->equalToMatrixDividedByDouble(_W, totalWeight);
	}
}

/*************************/
/* computeTabSigma_Lkk_Ck */
/*************************/
void GaussianGeneralParameter::computeTabSigma_Lk_Ck() {
	int64_t k;
	double * tabNk = _model->getTabNk();
	for (k = 0; k < _nbCluster; k++) {
		// Sigma_k = W_k / n_k
		_tabSigma[k]->equalToMatrixDividedByDouble(_tabWk[k], tabNk[k]);
	}
}

/************************/
/* computeTabSigma_L_Ck */
/************************/
void GaussianGeneralParameter::computeTabSigma_L_Ck() {
	double lambda = 0.0;
	int64_t k;
	double logDet;
	GaussianData * data = _model->getGaussianData();
	double totalWeight = data->_weightTotal;
	double detWk_k, tmp;
	double * detWk = new double[_nbCluster];

	try {
		for (k = 0; k < _nbCluster; k++) {
			NumericException error = NumericException(minDeterminantWkValueError);
			logDet = _tabWk[k]->determinant(error);
			detWk_k = powAndCheckIfNotNull(logDet, 1.0 / _pbDimension);
			lambda += detWk_k;
			detWk[k] = detWk_k;
		}

		lambda /= totalWeight; //lambda = somme sur k detWk^(1/pbDimension) / totalWeight
		if (lambda < minOverflow) {
			THROW(NumericException, errorSigmaConditionNumber);
		}

		for (k = 0; k < _nbCluster; k++) {
			_tabLambda[k] = lambda;
			tmp = detWk[k] / lambda;
			_tabSigma[k]->equalToMatrixDividedByDouble(_tabWk[k], tmp); //Sigma_k = Wk / tmp
		}

		delete [] detWk;
		detWk = NULL;
	}
	catch (Exception&errorType) {
		delete [] detWk;
		detWk = NULL;
		throw;
	}
}

/*****************************/
/* computeTabSigma_L_Dk_A_Dk */
/*****************************/
void GaussianGeneralParameter::computeTabSigma_L_Dk_A_Dk() {

	//Ma modification
	int64_t k;
	GaussianData * data = _model->getGaussianData();
	DiagMatrix* Shape = new DiagMatrix(_pbDimension, 0.0); //0.0
	double logDet;
	double detShape;
	double totalWeight = data->_weightTotal;
	// SVD Decomposition of Cluster Scattering Matrix (Symmetric) Wk: Wk=U*S*U'

	try {

		for (k = 0; k < _nbCluster; k++) {
			_tabWk[k]->computeSVD(_tabShape[k], _tabOrientation[k]);
			//      tabShape_k_store = _tabShape[k]->getStore();
			(*Shape) += _tabShape[k];
		}

		// Compute determinant of Shape matrix 
		NumericException error = NumericException(minDeterminantShapeValueError);
		logDet = Shape->determinant(error);
		detShape = powAndCheckIfNotNull(logDet, 1.0 / _pbDimension);

		for (k = 0; k < _nbCluster; k++) {
			//determinant de la somme sur k de _tabShape[k] / totalWeight
			_tabLambda[k] = detShape / totalWeight;
			if (_tabLambda[k] < minOverflow) {
				THROW(NumericException, errorSigmaConditionNumber);
			}

			// _tabShape[k] = somme sur k de _tabShape / det. de la somme sur k de _tabShape[k]
			_tabShape[k]->equalToMatrixDividedByDouble(Shape, detShape);
			_tabSigma[k]->compute_as__multi_O_S_O(_tabLambda[k], _tabOrientation[k], _tabShape[k]);
		}

		delete Shape;
	}
	catch (Exception&errorType) {
		delete Shape;
		throw;
	}
}

/********************************/
/* computeTabSigma_Lk_Dk_A_Dk() */
/********************************/
void GaussianGeneralParameter::computeTabSigma_Lk_Dk_A_Dk() {
	int64_t k;
	int64_t iter = 5;

	DiagMatrix* Shape = new DiagMatrix(_pbDimension); //Id
	DiagMatrix* S = new DiagMatrix(_pbDimension); //Id
	double * tabNk = _model->getTabNk();

	DiagMatrix ** tabS = new DiagMatrix*[_nbCluster];
	GeneralMatrix ** tabU = new GeneralMatrix*[_nbCluster];
	for (k = 0; k < _nbCluster; k++) {
		tabS[k] = new DiagMatrix(_pbDimension); //Id
		tabU[k] = new GeneralMatrix(_pbDimension); //Id
	}

	double logDet;
	double detShape;

	try {
		// SVD Decomposition of Cluster Scattering Matrix (Symmetric) Wk: Wk=U*S*U'
		for (k = 0; k < _nbCluster; k++) {
			_tabWk[k] ->computeSVD(tabS[k], tabU[k]);
		}
		while (iter) {

			(*Shape) = 0.0;
			for (k = 0; k < _nbCluster; k++) {
				(S)->equalToMatrixDividedByDouble(tabS[k], _tabLambda[k]);
				(*Shape) += S;
			}

			//Compute determinant of Shape matrix 
			NumericException error = NumericException(minDeterminantShapeValueError);
			logDet = Shape->determinant(error);
			detShape = powAndCheckIfNotNull(logDet, 1.0 / _pbDimension);

			for (k = 0; k < _nbCluster; k++) {
				//A=_tabShape = somme sur k 1/lambdak*Shape
				_tabShape[k]->equalToMatrixDividedByDouble(Shape, detShape);
				_tabLambda[k] = _tabWk[k]->trace_this_O_Sm1_O(tabU[k], _tabShape[k]);
				_tabLambda[k] /= (_pbDimension * tabNk[k]);

				if (_tabLambda[k] < minOverflow) {
					THROW(NumericException, errorSigmaConditionNumber);
				}
			}
			iter--;
		}

		for (k = 0; k < _nbCluster; k++) {

			//trace(Wk Dk A^(-1) Dk') / pbDimension*tabNk
			_tabLambda[k] = _tabWk[k]->trace_this_O_Sm1_O(_tabOrientation[k], _tabShape[k]);
			_tabLambda[k] /= (_pbDimension * tabNk[k]);

			if (_tabLambda[k] < minOverflow) {
				THROW(NumericException, errorSigmaConditionNumber);
			}

			(*_tabOrientation[k]) = (tabU[k]); //Dk = tabOrientationk = tabUk
			_tabSigma[k]->compute_as__multi_O_S_O(_tabLambda[k], _tabOrientation[k], _tabShape[k]);
		}

		for (k = 0; k < _nbCluster; k++) {
			delete tabS[k];
			tabS[k] = NULL;
			delete tabU[k];
			tabU[k] = NULL;
		}
		delete [] tabU;
		delete [] tabS;
		delete Shape;
		delete S;
	}

	catch (Exception&errorType) {
		for (k = 0; k < _nbCluster; k++) {
			delete tabS[k];
			tabS[k] = NULL;
			delete tabU[k];
			tabU[k] = NULL;
		}
		delete [] tabU;
		delete [] tabS;
		delete Shape;
		delete S;
		throw;
	}
}

/************************/
/* computeTabSigma_Lk_C */
/************************/
void GaussianGeneralParameter::computeTabSigma_Lk_C() {

	int64_t k;
	double * tabNk = _model->getTabNk();
	SymmetricMatrix * C = new SymmetricMatrix(_pbDimension); //Id
	Matrix * R = new SymmetricMatrix(_pbDimension);
	Matrix * Inv;
	Inv = new SymmetricMatrix(_pbDimension);
	double logDet, detR;
	int64_t iter = 5;

	try {

		while (iter) {
			// Inv = NULL;
			(*R) = 0.0;
			for (k = 0; k < _nbCluster; k++) {
				R->compute_product_Lk_Wk(_tabWk[k], _tabLambda[k]); // on calcule Wk/Lk
			}

			NumericException error = NumericException(minDeterminantRValueError);
			logDet = R->determinant(error);
			detR = powAndCheckIfNotNull(logDet, 1.0 / _pbDimension);
			// C = somme sur k de 1/lambdak * Wk / det (somme sur k de 1/lambdak * Wk )^(1/d)
			C->equalToMatrixDividedByDouble(R, detR);

			C->inverse(Inv);
			for (k = 0; k < _nbCluster; k++) {
				//lambdak = trace(Wk C^(-1)) / pbDimension * tabNk
				_tabLambda[k] = _tabWk[k]->compute_trace_W_C(Inv);
				_tabLambda[k] /= (_pbDimension * tabNk[k]);

				if (_tabLambda[k] < minOverflow) {
					THROW(NumericException, errorSigmaConditionNumber);
				}
			}
			iter--;
		}

		for (k = 0; k < _nbCluster; k++) {
			(*_tabSigma[k]) = (C);
			(*_tabSigma[k]) *= _tabLambda[k];
		}

		delete C;
		C = NULL;
		delete R;
		R = NULL;
		delete Inv;
		Inv = NULL;
	}
	catch (Exception&errorType) {
		delete C;
		C = NULL;
		delete R;
		R = NULL;
		delete Inv;
		Inv = NULL;
		throw;
	}
}

/******************************/
/* computeTabSigma_L_D_Ak_D() */
/******************************/
void GaussianGeneralParameter::computeTabSigma_L_D_Ak_D() {
	GaussianData * data = _model->getGaussianData();
	int64_t k;
	int64_t iter = 5;
	double diff, FOld, F;
	double sumTraceM = 0.0;
	double lambda, logDet, DetDiagQtmp;
	DiagMatrix * Shape_0 = new DiagMatrix(_pbDimension); //Id

	/* Flury algorithm                                                          */
	/* SVD Decomposition of Cluster Scattering Matrix (Symmetric) Wk: Wk=U*S*U' */
	/* Orientation matrix initialisation: _tabOrientation                       */
	try {
		(*Shape_0) = (_tabShape[0]);
		_tabWk[0]-> computeSVD(_tabShape[0], _tabOrientation[0]);
		(*_tabShape[0]) = (Shape_0);

		/* Compute Lambda: _lambda */
		/* Mtmp = D*Ak^-1*D'*Wk    */

		for (k = 0; k < _nbCluster; k++) {
			sumTraceM += _tabWk[k]->trace_this_O_Sm1_O(_tabOrientation[0], _tabShape[k]);
		}
		// lambda = somme sur k de trace(D Ak^(-1) D' Wk) / weightTotal
		lambda = sumTraceM / (_pbDimension * data->_weightTotal);

		diff = 10.0;
		FOld = 0.0;
		F = 0.0;
		/* Iterative procedure 1 */
		while ((iter) && (diff > defaultFluryEpsilon)) {
			/* Compute Shape matrix: _tabShape */
			/* Mtmp = D*Wk^-1*D'               */

			F = 0.0;

			for (k = 0; k < _nbCluster; k++) {
				F += _tabWk[k]->trace_this_O_Sm1_O(_tabOrientation[0], _tabShape[k]);

				_tabWk[k]->computeShape_as__diag_Ot_this_O(_tabShape[k], _tabOrientation[0]);
				NumericException error = NumericException(minDeterminantDiagQtmpValueError);
				logDet = _tabShape[k]->determinant(error);
				DetDiagQtmp = powAndCheckIfNotNull(logDet, 1.0 / _pbDimension);
				// tabShpaek = Ak = diag(D' Wk D) / det(diag(D' Wk D))^(1/d)
				(*_tabShape[k]) /= DetDiagQtmp;
			}

			FOld = F;
			F = flury(F);
			diff = fabs(F - FOld);
			iter--;
		}

		if (lambda < minOverflow) {
			THROW(NumericException, errorSigmaConditionNumber);
		}

		for (k = 0; k < _nbCluster; k++) {
			_tabLambda[k] = lambda;
			_tabSigma[k]->compute_as__multi_O_S_O(lambda, _tabOrientation[0], _tabShape[k]);
		}
		delete Shape_0;

	}
	catch (Exception&errorType) {
		delete Shape_0;
		Shape_0 = NULL;
		throw;
	}
}

/*****************************/
/* computeTabSigma_Lk_D_Ak_D */
/*****************************/
void GaussianGeneralParameter::computeTabSigma_Lk_D_Ak_D() {

	double * tabNk = _model->getTabNk();
	int64_t k, iter = 5;

	/* Flury algorithm                                                          */
	/* SVD Decomposition of Cluster Scattering Matrix (Symmetric) Wk: Wk=U*S*U' */
	//((SymmetricMatrix *) _tabWk[0]) -> computeSVD(_tabShape[0], _tabOrientation[0]);

	_tabWk[0] -> computeSVD(_tabShape[0], _tabOrientation[0]);

	double diff = 10.0;
	double F = 0.0;
	double FOld;

	/* Iterative procedure 1 */
	while ((iter) && (diff > defaultFluryEpsilon)) {
		/* Compute Shape matrix: _tabShape[k] */
		/* Qtmp = D*Wk^-1*D'                  */
		for (k = 0; k < _nbCluster; k++) {
			_tabWk[k]->computeShape_as__diag_Ot_this_O(
					_tabShape[k], _tabOrientation[0], tabNk[k]);
			//logDet = (_tabShape[k])->determinant(minDeterminantShapeValueError);
		}
		FOld = F;
		F = flury(F);
		diff = fabs(F - FOld);
		iter--;
	}

	for (k = 0; k < _nbCluster; k++) {
		(*_tabOrientation[k]) = _tabOrientation[0];
		//((SymmetricMatrix *)_tabSigma[k])
		//    ->compute_as__multi_O_S_O(1.0, _tabOrientation[k], _tabShape[k]);
		_tabSigma[k]->compute_as__multi_O_S_O(1.0, _tabOrientation[k], _tabShape[k]);
	}
}

/*******************/
/* computeTabSigma */
/*******************/
//------------------------------------------------/
/* Variance estimator for each of general model  */
//------------------------------------------------/
void GaussianGeneralParameter::computeTabSigma() {
	switch (_modelType->_nameModel) {

	case (Gaussian_p_L_C):
	case (Gaussian_pk_L_C):
		computeTabSigma_L_C();
		break;
	case (Gaussian_p_Lk_Ck):
	case (Gaussian_pk_Lk_Ck):
		computeTabSigma_Lk_Ck();
		break;
	case (Gaussian_p_L_Ck):
	case (Gaussian_pk_L_Ck):
		computeTabSigma_L_Ck();
		break;
	case (Gaussian_p_L_Dk_A_Dk):
	case (Gaussian_pk_L_Dk_A_Dk):
		computeTabSigma_L_Dk_A_Dk();
		break;
	case (Gaussian_p_Lk_Dk_A_Dk):
	case (Gaussian_pk_Lk_Dk_A_Dk):
		computeTabSigma_Lk_Dk_A_Dk();
		break;
	case (Gaussian_p_Lk_C):
	case (Gaussian_pk_Lk_C):
		computeTabSigma_Lk_C();
		break;
	case (Gaussian_p_L_D_Ak_D):
	case (Gaussian_pk_L_D_Ak_D):
		computeTabSigma_L_D_Ak_D();
		break;
	case (Gaussian_p_Lk_D_Ak_D):
	case (Gaussian_pk_Lk_D_Ak_D):
		computeTabSigma_Lk_D_Ak_D();
		break;
	default:
		THROW(OtherException, internalMixmodError);
		break;
	}

	updateTabInvSigmaAndDet();
}

//----------------
// Flury Algorithm
//----------------
// TODO: provide references
double GaussianGeneralParameter::flury(double F) {

	// coefficients de la matrice a diagonaliser
	//       [ a     b ]
	//  M  = [         ]
	//       [ b     c ]
	double a, b, c;
	// plus petite valeur propre
	double eigenValueMoins;
	// vecteur propre associe a eigenValueMoins
	double eigenVectorMoins_1;
	double eigenVectorMoins_2;
	int64_t k, il, im, iter = 0;
	double diff = 10;
	double FOld, tmp, tmp2, termesHorsDiag, termesDiag;
	double * Wk_store;
	int64_t p, q, r;
	double * tabShape_k_store;
	double * Ori = _tabOrientation[0]->getStore();
	GaussianData * data = _model->getGaussianData();
	double * D_im = data->getTmpTabOfSizePbDimension();
	//double * D_il = new double[_pbDimension];
	std::unique_ptr<double[]> D_il(new double[_pbDimension]);    
	int64_t il_p, im_p;
	double iSim_iSil_k;

	while ((iter < maxFluryIter) && (diff > defaultFluryEpsilon)) {

		/* Compute orientation matrix: _tabOrientation */
		for (il = 0; il < _pbDimension; il++) {
			for (im = il + 1; im < _pbDimension; im++) {

				/* Extraction of colums dl and dm */
				il_p = il;
				im_p = im;
				for (p = 0; p < _pbDimension; p++) {
					D_il[p] = Ori[il_p]; // remplir la colonne 1 :         
					D_im[p] = Ori[im_p]; // remplir la colonne 2 :

					il_p += _pbDimension;
					im_p += _pbDimension;
				}

				/* Compute matrix R where q1 is the second eigen vector */
				//RR_store = RR.Store();
				//initToZero(RR_store,4);

				a = 0.0;
				b = 0.0;
				c = 0.0;

				for (k = 0; k < _nbCluster; k++) {

					//Wk_store         = ((SymmetricMatrix *)_tabWk[k])->getStore();
					Wk_store = _tabWk[k]->getSymmetricStore();
					tabShape_k_store = _tabShape[k]->getStore();
					iSim_iSil_k = 1.0 / tabShape_k_store[il] - 1.0 / tabShape_k_store[im];

					// remplir RR(1,1)
					termesHorsDiag = termesDiag = 0.0;
					for (p = 0, r = 0; p < _pbDimension; p++, r++) {
						tmp = D_il[p]; // Dlm(p+1,1);
						for (q = 0; q < p; q++, r++) {
							termesHorsDiag += Wk_store[r] * tmp * D_il[q]; //(q+1,1);
						}
						termesDiag += Wk_store[r] * tmp * tmp;
					}
					termesDiag += 2 * termesHorsDiag;
					termesDiag *= iSim_iSil_k;
					a += termesDiag;

					// remplir RR(1,2)
					termesHorsDiag = termesDiag = 0.0;
					for (p = 0, r = 0; p < _pbDimension; p++, r++) {
						tmp = D_il[p];
						tmp2 = D_im[p];
						for (q = 0; q < p; q++, r++) {
							termesHorsDiag += Wk_store[r] * (tmp * D_im[q] + D_il[q] * tmp2);
						}
						termesDiag += Wk_store[r] * tmp * tmp2;
					}
					termesDiag += termesHorsDiag;
					termesDiag *= iSim_iSil_k;
					b += termesDiag;

					// remplir RR(2,2)   
					termesHorsDiag = termesDiag = 0.0;
					for (p = 0, r = 0; p < _pbDimension; p++, r++) {
						tmp = D_im[p];
						for (q = 0; q < p; q++, r++) {
							termesHorsDiag += Wk_store[r] * tmp * D_im[q]; //(q+1,2);
						}
						termesDiag += Wk_store[r] * tmp * tmp;
					}
					termesDiag += 2 * termesHorsDiag;
					termesDiag *= iSim_iSil_k;
					c += termesDiag;

				}

				tmp = a - c;
				eigenValueMoins = (a + c - sqrt(tmp * tmp + 4 * b * b)) / 2.0;

				eigenVectorMoins_1 = b 
						/ sqrt(b * b + (a - eigenValueMoins)*(a - eigenValueMoins));
				eigenVectorMoins_2 = -(a - eigenValueMoins) 
						/ sqrt(b * b + (a - eigenValueMoins)*(a - eigenValueMoins));

				// remplir la colonne d'indice il avec Dlm * q1
				il_p = il;
				im_p = im;
				for (p = 0; p < _pbDimension; p++) {
					Ori[il_p] = D_il[p] * eigenVectorMoins_1 + D_im[p] * eigenVectorMoins_2;
					Ori[im_p] = D_il[p] * eigenVectorMoins_2 - D_im[p] * eigenVectorMoins_1;
					il_p += _pbDimension;
					im_p += _pbDimension;
				}

			}
		}

		FOld = F;
		F = 0.0;

		for (k = 0; k < _nbCluster; k++) {
			// W_k    = (SymmetricMatrix *)(_tabWk[k]);
			// F     += W_k->trace_this_O_Sm1_O(_tabOrientation[0], _tabShape[k]);
			F += _tabWk[k]->trace_this_O_Sm1_O(_tabOrientation[0], _tabShape[k]);
		}

		diff = fabs(F - FOld);

		iter++;
	}

	//delete[] D_il;
	return F;
}

//-------------------
//getLogLikelihoodOne
//-------------------
double GaussianGeneralParameter::getLogLikelihoodOne() const {
	/* Compute log-likelihood for one cluster
	   useful for NEC criterion */

	/* Initialization */
	int64_t nbSample = _model->getNbSample();
	int64_t i;
	GaussianData * data = _model->getGaussianData();
	double logLikelihoodOne; // Log-likelihood for k=1
	//double * Mean = new double[_pbDimension];
	std::unique_ptr<double[]> Mean(new double[_pbDimension]);
	double ** y = data->_yStore;
	double * yi;
	//SymmetricMatrix * Sigma = new SymmetricMatrix(_pbDimension); //Id
	//SymmetricMatrix * W = new SymmetricMatrix(_pbDimension, 0.0); // 0.0
	std::unique_ptr<SymmetricMatrix> Sigma(new SymmetricMatrix(_pbDimension)); //Id
	std::unique_ptr<SymmetricMatrix> W(new SymmetricMatrix(_pbDimension, 0.0)); // 0.0
	double norme;
	double * weight = data->_weight;

	//  Mean Estimator (empirical estimator)
	double totalWeight = data->_weightTotal;
	computeMeanOne(Mean.get(), weight, y, nbSample, totalWeight);
	weight = data->_weight;

	/* Compute the Cluster Scattering Matrix W */

	int64_t p; // parcours
	double * xiMoinsMuk = data->getTmpTabOfSizePbDimension();

	for (i = 0; i < nbSample; i++) {
		yi = y[i];
		for (p = 0; p < _pbDimension; p++) {
			xiMoinsMuk[p] = yi[p] - Mean[p];
		}
		W->add(xiMoinsMuk, weight[i]);
	}
	/* Compute determinant of diag(W) */
	//  logDet   = W->detDiag(minDeterminantDiagWValueError);  // virtual
	//  detDiagW = powAndCheckIfNotNull(logDet ,1.0/_pbDimension);

	Sigma->equalToMatrixDividedByDouble(W.get(), totalWeight); // virtual

	// inverse of Sigma

	//SymmetricMatrix * SigmaMoins1 = new SymmetricMatrix(_pbDimension);
	Matrix* SigmaMoins1_p = NULL;
	//SigmaMoins1->inverse(Sigma);// virtual
	Sigma->inverse(SigmaMoins1_p);
    std::unique_ptr<Matrix> SigmaMoins1(SigmaMoins1_p); //only for a safe deletion of SigmaMoins1_p
	NumericException error = NumericException(minDeterminantSigmaValueError);
	double detSigma = Sigma->determinant(error); // virtual

	// Compute the log-likelihood for one cluster (k=1)
	logLikelihoodOne = 0.0;
	for (i = 0; i < nbSample; i++) {
		yi = y[i];
		for (p = 0; p < _pbDimension; p++) {
			xiMoinsMuk[p] = yi[p] - Mean[p];
		}

		norme = SigmaMoins1->norme(xiMoinsMuk); // virtual
		logLikelihoodOne += norme * weight[i];
	}

	logLikelihoodOne += totalWeight * (data->getPbDimensionLog2Pi() + log(detSigma));
	logLikelihoodOne *= -0.5;

	//delete W;
	//delete Sigma;
	//delete SigmaMoins1;

	//delete[] Mean;

	return logLikelihoodOne;
}

//----------------
//getFreeParameter
//----------------
int64_t GaussianGeneralParameter::getFreeParameter() const {
	int64_t nbParameter; // Number of parameters //
	int64_t k = _nbCluster; // Sample size          //
	int64_t d = _pbDimension; // Sample dimension     //

	int64_t alphaR = k*d; // alpha for for models with Restrainct proportions (Gaussian_p_...)
	int64_t alphaF = (k * d) + k - 1; // alpha for models with Free proportions (Gaussian_pk_...)
	int64_t beta = d * (d + 1) / 2;

	switch (_modelType->_nameModel) {
	case (Gaussian_p_L_C):
		nbParameter = alphaR + beta;
		break;
	case (Gaussian_p_Lk_C):
		nbParameter = alphaR + beta + (k - 1);
		break;
	case (Gaussian_p_L_D_Ak_D):
		nbParameter = alphaR + beta + (k - 1)*(d - 1);
		break;
	case (Gaussian_p_Lk_D_Ak_D):
		nbParameter = alphaR + beta + (k - 1) * d;
		break;
	case (Gaussian_p_L_Dk_A_Dk):
		nbParameter = alphaR + (k * beta) - (k - 1) * d;
		break;
	case (Gaussian_p_Lk_Dk_A_Dk):
		nbParameter = alphaR + (k * beta) - (k - 1)*(d - 1);
		break;
	case (Gaussian_p_L_Ck):
		nbParameter = alphaR + (k * beta) - (k - 1);
		break;
	case (Gaussian_p_Lk_Ck):
		nbParameter = alphaR + (k * beta);
		break;
	case (Gaussian_pk_L_C):
		nbParameter = alphaF + beta;
		break;
	case (Gaussian_pk_Lk_C):
		nbParameter = alphaF + beta + (k - 1);
		break;
	case (Gaussian_pk_L_D_Ak_D):
		nbParameter = alphaF + beta + (k - 1)*(d - 1);
		break;
	case (Gaussian_pk_Lk_D_Ak_D):
		nbParameter = alphaF + beta + (k - 1) * d;
		break;
	case (Gaussian_pk_L_Dk_A_Dk):
		nbParameter = alphaF + (k * beta) - (k - 1) * d;
		break;
	case (Gaussian_pk_Lk_Dk_A_Dk):
		nbParameter = alphaF + (k * beta) - (k - 1)*(d - 1);
		break;
	case (Gaussian_pk_L_Ck):
		nbParameter = alphaF + (k * beta) - (k - 1);
		break;
	case (Gaussian_pk_Lk_Ck):
		nbParameter = alphaF + (k * beta);
		break;
	default:
		THROW(OtherException, internalMixmodError);
		break;
	}
	return nbParameter;
}

}
