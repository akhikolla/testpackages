/***************************************************************************
                             SRC/mixmod/Kernel/Parameter/GaussianDiagParameter.cpp  description
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
#include "mixmod/Kernel/Parameter/GaussianDiagParameter.h"
#include "mixmod/Kernel/Parameter/GaussianGeneralParameter.h"
#include "mixmod/Kernel/Parameter/GaussianEDDAParameter.h"
#include "mixmod/Kernel/IO/GaussianData.h"
#include "mixmod/Kernel/Model/Model.h"
#include "mixmod/Kernel/Model/ModelType.h"
#include "mixmod/Matrix/DiagMatrix.h"

#include <memory>

namespace XEM {

/****************/
/* Constructors */
/****************/

GaussianDiagParameter::GaussianDiagParameter() {
	THROW(OtherException, wrongConstructorType);
}

//-------------------------------------------------------------------------------------
// constructor called by XEMModel
//-------------------------------------------------------------------------------------
GaussianDiagParameter::GaussianDiagParameter(Model * iModel, ModelType * iModelType) 
: GaussianEDDAParameter(iModel, iModelType) 
{
	int64_t k;
	_tabLambda = new double [_nbCluster];
	_tabShape = new DiagMatrix*[_nbCluster];
	_W = new DiagMatrix(_pbDimension); //Id
	for (k = 0; k < _nbCluster; k++) {
		_tabLambda[k] = 1.0;
		_tabShape[k] = new DiagMatrix(_pbDimension); //Id 

		_tabSigma[k] = new DiagMatrix(_pbDimension); //Id
		_tabInvSigma[k] = new DiagMatrix(_pbDimension); //Id
		_tabWk[k] = new DiagMatrix(_pbDimension); //Id
	}
}
//-------------------------------------------------------------------------------------
// constructor for ParameterDescription (C.Poli)
//-------------------------------------------------------------------------------------
GaussianDiagParameter::GaussianDiagParameter(int64_t iNbCluster, int64_t iPbDimension, ModelType * iModelType) 
: GaussianEDDAParameter(iNbCluster, iPbDimension, iModelType) 
{
	int64_t k;
	_tabLambda = new double [_nbCluster];
	_tabShape = new DiagMatrix*[_nbCluster];
	_W = new DiagMatrix(_pbDimension); //Id
	for (k = 0; k < _nbCluster; k++) {
		_tabLambda[k] = 1.0;
		_tabShape[k] = new DiagMatrix(_pbDimension); //Id 

		_tabSigma[k] = new DiagMatrix(_pbDimension); //Id
		_tabInvSigma[k] = new DiagMatrix(_pbDimension); //Id
		_tabWk[k] = new DiagMatrix(_pbDimension); //Id
	}
}

//---------------------------------------------------------------------
// copy Constructor
//---------------------------------------------------------------------
GaussianDiagParameter::GaussianDiagParameter(const GaussianDiagParameter * iParameter) 
: GaussianEDDAParameter(iParameter) 
{
	int64_t k;
	_tabLambda = copyTab(iParameter->getTabLambda(), _nbCluster);
	_tabShape = new DiagMatrix*[_nbCluster];
	_W = new DiagMatrix(_pbDimension);
	*_W = iParameter->getW();

	Matrix ** iTabSigma = iParameter->getTabSigma();
	Matrix ** iTabInvSigma = iParameter->getTabInvSigma();
	Matrix ** iTabWk = iParameter->getTabWk();
	DiagMatrix ** iTabShape = iParameter->getTabShape();

	for (k = 0; k < _nbCluster; k++) {
		_tabSigma[k] = new DiagMatrix(_pbDimension);
		(* _tabSigma[k]) = iTabSigma[k];
		_tabInvSigma[k] = new DiagMatrix(_pbDimension);
		(* _tabInvSigma[k]) = iTabInvSigma[k];
		_tabWk[k] = new DiagMatrix(_pbDimension);
		(* _tabWk[k]) = iTabWk[k];
		_tabShape[k] = new DiagMatrix(_pbDimension);
		(* _tabShape[k]) = iTabShape[k];
	}
}

/**************/
/* Destructor */
/**************/
GaussianDiagParameter::~GaussianDiagParameter() {
	int64_t k;

	if (_tabLambda) {
		delete[] _tabLambda;
		_tabLambda = NULL;
	}

	if (_tabShape) {
		for (k = 0; k < _nbCluster; k++) {
			delete _tabShape[k];
			_tabShape[k] = NULL;
		}
		delete[] _tabShape;
		_tabShape = NULL;
	}

	if (_tabInvSigma) {
		for (k = 0; k < _nbCluster; k++) {
			delete _tabInvSigma[k];
			_tabInvSigma[k] = NULL;
		}
	}

	if (_tabSigma) {
		for (k = 0; k < _nbCluster; k++) {
			delete _tabSigma[k];
			_tabSigma[k] = NULL;
		}
	}
}

//------------------------
// reset to default values
//------------------------
void GaussianDiagParameter::reset() {
	int64_t k;
	for (k = 0; k < _nbCluster; k++) {
		_tabLambda[k] = 1.0;
		*(_tabShape[k]) = 1.0;
	}
	GaussianEDDAParameter::reset();
}

/*********/
/* clone */
/*********/
Parameter * GaussianDiagParameter::clone() const {
	GaussianDiagParameter * newParam = new GaussianDiagParameter(this);
	return (newParam);
}

/************/
/* initUSER */
/************/
void GaussianDiagParameter::initUSER(Parameter * iParam) {
	GaussianEDDAParameter::initUSER(iParam);
	updateTabInvSigmaAndDet();
}

/*******************/
/* computeTabSigma */
/*******************/
void GaussianDiagParameter::computeTabSigma() {
	/* Initialization */
	GaussianData * data = _model->getGaussianData();
	double * tabNk = _model->getTabNk();
	int64_t k;
	//DiagMatrix * B = new DiagMatrix(_pbDimension);
	//DiagMatrix * Bk = new DiagMatrix(_pbDimension);
	std::unique_ptr<DiagMatrix> B(new DiagMatrix(_pbDimension));
	std::unique_ptr<DiagMatrix>Bk(new DiagMatrix(_pbDimension));
	double detB = 0.0; // Determinant of matrix B
	int64_t iter = 5; // Number of iterations in iterative procedure 
	//  double detDiagW;         // Determinant of diagonal matrix W
	double detDiagWk; // Determinant of diagonal matrix Wk
	double lambda = 0.0; // Volume
	double weightTotal = data->_weightTotal;
	double power = 1.0 / _pbDimension;
	double det = 0.0;
	double logDet = 0.0;
	//double * W_k = new double[_pbDimension];
	std::unique_ptr<double[]> W_k(new double[_pbDimension]);    
	double * Shape_k;
	double tmp;
	int64_t p;

	// Compute det[diag(W)]
	NumericException error = NumericException(minDeterminantDiagWValueError);
	logDet = _W->determinant(error);
	//  detDiagW = powAndCheckIfNotNull(logDet,power);

	// Variance estimator for each of diagonal model
	switch (_modelType->_nameModel) {

	case (Gaussian_p_L_B):
	case (Gaussian_pk_L_B):

		for (k = 0; k < _nbCluster; k++) {
			/*_tabLambda[k]  = detDiagW / weightTotal;
			if ( _tabLambda[k]< minOverflow) {
			  THROW(NumericException,errorSigmaConditionNumber);
		   }*/
			//_tabSigma[k] = _W/weightTotal
			_tabSigma[k]->equalToMatrixDividedByDouble(_W, weightTotal);
		}
		break;

	case (Gaussian_p_L_Bk):
	case (Gaussian_pk_L_Bk):

		for (k = 0; k < _nbCluster; k++) {
			NumericException error = NumericException(minDeterminantDiagWkValueError);
			det = _tabWk[k]->determinant(error);
			detDiagWk = powAndCheckIfNotNull(det, power);

			//_tabShape[k] = _tabW[k]/detWk
			_tabShape[k]->equalToMatrixDividedByDouble(_tabWk[k], detDiagWk);
			lambda += detDiagWk;
		}

		for (k = 0; k < _nbCluster; k++) {
			_tabLambda[k] = lambda / weightTotal;
			if (_tabLambda[k] < minOverflow)
				THROW(NumericException, errorSigmaConditionNumber);

			//tabSigma[k] = tabShape[k]*somme(lambda[k])
			_tabSigma[k]->equalToMatrixMultiplyByDouble(_tabShape[k], _tabLambda[k]);
		}
		break;

	case (Gaussian_p_Lk_Bk):
	case (Gaussian_pk_Lk_Bk):
		
		for (k = 0; k < _nbCluster; k++) {
			NumericException error = NumericException(minDeterminantDiagWkValueError);
			logDet = _tabWk[k]->determinant(error);
			detDiagWk = powAndCheckIfNotNull(logDet, power);

			_tabLambda[k] = detDiagWk / tabNk[k];
			if (_tabLambda[k] < minOverflow)
				THROW(NumericException, errorSigmaConditionNumber);

			//_tabShape[k] = _tabW[k]/detWk
			_tabShape[k]->equalToMatrixDividedByDouble(_tabWk[k], detDiagWk);

			//tabSigma[k] = tabShape[k]*tabLambda[k]
			_tabSigma[k]->equalToMatrixMultiplyByDouble(_tabShape[k], _tabLambda[k]);
		}
		break;

	case (Gaussian_p_Lk_B):
	case (Gaussian_pk_Lk_B):
		
		while (iter) {
			/* Pb Overflow */
			for (k = 0; k < _nbCluster; k++) {
				if (_tabLambda[k] < minOverflow)
					THROW(NumericException, errorSigmaConditionNumber);
			}

			/* Compute matrix B */
			(*B) = 0.0;
			for (k = 0; k < _nbCluster; k++) {
				//_tabShape[k] = _tabW[k]/_tabLambda[k]
				Bk->equalToMatrixDividedByDouble(_tabWk[k], _tabLambda[k]);
				(*B) += Bk.get(); // B->operator+=(Bk) : B=somme sur k des Bk
			}
			/* Compute det(B) */
			NumericException error = NumericException(minDeterminantBValueError);
			logDet = B->determinant(error);
			detB = powAndCheckIfNotNull(logDet, power);

			/* Compute Shape[k] and Lambda[k] */
			for (k = 0; k < _nbCluster; k++) {
				// W_k      = _tabWk[k]->getDiagonalValue();
              _tabWk[k]->putDiagonalValueInStore(W_k.get());
				_tabShape[k]->equalToMatrixDividedByDouble(B.get(), detB);
				Shape_k = _tabShape[k]->getStore();
				tmp = 0.0;
				for (p = 0; p < _pbDimension; p++) {
					tmp += W_k[p] / Shape_k[p];
				}
				tmp /= (_pbDimension * tabNk[k]);

				_tabLambda[k] = tmp;
				if (_tabLambda[k] < minOverflow)
					THROW(NumericException, errorSigmaConditionNumber);
			}
			iter--;
		}

		for (k = 0; k < _nbCluster; k++) {
			_tabSigma[k]->equalToMatrixMultiplyByDouble(_tabShape[k], _tabLambda[k]);
		}

		break;

	default:
		//------
		THROW(OtherException, internalMixmodError);
		break;
	}

	updateTabInvSigmaAndDet();
	//delete Bk;
	//delete B;
	//delete [] W_k;
}

//-------------------
//getLogLikelihoodOne
//-------------------
double GaussianDiagParameter::getLogLikelihoodOne() const {

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
	//DiagMatrix * Sigma = new DiagMatrix(_pbDimension);
	//DiagMatrix * W = new DiagMatrix(_pbDimension, 0.0);
	std::unique_ptr<DiagMatrix> Sigma(new DiagMatrix(_pbDimension));
	std::unique_ptr<DiagMatrix> W(new DiagMatrix(_pbDimension, 0.0));
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
	//logDet   = W->detDiag(minDeterminantDiagWValueError); //virtual
	Sigma->equalToMatrixDividedByDouble(W.get(), totalWeight); //virtual

	// inverse of Sigma
	//DiagMatrix * SigmaMoins1 = new DiagMatrix(_pbDimension);
	Matrix * SigmaMoins1_p = NULL;
	//SigmaMoins1->inverse(Sigma);// virtual
	Sigma->inverse(SigmaMoins1_p); //construction de SigmaMoins1 dans la fonction inverse
    std::unique_ptr<Matrix> SigmaMoins1(SigmaMoins1_p);
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
int64_t GaussianDiagParameter::getFreeParameter() const {
	int64_t nbParameter; // Number of parameters
	int64_t k = _nbCluster; // Sample size
	int64_t d = _pbDimension; // Sample dimension

	int64_t alphaR = k*d; //alpha for for models with Restrainct proportions (Gaussian_p_...)
	int64_t alphaF = (k * d) + k - 1; //alpha for models with Free proportions (Gaussian_pk_...)

	switch (_modelType->_nameModel) {
	case (Gaussian_p_L_B):
		nbParameter = alphaR + d;
		break;
	case (Gaussian_p_Lk_B):
		nbParameter = alphaR + d + k - 1;
		break;
	case (Gaussian_p_L_Bk):
		nbParameter = alphaR + (k * d) - k + 1;
		break;
	case (Gaussian_p_Lk_Bk):
		nbParameter = alphaR + (k * d);
		break;
	case (Gaussian_pk_L_B):
		nbParameter = alphaF + d;
		break;
	case (Gaussian_pk_Lk_B):
		nbParameter = alphaF + d + k - 1;
		break;
	case (Gaussian_pk_L_Bk):
		nbParameter = alphaF + (k * d) - k + 1;
		break;
	case (Gaussian_pk_Lk_Bk):
		nbParameter = alphaF + (k * d);
		break;
	default:
		THROW(OtherException, internalMixmodError);
		break;
	}
	return nbParameter;
}

}
