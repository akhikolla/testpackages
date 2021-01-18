/***************************************************************************
                             SRC/mixmod/Kernel/Parameter/GaussianSphericalParameter.cpp  description
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
#include "mixmod/Kernel/Parameter/GaussianSphericalParameter.h"
#include "mixmod/Kernel/Parameter/GaussianGeneralParameter.h"
#include "mixmod/Kernel/Parameter/GaussianEDDAParameter.h"
#include "mixmod/Kernel/IO/GaussianData.h"
#include "mixmod/Kernel/Model/Model.h"
#include "mixmod/Kernel/Model/ModelType.h"
#include "mixmod/Matrix/SphericalMatrix.h"

namespace XEM {

/****************/
/* Constructors */
/****************/

GaussianSphericalParameter::GaussianSphericalParameter() {
	THROW(OtherException, wrongConstructorType);
}

//-------------------------------------------------------------------------------------
// constructor called by XEMModel
//-------------------------------------------------------------------------------------
GaussianSphericalParameter::GaussianSphericalParameter(Model * iModel, ModelType * iModelType) 
: GaussianEDDAParameter(iModel, iModelType) 
{
	int64_t k;
	_W = new SphericalMatrix(_pbDimension);
	for (k = 0; k < _nbCluster; k++) {
		_tabSigma[k] = new SphericalMatrix(_pbDimension); //Id
		_tabInvSigma[k] = new SphericalMatrix(_pbDimension); //Id
		_tabWk[k] = new SphericalMatrix(_pbDimension); //Id
	}
}
//-------------------------------------------------------------------------------------
// constructor for ParameterDescription (C.Poli)
//-------------------------------------------------------------------------------------
GaussianSphericalParameter::GaussianSphericalParameter(int64_t iNbCluster, int64_t iPbDimension, ModelType * iModelType) 
: GaussianEDDAParameter(iNbCluster, iPbDimension, iModelType) 
{
	int64_t k;
	_W = new SphericalMatrix(_pbDimension);
	for (k = 0; k < _nbCluster; k++) {
		_tabSigma[k] = new SphericalMatrix(_pbDimension); //Id
		_tabInvSigma[k] = new SphericalMatrix(_pbDimension); //Id
		_tabWk[k] = new SphericalMatrix(_pbDimension); //Id
	}
}

//---------------------------------------------------------------------
// copy constructor
//---------------------------------------------------------------------
GaussianSphericalParameter::GaussianSphericalParameter(
		const GaussianSphericalParameter * iParameter) 
: GaussianEDDAParameter(iParameter) 
{
	int64_t k;

	_W = new SphericalMatrix((SphericalMatrix *) (iParameter->getW())); // copy constructor

	Matrix ** iTabWk = iParameter->getTabWk();
	Matrix ** iTabSigma = iParameter->getTabSigma();
	Matrix ** iTabInvSigma = iParameter->getTabInvSigma();

	for (k = 0; k < _nbCluster; k++) {
		_tabWk[k] = new SphericalMatrix(_pbDimension);
		(* _tabWk[k]) = iTabWk[k];
		_tabSigma[k] = new SphericalMatrix(_pbDimension);
		(* _tabSigma[k]) = iTabSigma[k];
		_tabInvSigma[k] = new SphericalMatrix(_pbDimension);
		(*_tabInvSigma[k]) = iTabInvSigma[k];
	}
}

/**************/
/* Destructor */
/**************/
GaussianSphericalParameter::~GaussianSphericalParameter() {

	if (_tabSigma) {
		for (int64_t k = 0; k < _nbCluster; k++) {
			delete _tabSigma[k];
			//_tabSigma[k] = NULL;
		}
	}

	if (_tabInvSigma) {
		for (int64_t k = 0; k < _nbCluster; k++) {
			delete _tabInvSigma[k];
			// _tabInvSigma[k] = NULL;
		}
	}
}

/*********/
/* clone */
/*********/
Parameter * GaussianSphericalParameter::clone() const {
	GaussianSphericalParameter * newParam = new GaussianSphericalParameter(this);
	return (newParam);
}

/************/
/* initUSER */
/************/
void GaussianSphericalParameter::initUSER(Parameter * iParam) {
	GaussianEDDAParameter::initUSER(iParam);
	updateTabInvSigmaAndDet();
}

/*******************/
/* computeTabSigma */
/*******************/
void GaussianSphericalParameter::computeTabSigma() {

	// Initialization
	GaussianData * data = _model->getGaussianData();
	double * tabNk = _model->getTabNk();
	int64_t k;
	double sigmaValue;

	double totalWeight = data->_weightTotal;

	// Variance estimator for each of spherical model
	switch (_modelType->_nameModel) {

	case (Gaussian_p_L_I):
	case (Gaussian_pk_L_I):
		//---------------------
		_W->putSphericalValueInStore(sigmaValue); // / totalWeight;
		sigmaValue /= totalWeight;
		if (sigmaValue < minOverflow) {
			THROW(NumericException, errorSigmaConditionNumber);
		}
		for (k = 0; k < _nbCluster; k++) {
			*(_tabSigma[k]) = sigmaValue;
		}

		break;

	case (Gaussian_p_Lk_I):
	case (Gaussian_pk_Lk_I):
		//----------------------
		for (k = 0; k < _nbCluster; k++) {
			_tabWk[k]->putSphericalValueInStore(sigmaValue);
			sigmaValue /= tabNk[k];
			if (sigmaValue < minOverflow)
				THROW(NumericException, errorSigmaConditionNumber);
			*(_tabSigma[k]) = sigmaValue;
		}

		break;

	default:
		THROW(OtherException, internalMixmodError);
		break;
	}

	updateTabInvSigmaAndDet();
}

//-------------------
//getLogLikelihoodOne
//-------------------
double GaussianSphericalParameter::getLogLikelihoodOne() const {

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
	//SphericalMatrix * Sigma = new SphericalMatrix(_pbDimension);
	//SphericalMatrix * W = new SphericalMatrix(_pbDimension);
	std::unique_ptr<SphericalMatrix> Sigma(new SphericalMatrix(_pbDimension));
	std::unique_ptr<SphericalMatrix> W(new SphericalMatrix(_pbDimension));
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
          xiMoinsMuk[p] = yi[p] - Mean.get()[p];
		}
		W->add(xiMoinsMuk, weight[i]);
	}

	/* Compute determinant of diag(W) */
	//  logDet   = W->detDiag(minDeterminantDiagWValueError);  // virtual
	//  detDiagW = powAndCheckIfNotNull(logDet ,1.0/_pbDimension);
	Sigma->equalToMatrixDividedByDouble(W.get(), totalWeight); // virtual

	// inverse of Sigma
	//XEMSphericalMatrix * SigmaMoins1 = new XEMSphericalMatrix( _pbDimension);
	Matrix * SigmaMoins1_p = NULL;
	//SigmaMoins1->inverse(Sigma);// virtual
	//cout<<"S"<<endl;
	//Sigma->edit(cout,"");
	Sigma->inverse(SigmaMoins1_p);//actually SigmaMoins1 are created by inverse() with "New Matrix..."
    std::unique_ptr<Matrix> SigmaMoins1(SigmaMoins1_p); //instead of the final delete (still ugly, althtough better than the final delete...)
	//cout<<"S-1"<<endl;
	//SigmaMoins1->edit(cout,"");
	NumericException error = NumericException(minDeterminantSigmaValueError);
	double detSigma = Sigma->determinant(error); // virtual

	// Compute the log-likelihood for one cluster (k=1)
	logLikelihoodOne = 0.0;
	for (i = 0; i < nbSample; i++) {
		yi = y[i];
		for (p = 0; p < _pbDimension; p++) {
          xiMoinsMuk[p] = yi[p] - Mean.get()[p];
		}
		norme = SigmaMoins1->norme(xiMoinsMuk); // virtual
		logLikelihoodOne += norme * weight[i];
	}

	logLikelihoodOne += totalWeight * (data->getPbDimensionLog2Pi() + log(detSigma));
	logLikelihoodOne *= -0.5;

	//delete[] Mean;

	//delete W;
	//delete Sigma;
	//delete SigmaMoins1;

	return logLikelihoodOne;
}

//----------------
//getFreeParameter
//----------------
int64_t GaussianSphericalParameter::getFreeParameter() const {
	int64_t nbParameter; // Number of parameters //
	int64_t k = _nbCluster; // Sample size          //
	int64_t d = _pbDimension; // Sample dimension     //

	int64_t alphaR = k*d; // alpha for for models with Restrainct proportions (Gaussian_p_...)
	int64_t alphaF = (k * d) + k - 1; // alpha for models with Free proportions (Gaussian_pk_...)

	switch (_modelType->_nameModel) {
	case (Gaussian_p_L_I):
		nbParameter = alphaR + 1;
		break;
	case (Gaussian_p_Lk_I):
		nbParameter = alphaR + k; // correction 09/09/2009 (old version was alphaR + d)
		break;
	case (Gaussian_pk_L_I):
		nbParameter = alphaF + 1;
		break;
	case (Gaussian_pk_Lk_I):
		nbParameter = alphaF + k; // correction 09/09/2009 (old version was alphaR + d)
		break;
	default:
		THROW(OtherException, internalMixmodError);
		break;
	}
	return nbParameter;
}

}
