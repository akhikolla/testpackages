/***************************************************************************
                             SRC/mixmod/Kernel/Parameter/GaussianParameter.cpp  description
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
#include "mixmod/Kernel/Parameter/GaussianParameter.h"
#include "mixmod/Kernel/Parameter/GaussianGeneralParameter.h"
#include "mixmod/Utilities/Util.h"
#include "mixmod/Kernel/Model/Model.h"
#include "mixmod/Kernel/Model/ModelType.h"
#include "mixmod/Kernel/IO/GaussianData.h"
#include "mixmod/Kernel/IO/GaussianSample.h"
#include "mixmod/Kernel/IO/Partition.h"
#include "mixmod/Matrix/Matrix.h"
#include "mixmod/Matrix/DiagMatrix.h"
#include "mixmod/Matrix/SphericalMatrix.h"
#include "mixmod/Matrix/GeneralMatrix.h"
#include "mixmod/Matrix/SymmetricMatrix.h"

namespace XEM {

//------------
// Constructor
//------------
// Default constructor
GaussianParameter::GaussianParameter() {
	THROW(OtherException, wrongConstructorType);
}

//------------
// Constructor
//------------
GaussianParameter::GaussianParameter(Model * iModel, ModelType * iModelType) 
: Parameter(iModel, iModelType) 
{
	int64_t k, j;
	_pbDimension = iModel->getGaussianData()->getPbDimension();
	_tabMean = new double*[_nbCluster];

	_tabWk = new Matrix* [_nbCluster];
	// _tabSigma           = new XEMMatrix* [_nbCluster];

	for (k = 0; k < _nbCluster; k++) {
		_tabMean[k] = new double[_pbDimension];
		for (j = 0; j < _pbDimension; j++) {
			_tabMean[k][j] = 0;
		}
	}

	initFreeProportion(iModelType);
}

//------------
// Constructor
// called by XEMGaussianEDDAParameter if initialization is USER
//------------
GaussianParameter::GaussianParameter(
		int64_t iNbCluster, int64_t iPbDimension, ModelType * iModelType) 
: Parameter(iNbCluster, iPbDimension, iModelType) 
{
	int64_t k;
	_tabMean = new double*[_nbCluster];

	_tabWk = new Matrix* [_nbCluster];
	// _tabSigma           = new XEMMatrix* [_nbCluster];

	for (k = 0; k < _nbCluster; k++) {
		_tabMean[k] = new double[_pbDimension];
		for (int64_t j = 0; j < _pbDimension; j++) {
			_tabMean[k][j] = 0;
		}
	}

	initFreeProportion(iModelType);
}

//------------
// Constructor (copy)
//------------
GaussianParameter::GaussianParameter(const GaussianParameter * iParameter) 
: Parameter(iParameter) 
{
	int64_t k;
	_tabMean = new double*[_nbCluster];
	double ** iTabMean = iParameter->getTabMean();
	for (k = 0; k < _nbCluster; k++) {
		_tabMean[k] = copyTab(iTabMean[k], _pbDimension);
	}
	_tabWk = new Matrix* [_nbCluster];
	//_tabSigma           = new XEMMatrix* [_nbCluster];
}

//-----------
// Destructor
//-----------
GaussianParameter::~GaussianParameter() {
	int64_t k;

	if (_tabMean) {
		for (k = 0; k < _nbCluster; k++) {
			delete[] _tabMean[k];
			_tabMean[k] = NULL;
		}
		delete[] _tabMean;
		_tabMean = NULL;
	}

	if (_W) {
		delete _W;
		_W = NULL;
	}
	if (_tabWk) {
		for (k = 0; k < _nbCluster; k++) {
			delete _tabWk[k];
		}
		delete[] _tabWk;
		_tabWk = NULL;
	}
	/* if(_tabSigma){
	   for(k=0; k<_nbCluster; k++){
		 delete _tabSigma[k];
	   }
	   delete[] _tabSigma;
	   _tabSigma = NULL;
	 }*/
}

//---------------------
/// Comparison operator
//---------------------
bool GaussianParameter::operator ==(const GaussianParameter & param) const {
	if (!Parameter::operator==(param)) return false;
	for (int64_t k = 0; k < _nbCluster; k++) {
		for (int64_t j = 0; j < _pbDimension; j++) {
			if (_tabMean[k][j] != param.getTabMean()[k][j]) return false;
		}
	}
	return true;
}

//------------------------
// reset to default values
//------------------------
void GaussianParameter::reset() {
	int64_t k, j;

	for (k = 0; k < _nbCluster; k++) {
		*(_tabWk[k]) = 1.0;
		for (j = 0; j < _pbDimension; j++) {
			_tabMean[k][j] = 0;
		}
	}

	*(_W) = 1.0;

	Parameter::reset();
}

//-----------------------------
// updateForOneRunOfInitRANDOM
//-----------------------------
void GaussianParameter::updateForInitRANDOMorUSER_PARTITION(
		Sample ** tabSampleForInit, bool * tabClusterToInitialize) 
{
	double * sampleValue = NULL;
	Sample * curSample = NULL;
	for (int64_t k = 0; k < _nbCluster; k++) {
		if (tabClusterToInitialize[k]) {
			curSample = tabSampleForInit[k];
			sampleValue = (curSample->getGaussianSample())->getTabValue();
			recopyTab(sampleValue, _tabMean[k], _pbDimension);
		}
	}
}

//-----------------------------
//computeGlobalDiagDataVariance
//-----------------------------
void GaussianParameter::computeGlobalDiagDataVariance(DiagMatrix * matrixDiagDataVar) {
	int64_t nbSample = _model->getNbSample();
	int64_t p;
	int64_t i;
	GaussianData * data = _model->getGaussianData();
	double totalWeight = data->_weightTotal;
	double ** p_yStore = data->getYStore();
	double * p_yStore_i;
	double * p_weight = data->_weight;
	//double * Mean = new double[_pbDimension];
	std::unique_ptr<double[]> Mean(new double[_pbDimension]);    
	double wi;
	double * xiMoinsMean = data->getTmpTabOfSizePbDimension();

	computeMeanOne(Mean.get(), p_weight, p_yStore, nbSample, totalWeight);

	// compute W diagonal
	p_yStore = data->getYStore();
	p_weight = data->_weight;
	*matrixDiagDataVar = 0.0;
	for (i = 0; i < nbSample; i++) {
		p_yStore_i = *p_yStore;
		wi = *p_weight;
		for (p = 0; p < _pbDimension; p++) {
			xiMoinsMean[p] = p_yStore_i[p] - Mean[p];
		}
		matrixDiagDataVar->add(xiMoinsMean, wi); // virtual
		p_yStore++;
		p_weight++;
	}
	*matrixDiagDataVar /= totalWeight; // virtual operator

	//delete[] Mean;
}

//--------------
// Compute TabWk
//--------------
void GaussianParameter::computeTabWkW() {
	// Compute the cluster scattering matrices Wk and W

	// NB: _tabMean and _zik must be updated
	// this method must be called to update Wk and W if _tabMean or _zik have changed
	// makes calls to virtual methods

	//recuperation des différents paramètres

	double ** tabCik = _model->getTabCik();
	double ** p_cik;
	int64_t nbSample = _model->getNbSample();
	GaussianData * data = _model->getGaussianData();
	double ** p_tabMean = _tabMean;
	double * weight = data->_weight;
	int64_t i;
	int64_t k;
	// storage
	double ** matrix = data->getYStore(); // to store x_i i=1,...,n
	double * xiMoinsMuk = data->getTmpTabOfSizePbDimension(); //to store x_i - mu_k at each step
	double* muk;
	int64_t p;
	double cik;
	*(_W) = 0.0;

	for (k = 0; k < _nbCluster; k++) {
		muk = *p_tabMean;
		(*_tabWk[k]) = 0.0; // virtual method
		p_cik = tabCik;
		matrix = data->getYStore();
		for (i = 0; i < nbSample; i++) {
			cik = (*p_cik)[k] * weight[i];
			for (p = 0; p < _pbDimension; p++) {
				xiMoinsMuk[p] = (*matrix)[p] - muk[p];
			}
			// W_k += cik * xiMoinsMuk * xiMoinsMuk'
			_tabWk[k]->add(xiMoinsMuk, cik); // virtual method

			matrix++;
			p_cik++;
		}
		p_tabMean++;
		(*_W) += _tabWk[k];
	} //end for k
}

void GaussianParameter::initFreeProportion(ModelType * iModelType) {

	if (hasFreeProportion(iModelType->_nameModel)) {
		_freeProportion = true;
	}
	else {
		_freeProportion = false;
	}
}

//----------------------------------------------
// compute class assigment of idxSample element (idxSample : 0->_nbSample-1)
// Note : distance euclidienne 
//----------------------------------------------
int64_t GaussianParameter::computeClassAssigment(int64_t idxSample) {
	GaussianData * data = _model->getGaussianData();

	int64_t p, k, k0 = 0;
	double bestDist = 0.0;
	double * x_idxSample = (data->getYStore())[idxSample];

	double dist, tmp;

	double ** p_tabMean = _tabMean; // pointeur pour parcourir _tabMean
	double * p_tabMean_k;

	for (k = 0; k < _nbCluster; k++) {
		p_tabMean_k = *(p_tabMean);
		dist = 0.0;

		// calcul
		for (p = 0; p < _pbDimension; p++) {
			tmp = x_idxSample[p] - p_tabMean_k[p];
			dist += tmp * tmp;
		}

		// recherche de la plus grande distance et de l'index k0 correspondant
		if (dist > bestDist) {
			bestDist = dist;
			k0 = k;
		}

		p_tabMean++; // next mean vector
	}

	return k0;
}

/****************************************************/
/*  computeTabMean  in USER_PARTITION initialization*/
/***************************************************/
void GaussianParameter::computeTabMeanInitUSER_PARTITION(int64_t & nbInitializedCluster, 
		bool * tabNotInitializedCluster, Partition * initPartition) 
{
	int64_t k;
	int64_t i;
	int64_t ** initPartitionValue = initPartition->_tabValue;
	int64_t nbSample = _model->getNbSample();
	GaussianData * data = _model->getGaussianData();
	double ** matrix = data->getYStore();
	double cik = 0.0;
	double * tabWeight = data->_weight;
	//double * tabWeightK = new double[_nbCluster];
	std::unique_ptr<double[]> tabWeightK(new double[_nbCluster]);    
	int64_t p; // parcours 0,..., _pbDimension

	for (k = 0; k < _nbCluster; k++) {
		tabWeightK[k] = 0;
		for (p = 0; p < _pbDimension; p++) {
			_tabMean[k][p] = 0.0;
		}
		for (i = 0; i < nbSample; i++) {
			// calcul de cik * wi
			if (initPartitionValue[i][k] == 1) {
				cik = tabWeight[i];
				tabWeightK[k] += cik;
				// ajout de la valeur de l'individu courant
				for (p = 0; p < _pbDimension; p++) {
					_tabMean[k][p] += matrix[i][p] * cik;
				}
			}
		} //end for i

		if (tabWeightK[k] != 0) {
			//---------------
			// redivision par le poids de la classe k
			for (p = 0; p < _pbDimension; p++) {
				_tabMean[k][p] /= tabWeightK[k];
			}
		}
		else {
			// initialisation RANDOM for this cluster
			//----------------------------------------
		}
	}

	nbInitializedCluster = 0;
	for (int64_t k = 0; k < _nbCluster; k++) {
		if (tabWeightK[k] == 0) {
			tabNotInitializedCluster[k] = true;
			//cout<<"classe vide � l'initialisation"<<endl;
		}
		else {
			tabNotInitializedCluster[k] = false;
			nbInitializedCluster++;
		}
	}

	//delete [] tabWeightK;
}

void GaussianParameter::computeTabMean() {

	int64_t k;
	int64_t i;
	double ** tabCik = _model->getTabCik();
	double * tabNk = _model->getTabNk();
	int64_t nbSample = _model->getNbSample();
	GaussianData * data = _model->getGaussianData();
	double ** matrix = data->getYStore();
	double * muk;
	double cik = 0.0;
	// pointeurs
	double ** p_matrix;
	double ** p_cik;
	double ** p_tabMean = _tabMean;
	double * p_weight;
	int64_t p; // parcours 0,..., _pbDimension

	for (k = 0; k < _nbCluster; k++) {
		// refresh
		muk = *p_tabMean;
		for (p = 0; p < _pbDimension; p++)
			muk[p] = 0.0;

		// initialisations
		p_matrix = matrix;
		p_cik = tabCik;
		p_weight = data->_weight;

		for (i = 0; i < nbSample; i++) {
			// calcul de cik * wi
			cik = (*p_cik)[k] * (*p_weight);


			// ajout de la valeur de l'individu courant
			for (p = 0; p < _pbDimension; p++) {
				muk[p] += (*p_matrix)[p] * cik;
			}

			p_matrix++; // individu suivant
			p_cik++; // ligne suivante de Tik
			p_weight++;
		}

		// redivision par le poids total
		for (p = 0; p < _pbDimension; p++) {
			muk[p] /= tabNk[k];
		}

		p_tabMean++;
	}
}

/// compute Mean when there is only one cluster
/// called by initRANDOM, getLogLikelihoodOne
void GaussianParameter::computeMeanOne(double * Mean, double * weight, 
		double** y_Store, int64_t nbSample, double totalWeight) const 
{
	double ** p_yStore = y_Store;
	double * p_weight = weight;
	double * p_yStore_i;
	int64_t p;
	int64_t i;
	double wi;

	initToZero(Mean, _pbDimension);

	for (i = 0; i < nbSample; i++) {
		p_yStore_i = *p_yStore;
		wi = *p_weight;
		for (p = 0; p < _pbDimension; p++) {
			Mean[p] += p_yStore_i[p] * wi;
		}
		p_weight++;
		p_yStore++;
	}
	for (p = 0; p < _pbDimension; p++) {
		Mean[p] /= totalWeight;
	}
}


/*------------------------*/
/*------------------------*/
/* Initialization Methods */
/*------------------------*/
/*------------------------*/


//------------
//------------
// Algorithms
//------------
//------------

/*-------------------------------------------------------------------------------------------
	M step
	------
	
	already updated :
	- _model (_tabFik, _tabTik, _tabCik, _tabNk) if Estep is done before
	- _model->_tabCik, _model->_tabNk if USER_PARTITION is done before (Disciminant analysis)
	In all cases, only  _tabCik and _tabNk are needed 
	
	updated in this method :
	- _tabProportion
	- _tabMean
	- _tabWk 
	- _W
--------------------------------------------------------------------------------------------*/
void GaussianParameter::MStep() {

	/* Proportion estimator (if proportions free) */
	computeTabProportion();

	/* Centers mean estimator */
	computeTabMean();

	/* Compute Wk W */
	computeTabWkW();
}

double GaussianParameter::determinantDiag(double * mat_store, Exception& errorType) {
	int64_t p;
	double det = mat_store[0];
	for (p = 1; p < _pbDimension; p++) {
		det *= mat_store[p];
	}
	if (det < minDeterminantValue)
		throw errorType;
	return det;
}

//-------------
// updateForCV
//------------
//updating the gaussian parameter from the original Model 
//without all samples belonging to the CVblock
void GaussianParameter::updateForCV(Model * originalModel, CVBlock & CVBlock) {
	GaussianParameter * oParam = (originalModel->getGaussianParameter());
	Matrix ** oTabWk = oParam->getTabWk();
	double ** oTabMean = oParam->getTabMean();
	double * oTabNk = originalModel->getTabNk();

	int64_t k;
	GaussianData * oData = originalModel->getGaussianData();
	double ** matrix = oData->getYStore();
	double * x_i;

	double * tabNk = _model->getTabNk();

	double ** tabCik = _model->getTabCik();

	//--------------------------------- updates the proportions
	computeTabProportion();

	//--------------------------------- updates the means
	//  mu*_k = (1/n*_k).(n_k.mu_k - sum_{test}(cik wi xi))

	int64_t p;
	int64_t i;
	double cik;

	for (k = 0; k < _nbCluster; k++) {
		for (p = 0; p < _pbDimension; p++) {
			_tabMean[k][p] = oTabMean[k][p] * oTabNk[k];
		}
		for (int64_t ii = 0; ii < CVBlock._nbSample; ii++) {
			i = CVBlock._tabWeightedIndividual[ii].val;
			cik = tabCik[i][k] * CVBlock._tabWeightedIndividual[ii].weight;
			x_i = matrix[i];
			for (p = 0; p < _pbDimension; p++) {
				_tabMean[k][p] -= cik * x_i[p];

			}
		}
		for (p = 0; p < _pbDimension; p++) {
			_tabMean[k][p] /= tabNk[k];
		}
	}

	//--------------------------------- updates the W_k

	/*
	_tabWk*[k] = _tabWk[k] - Sum(i in CVBlock)cik*(xi-mu*k)(xi-mu*k)' + nk*(muk-oMuk)(muk-oMuk)'
	 */

	double * tmp = oData->getTmpTabOfSizePbDimension();
	//double * mukMoinsoMuk = new double[_pbDimension];
	std::unique_ptr<double[]> mukMoinsoMuk(new double[_pbDimension]);    

	double * xi;
	*(_W) = 0.0;
	for (k = 0; k < _nbCluster; k++) {
		// _tabWk[k]->recopy(oTabWk[k]);
		(* _tabWk[k]) = oTabWk[k];
		for (int64_t ii = 0; ii < CVBlock._nbSample; ii++) {
			i = CVBlock._tabWeightedIndividual[ii].val;
			xi = matrix[i];
			for (p = 0; p < _pbDimension; p++) {
				tmp[p] = xi[p] - _tabMean[k][p];
			}
			cik = -tabCik[i][k] * CVBlock._tabWeightedIndividual[ii].weight;
			_tabWk[k]->add(tmp, cik);
		}

		for (p = 0; p < _pbDimension; p++) {
			mukMoinsoMuk[p] = _tabMean[k][p] - oTabMean[k][p];
		}
		_tabWk[k]->add(mukMoinsoMuk.get(), oTabNk[k]); // change : 17/03/06
		(*_W) += _tabWk[k];
	}
	//delete [] mukMoinsoMuk;
}

}
