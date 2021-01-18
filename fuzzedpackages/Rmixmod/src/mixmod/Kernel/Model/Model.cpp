/***************************************************************************
                             SRC/mixmod/Kernel/Model/Model.cpp  description
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

#include "mixmod/Kernel/Model/Model.h"
#include "mixmod/Kernel/Model/ModelType.h"
#include "mixmod/Kernel/Parameter/GaussianDiagParameter.h"
#include "mixmod/Kernel/Parameter/GaussianSphericalParameter.h"
#include "mixmod/Kernel/Parameter/GaussianGeneralParameter.h"
#include "mixmod/Kernel/Parameter/GaussianHDDAParameter.h"
#include "mixmod/Kernel/Parameter/BinaryEParameter.h"
#include "mixmod/Kernel/Parameter/BinaryEkParameter.h"
#include "mixmod/Kernel/Parameter/BinaryEjParameter.h"
#include "mixmod/Kernel/Parameter/BinaryEkjParameter.h"
#include "mixmod/Kernel/Parameter/BinaryEkjhParameter.h"
#include "mixmod/Kernel/Parameter/CompositeParameter.h"
#include "mixmod/Utilities/Random.h"
#include "mixmod/Kernel/IO/GaussianData.h"
#include "mixmod/Kernel/IO/BinaryData.h"
#include "mixmod/Kernel/IO/CompositeData.h"
#include "mixmod/Kernel/IO/Partition.h"
#include "mixmod/Kernel/IO/LabelDescription.h"
#include "mixmod/Matrix/Matrix.h"
#include "mixmod/Kernel/IO/Sample.h"
#include "mixmod/Kernel/IO/BinarySample.h"
#include <string.h>

namespace XEM {

//------------
// Constructor
//------------
Model::Model() {
	THROW(OtherException, internalMixmodError);
}

Model* Model::clone() {
	return new Model(this);
}

//------------
// Constructor
//------------
Model::Model(Model * iModel)
: _modelType(iModel->getModelType()), _nbCluster(iModel->getNbCluster()),
  _nbSample(iModel->getNbSample()), _deleteData(true), _parameter(iModel->getParameter()->clone()),
  _tabFik(copyTab(iModel->getTabFik(), _nbSample, _nbCluster)),
  _tabSumF(copyTab(iModel->getTabSumF(), _nbSample)), _tabTik(copyTab(iModel->getTabTik(),
  _nbSample, _nbCluster)), _tabZikKnown(copyTab(iModel->getTabZikKnown(), _nbSample,
  _nbCluster)), _tabCik(copyTab(iModel->getTabCik(), _nbSample, _nbCluster)),
  _tabZiKnown(copyTab(iModel->getTabZiKnown(), _nbSample)), _tabNk(copyTab(iModel->getTabNk(),
  _nbCluster)), _algoName(iModel->getAlgoName())
{
	if (isHeterogeneous(_modelType->_nameModel)) {
		CompositeData * cD = (CompositeData *) iModel->getData();
		_data = new CompositeData(cD);
	}
	else {
		if (isBinary(_modelType->_nameModel)) {
			BinaryData * bD = (iModel->getBinaryData());
			_data = new BinaryData(*bD);
		}
		else {
			GaussianData * gD = iModel->getGaussianData();
			_data = new GaussianData(*gD);
		}
	}
	// set model into parameter object
	_parameter->setModel(this);
}

//------------
// Constructor
//------------
Model::Model(ModelType * modelType, int64_t nbCluster, Data *& data, Partition * knownPartition)
: _modelType(modelType), _nbCluster(nbCluster), _nbSample(data->_nbSample),
  _data(data), _deleteData(false), _parameter(0), _algoName(UNKNOWN_ALGO_NAME)
{
	// initialize probabilities
	int64_t k, i;

	_tabFik = new double*[_nbSample];
	_tabCik = new double*[_nbSample];
	_tabSumF = new double[_nbSample];
	_tabTik = new double*[_nbSample];
	_tabZikKnown = new int64_t *[_nbSample];
	_tabZiKnown = new bool[_nbSample];
	_tabNk = new double[_nbCluster];

	for (i = 0; i < _nbSample; i++) {
		_tabFik[i] = new double[_nbCluster];
		_tabTik[i] = new double[_nbCluster];
		_tabZikKnown[i] = new int64_t [_nbCluster];
		_tabCik[i] = new double[_nbCluster];
		for (k = 0; k < _nbCluster; k++) {
			_tabFik[i][k] = 0.0;
			_tabTik[i][k] = 0.0;
			_tabZikKnown[i][k] = 0;
			_tabCik[i][k] = 0.0;
		}
		_tabZiKnown[i] = false;
		_tabSumF[i] = 0.0;
	}

	// _tabNk[k] = 0 even if knownPartition because this partition could be partial
	for (k = 0; k < _nbCluster; k++) {
		_tabNk[k] = 0.0;
	}

	// save the partition
	FixKnownPartition(knownPartition);

	// set the parameters
	ModelName modelName = _modelType->_nameModel;
	// create Param
	if (isSpherical(modelName)) {
		_parameter = new GaussianSphericalParameter(this, _modelType);
	}
	else if (isDiagonal(modelName)) {
		_parameter = new GaussianDiagParameter(this, _modelType);
	}
	else if (isGeneral(modelName)) {
		_parameter = new GaussianGeneralParameter(this, _modelType);
	}
		//HDDA models
	else if (isHD(modelName)) {
		_parameter = new GaussianHDDAParameter(this, _modelType);
	}
	// Binary models
	else {
		if (modelName == Binary_p_E) {
			_parameter = new BinaryEParameter(
					this, _modelType, ((BinaryData*) _data)->getTabNbModality());
		}
		else if (modelName == Binary_p_Ek) {
			_parameter = new BinaryEkParameter(
					this, _modelType, ((BinaryData*) _data)->getTabNbModality());
		}
		else if (modelName == Binary_p_Ej) {
			_parameter = new BinaryEjParameter(
					this, _modelType, ((BinaryData*) _data)->getTabNbModality());
		}
		else if (modelName == Binary_p_Ekj) {
			_parameter = new BinaryEkjParameter(
					this, _modelType, ((BinaryData*) _data)->getTabNbModality());
		}
		else if (modelName == Binary_p_Ekjh) {
			_parameter = new BinaryEkjhParameter(
					this, _modelType, ((BinaryData*) _data)->getTabNbModality());
		}
		else if (modelName == Binary_pk_E) {
			_parameter = new BinaryEParameter(
					this, _modelType, ((BinaryData*) _data)->getTabNbModality());
		}
		else if (modelName == Binary_pk_Ek) {
			_parameter = new BinaryEkParameter(
					this, _modelType, ((BinaryData*) _data)->getTabNbModality());
		}
		else if (modelName == Binary_pk_Ej) {
			_parameter = new BinaryEjParameter(
					this, _modelType, ((BinaryData*) _data)->getTabNbModality());
		}
		else if (modelName == Binary_pk_Ekj) {
			_parameter = new BinaryEkjParameter(
					this, _modelType, ((BinaryData*) _data)->getTabNbModality());
		}
		else if (modelName == Binary_pk_Ekjh) {
			_parameter = new BinaryEkjhParameter(
					this, _modelType, ((BinaryData*) _data)->getTabNbModality());
		}
		else if (isHeterogeneous(modelName)) {
			_parameter = new CompositeParameter(
					this, _modelType, (_data->getBinaryData())->getTabNbModality());
		}
	}
}

//-----------
// Destructor
//-----------
Model::~Model() {
	int64_t i;

	if (_tabFik) {
		for (i = 0; i < _nbSample; i++) {
			delete[] _tabFik[i];
			_tabFik[i] = NULL;
		}
		delete[] _tabFik;
		_tabFik = NULL;
	}

	if (_tabCik) {
		for (i = 0; i < _nbSample; i++) {
			delete[] _tabCik[i];
			_tabCik[i] = NULL;
		}
		delete[] _tabCik;
		_tabCik = NULL;
	}

	if (_tabTik) {
		for (i = 0; i < _nbSample; i++) {
			delete[] _tabTik[i];
			_tabTik[i] = NULL;
		}
		delete[] _tabTik;
		_tabTik = NULL;
	}

	if (_tabZikKnown) {
		for (i = 0; i < _nbSample; i++) {
			delete[] _tabZikKnown[i];
			_tabZikKnown[i] = NULL;
		}
		delete[] _tabZikKnown;
		_tabZikKnown = NULL;
	}

	if (_tabZiKnown) {
		delete[] _tabZiKnown;
		_tabZiKnown = NULL;
	}

	if (_tabNk) {
		delete[] _tabNk;
		_tabNk = NULL;
	}

	if (_tabSumF) {
		delete[] _tabSumF;
		_tabSumF = NULL;
	}

	if (_parameter) {
		delete _parameter;
		_parameter = NULL;
	}

	if (_deleteData) {
		delete _data;
		_data = NULL;
	}
}

//------------
// get the log of N
//------------
double Model::getLogN() {
	return log(_data->_weightTotal);
}

//------------
// updateForCV
//------------
void Model::updateForCV(Model * originalModel, CVBlock & CVBlock) {
	int64_t k, i;

	// _data
	//------
	// _data=originalData but weight will be different
	// _data has already been updated (copy constructor)
	_data->_weightTotal = originalModel->_data->_weightTotal - CVBlock._weightTotal;
	recopyTab(originalModel->_data->_weight, _data->_weight, _nbSample);
	for (int64_t ii = 0; ii < CVBlock._nbSample; ii++) {
		i = CVBlock._tabWeightedIndividual[ii].val;
		_data->_weight[i] -= CVBlock._tabWeightedIndividual[ii].weight;
	}

	/* new version : Dec 2006*/
	/* recopy fik, tik and zik*/
	recopyTab(originalModel->_tabFik, _tabFik, _nbSample, _nbCluster);
	recopyTab(originalModel->_tabSumF, _tabSumF, _nbSample);
	recopyTab(originalModel->_tabTik, _tabTik, _nbSample, _nbCluster);
	recopyTab(originalModel->_tabCik, _tabCik, _nbSample, _nbCluster);

	// already done : (with the copy constructor
	//  recopyTab(originalModel->_tabZikKnown, _tabZik, _nbSample);

	//--------------
	// update _tabNk
	//--------------
	recopyTab(originalModel->_tabNk, _tabNk, _nbCluster);
	for (int64_t ii = 0; ii < CVBlock._nbSample; ii++) {
		i = CVBlock._tabWeightedIndividual[ii].val;
		for (k = 0; k < _nbCluster; k++) {
			_tabNk[k] -= CVBlock._tabWeightedIndividual[ii].weight * _tabCik[i][k];
		}
	}
	//----------
	// parameter
	//----------
	_parameter->updateForCV(originalModel, CVBlock);
}

//--------------
// getKnownLabel
//--------------
// return a value in : 0, ..., K-1
// i = 0 ... nbSample-1
int64_t Model::getKnownLabel(int64_t i) {
	int64_t res = -1;
	int64_t k;
	if (_tabZiKnown[i]) {
		for (k = 0; k < _nbCluster; k++) {
			if (_tabZikKnown[i][k] == 1) {
				res = k;
			}
		}
	}
	else {
		THROW(OtherException, internalMixmodError);
	}
	return res;
}

//-------------------------------------------
// getLabelAndPartitionByMAPOrKnownPartition
// label[i] = 1 ... nbSample
//------------------------------------------
void Model::getLabelAndPartitionByMAPOrKnownPartition(int64_t * label, int64_t ** partition) {
	//if (!_isCikEqualToMapOrKnownPartition){
	if (_algoName == UNKNOWN_ALGO_NAME)
		throw;

	if (_algoName == MAP || _algoName == CEM || _algoName == M) {
		// _tabCik contains the result
		int64_t i;
		int64_t k;
		for (i = 0; i < _nbSample; i++) {
			for (k = 0; k < _nbCluster; k++) {
				// cast double to int64_t  (_tabCik[i][k] = 0.0 or 1.0)
				partition[i][k] = (int64_t) _tabCik[i][k];
				if (partition[i][k] == 1) {
					label[i] = k + 1;
				}
			}
		}
	}
	else {
		int64_t k, kMax;
		int64_t i;
		double tikMax = 0;
		for (i = 0; i < _nbSample; i++) {
			if (_tabZiKnown[i]) {
				//-----------------
				for (k = 0; k < _nbCluster; k++) {
					partition[i][k] = _tabZikKnown[i][k];
					if (_tabZikKnown[i][k] == 1) {
						label[i] = k + 1;
					}
				}
			}
			else {
				// !ziKnown
				//---------
				kMax = 0;
				tikMax = _tabTik[i][0];
				for (k = 1; k < _nbCluster; k++) {
					if (_tabTik[i][k] > tikMax) {
						tikMax = _tabTik[i][k];
						kMax = k;
					}
				}
				for (k = 0; k < _nbCluster; k++) {
					partition[i][k] = 0;
					partition[i][kMax] = 1;
					label[i] = kMax + 1;
				}
			}
		}
	}
}

//------------------------------
// getLabelByMAPOrKnownPartition
//------------------------------
// return a value in : 0, ..., K-1
// i = 0 ... nbSample-1
int64_t Model::getLabelByMAPOrKnownPartition(int64_t i) {

	int64_t k, kMax;
	double tikMax = 0;
	int64_t res = -1;

	if (_algoName == UNKNOWN_ALGO_NAME)
		throw;
	//
	if (_algoName == CEM || _algoName == MAP || _algoName == M) {
		//_tabCik[i] gives to result
		for (k = 0; k < _nbCluster; k++) {
			if (_tabCik[i][k] == 1) {
				res = k;
			}
		}
	}

	else {
		// must be computed
		// This method is called by getCompletedLogLikelihood
		// In this case (!_isCikEqualToMapOrKnownPartition), it's called
		// by ICLCriterion or LikelihoodOutput, so an Estep have been done before
		// to update fik and tik used in this section
		if (_tabZiKnown[i]) {
			for (k = 0; k < _nbCluster; k++) {
				if (_tabZikKnown[i][k] == 1) {
					res = k;
				}
			}
		}
		else {
			tikMax = _tabTik[i][0];
			kMax = 0;
			for (k = 0; k < _nbCluster; k++) {
				if (_tabTik[i][k] > tikMax) {
					tikMax = _tabTik[i][k];
					kMax = k;
				}
			}
			res = kMax;
		}
	}

	if (res == -1) {

		if (VERBOSE == 1)
			// label couldn't be found
			cout << "internalMixmodError in Model::getLabelByMAPOrKnownPartition, i=" << i << endl;

		THROW(OtherException, internalMixmodError);
	}

	return res;
}

//-------------------------------------
// get log-likelihood one (one cluster)
//-------------------------------------
double Model::getLogLikelihoodOne() {
	return _parameter->getLogLikelihoodOne();
}

//-------------------
// get log-likelihood
//-------------------
double Model::getLogLikelihood(bool fikMustBeComputed) {
	// Compute the log-likelihood (observed) //

	if (fikMustBeComputed) {
		computeFik();
	}

	int64_t i;
	double logLikelihood = 0.0;
	double ** p_tabFik;
	double * p_tabFik_i;
	double * weight = _data->_weight;
	p_tabFik = _tabFik;
	for (i = 0; i < _nbSample; i++) {
		p_tabFik_i = *p_tabFik;
		if (_tabZiKnown[i]) {
			int64_t ki = getKnownLabel(i); // la classe de l'individu i
			logLikelihood += log(p_tabFik_i[ki]) * weight[i];
		}
		else {
			if (_tabSumF[i] > 0)
				logLikelihood += log(_tabSumF[i]) * weight[i];
		}
		//cout<<"compute LL, with ind "<<i<<", LL = "<<logLikelihood<<endl;
		p_tabFik++;
	}

	return logLikelihood;
}

//-----------------------------------------
// get completed LL (if CEM) or LL (elseif)
//-----------------------------------------
double Model::getCompletedLogLikelihoodOrLogLikelihood() {
	if (_algoName == UNKNOWN_ALGO_NAME) {
		THROW(OtherException, internalMixmodError);
	}
	else {
		if (_algoName == CEM) {
			return getCompletedLogLikelihood();
		}
		else {
			return getLogLikelihood(true);
		}
	}
}

//-----------------------------
// get completed log-likelihood
//-----------------------------
double Model::getCompletedLogLikelihood() {
	// Compute the observed completed log-likelihood

	int64_t i;
	double cLogLikelihood = 0.0;
	int64_t ki; // la classe de l'individu i
	for (i = 0; i < _nbSample; i++) {
		ki = getLabelByMAPOrKnownPartition(i); // la classe de l'individu i
		if (_tabFik[i][ki] > 0) {
			cLogLikelihood += log(_tabFik[i][ki]) * _data->_weight[i];
		}
	}

	return cLogLikelihood;
}

//------------
// get entropy
//------------
double Model::getEntropy() {
	// Entropy

	// Initialization //
	int64_t i;
	int64_t k;
	double entropy = 0;

	// Compute entropy: sum[tik * log(tik)] //
	for (i = 0; i < _nbSample; i++) {
		// ajout du 16/06/2004 : ligne suivante : on ajoute uniqt si le label est inconnu
    if (!_tabZiKnown[i]) {
      for (k = 0; k < _nbCluster; k++) {
        if (_tabTik[i][k] > 0 && _tabTik[i][k] != 1) {
          entropy += _tabTik[i][k] * log(_tabTik[i][k]) * _data->_weight[i];
        }
      }
    }
  }

	return -entropy;
}

//------------
// get entropy matrix (for massiccc's visualization)
//------------
vector< vector<double> > Model::getEntropyMatrix() {

  // Initialization //
  int nbVariables = _parameter->getModel()->getData()->getPbDimension();
  vector< vector<double> > entropyM(nbVariables, vector<double>(_nbCluster));
  vector< vector<double> > tabTikj(_nbSample, vector<double>(_nbCluster));
  double sum = 0.0;
  double nlnk = _nbSample * log(_nbCluster);

  for (int64_t j = 0; j < nbVariables; j++) {
    if (isHeterogeneous(_modelType->_nameModel)) {
      CompositeParameter * cParameter = dynamic_cast<CompositeParameter*> (_parameter);
      //Binary Parameter
      BinaryParameter * bParameter = cParameter->getBinaryParameter();
      int64_t bPbDimension = bParameter->getPbDimension();
      if (j < bPbDimension) {
        double *** scatter = bParameter->scatterToArray();
        for (int64_t i = 0; i < _nbSample; i++) {
          for (int64_t k = 0; k < _nbCluster; k++) {
            int64_t tabNbModality = bParameter->getTabNbModality()[j];
            int64_t tabDataValue = bParameter->getModel()->getBinaryData()->getDataMatrix()[i]->getBinarySample()->getTabValue()[j];
            double prod = 1.0;
            for (int64_t m = 0; m < tabNbModality; m++) {
              if (tabDataValue == (m + 1)) {
                if (bParameter->getTabCenter()[k][j] == (m + 1)) {
                  prod = prod * ((1 - scatter[k][j][m]));
                }
                else {
                  prod = prod * (scatter[k][j][m] / (tabNbModality - 1));
                }
              }
            }
            tabTikj[i][k] = bParameter->getTabProportion()[k] * prod;
            sum += tabTikj[i][k];
          }
          for (int64_t k = 0; k < _nbCluster; k++) {
            tabTikj[i][k] = tabTikj[i][k] / sum;
          }
          sum = 0.0;
        }
        for (int64_t l = 0; l < _nbCluster; l++) {
          for (int64_t m = 0; m < bPbDimension; m++) {
            delete[] scatter[l][m];
          }
          delete[] scatter[l];
        }
        delete[] scatter;
      }
      else {
        //Gaussian Parameter
        GaussianEDDAParameter * gParameter = dynamic_cast<GaussianEDDAParameter*> (cParameter->getGaussianParameter());
        int gPbDimension = gParameter->getPbDimension();
        double ** tabData = gParameter->getModel()->getGaussianData()->getYStore();
        for (int64_t i = 0; i < _nbSample; i++) {
          for (int64_t k = 0; k < _nbCluster; k++) {
            double ** store = gParameter->getTabSigma()[k]->storeToArray();
            tabTikj[i][k] = gParameter->getTabProportion()[k] * 1.0/(sqrt(2.0*XEMPI)*sqrt(store[j - bPbDimension][j - bPbDimension])) * exp(-pow(tabData[i][j - bPbDimension] - gParameter->getTabMean()[k][j - bPbDimension],2)/(2*store[j - bPbDimension][j - bPbDimension]));
            sum += tabTikj[i][k];
            for (int64_t l = 0; l < gPbDimension; l++) {
              delete[] store[l];
            }
            delete[] store;
          }
          for (int k = 0; k < _nbCluster; k++) {
            tabTikj[i][k] = tabTikj[i][k] / sum;
          }
          sum = 0.0;
        }
      }
    }
    else {
      if (isBinary(_modelType->_nameModel)) {
        BinaryParameter * bParameter = dynamic_cast<BinaryParameter*> (_parameter);
        double *** scatter = bParameter->scatterToArray();
        for (int64_t i = 0; i < _nbSample; i++) {
          int tabNbModality = bParameter->getTabNbModality()[j];
          int tabDataValue = bParameter->getModel()->getBinaryData()->getDataMatrix()[i]->getBinarySample()->getTabValue()[j];
          for (int64_t k = 0; k < _nbCluster; k++) {
            double prod = 1.0;
            for (int64_t m = 0; m < tabNbModality; m++) {
              if (tabDataValue == (m + 1)) {
                if (bParameter->getTabCenter()[k][j] == (m + 1)) {
                  prod = prod * ((1 - scatter[k][j][m]));
                }
                else { 
                  prod = prod * (scatter[k][j][m] / (tabNbModality - 1));
                }
              }
            }
            tabTikj[i][k] = bParameter->getTabProportion()[k] * prod;
            sum += tabTikj[i][k];
          }
          for (int64_t k = 0; k < _nbCluster; k++) {
            tabTikj[i][k] = tabTikj[i][k] / sum;
          }
          sum = 0.0;
        }
        for (int64_t l = 0; l < _nbCluster; l++) {
          for (int64_t m = 0; m < nbVariables; m++) {
            delete[] scatter[l][m];
          }
          delete[] scatter[l];
        }
        delete[] scatter;

      }
      else {
        GaussianEDDAParameter * gParameter = dynamic_cast<GaussianEDDAParameter*> (_parameter);
        double ** tabData = gParameter->getModel()->getGaussianData()->getYStore();
        for (int64_t i = 0; i < _nbSample; i++) {
          for (int64_t k = 0; k < _nbCluster; k++) {
            double ** store = gParameter->getTabSigma()[k]->storeToArray();
            tabTikj[i][k] = gParameter->getTabProportion()[k] * (1.0/(sqrt(2.0*XEMPI*store[j][j])) * exp(-pow(tabData[i][j] - gParameter->getTabMean()[k][j],2)/(2*store[j][j])));
            sum += tabTikj[i][k];
            for (int64_t l = 0; l < nbVariables; l++) {
              delete[] store[l];
            }
            delete[] store;
          }
          for (int64_t k = 0; k < _nbCluster; k++) {
            tabTikj[i][k] = tabTikj[i][k] / sum;
          }
          sum = 0.0;
        }

      }
    }

    for (int64_t k = 0; k < _nbCluster; k++) {
      for (int64_t i = 0; i < _nbSample; i++) {
        if (tabTikj[i][k] == 0.0) tabTikj[i][k] = 1.0;
        entropyM[j][k] += - tabTikj[i][k] * log(tabTikj[i][k]);
      }
      entropyM[j][k] = entropyM[j][k] / nlnk;
    }
  }
  return entropyM;
}

// set parameters
void Model::setParameter(Parameter * parameter) {
	if (_parameter) delete _parameter;
	_parameter = parameter;
}

// set name of the algorithm
void Model::setAlgoName(AlgoName algoName) {
	_algoName = algoName;
}

AlgoName Model::getAlgoName() {
	return _algoName;
}

// set an error for the model
void Model::setError(Exception& errorType) {
	_error.setError(errorType);
}

//----------------------------------
// compute the probabilities _tabFik
//----------------------------------
void Model::computeFik() {
	// updates _tabSumF, _tabFik

	int64_t k;
	int64_t i;
	double * tabProportion = _parameter->getTabProportion();
	double * p_tabSumF = _tabSumF; // parcours de _tabSum
	double ** p_tabFik = _tabFik; // parcours de _tabFik

	_parameter->getAllPdf(_tabFik, tabProportion);

	for (i = 0; i < _nbSample; i++) {
		*p_tabSumF = 0.0;
		for (k = 0; k < _nbCluster; k++) {
			*p_tabSumF += (*p_tabFik)[k];
		}
		p_tabFik++;
		p_tabSumF++;
	}
}

//-------------------------------------
// compute number of element by cluster
//-------------------------------------
void Model::computeNk() {
	int64_t k;
	int64_t i;
	double ** p_tabCik = _tabCik; // parcours le tableau _tabCik
	double * w = _data->_weight; // parcours des poids
	double wi; // poids courant
	double * p_tabCik_i;

	// initialisation
	initToZero(_tabNk, _nbCluster);
	for (i = 0; i < _nbSample; i++) {
		wi = *w;
		p_tabCik_i = *p_tabCik;
		for (k = 0; k < _nbCluster; k++) {
			_tabNk[k] += p_tabCik_i[k] * wi;
		}

		p_tabCik++;
		w++;
	}
	// verification
	for (k = 0; k < _nbCluster; k++) {
		if (_tabNk[k] == 0.0) {
			THROW(NumericException, nullNk);
		}
	}
}

//----------------
//getFreeParameter
//----------------
int64_t Model::getFreeParameter() {
	return _parameter->getFreeParameter();
}

// compute label of samples in x  (res : 0 -> _nbCluster-1)
//------------------------------
int64_t Model::computeLabel(Sample * x) {
	int64_t k, res = 0;
	double sumfk = 0.0;
	double * tk = new double[_nbCluster];
	double * fk = new double[_nbCluster];
	double max = 0.0;
	double * tabProportion = _parameter->getTabProportion();
	double tmp;

	// Compute the probabilities _tabFik
	for (k = 0; k < _nbCluster; k++) {
	  tmp = tabProportion[k] * _parameter->getPdf(x, k);
		fk[k] = tmp;
		// Compute the sum of probabilities sumfk
		sumfk += tmp;
	}

	// Compute the conditional probabilities (posteriori probabilities)
	for (k = 0; k < _nbCluster; k++) {
		tk[k] = fk[k] / sumfk;
	}

	// get the max of tk[k]
	for (k = 0; k < _nbCluster; k++) {
		tmp = tk[k];
		if (tmp > max) {
			max = tmp;
			res = k;
		}
	}

	delete[] fk;
	delete[] tk;

	return res;
}

//-------------
// computeLabel
//-------------
// compute the label of the i0-th point of the sample
// i0 : 0 -> _nbSample-1  and res : 0 -> _nbCluster-1)
int64_t Model::computeLabel(int64_t i0) {
	int64_t k, res = 0;
	double max = 0.0;
	double * Tik_i0 = _tabTik[i0];

	for (k = 0; k < _nbCluster; k++) {
		if (Tik_i0[k] > max) {
			max = Tik_i0[k];
			res = k;
		}
	}

	return res;
}

//-----------------
// Fix label Known
//----------------
void Model::FixKnownPartition(Partition *& knownPartition) {
	// update Nk if knownLabel
	if (knownPartition != NULL) {
		int64_t ** knownPartitionValue = knownPartition->_tabValue;
		int64_t * knownPartition_i;
		double ** p_cik = _tabCik;
		int64_t ** p_zikKnown = _tabZikKnown;
		//double ** p_tik = _tabTik;
		int64_t k;
		int64_t i;
		double sumLabel = 0.0;
		for (i = 0; i < _nbSample; i++) {
			sumLabel = 0.0;
			knownPartition_i = *knownPartitionValue;

			for (k = 0; k < _nbCluster; k++) {
				sumLabel += knownPartition_i[k];
			}
			if (sumLabel != 0.0) {
				_tabZiKnown[i] = true;
				recopyTab(knownPartition_i, *p_cik, _nbCluster);
				recopyTab(knownPartition_i, *p_zikKnown, _nbCluster);
				//recopyTab(knownLabel_i,*p_tik,_nbCluster);
			}
			knownPartitionValue++;
			p_cik++;
			//p_tik++;
			p_zikKnown++;
		}

		computeNk();
	}//endif
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
// 										Initialization methods
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------

/*-------------------------------------------
		initRANDOM
		----------

  updated in this method :
	- _parameter
	- _tabFik and _tabSumF (because bestParameter is choose
      with the best LL which is computed with fik (and sumF)
	Note : _tabFik and sumF wil be 're'computed in the following EStep
	So only _parameter have to be updated in this method
-------------------------------------------*/
void Model::initRANDOM(int64_t nbTry) {
	// cout<<"init RANDOM, nbTryInInit="<<nbTry<<endl;
	_algoName = UNKNOWN_ALGO_NAME;
	int64_t i, k;
	double logLikelihood, bestLogLikelihood;
	Parameter * bestParameter = _parameter->clone();
	bool * tabIndividualCanBeUsedForInitRandom = new bool[_nbSample];
	for (i = 0; i < _nbSample; i++) {
		tabIndividualCanBeUsedForInitRandom[i] = true;
	}
	bool * tabClusterToInitialize = new bool[_nbCluster];
	for (k = 0; k < _nbCluster; k++) {
		tabClusterToInitialize[k] = true;
	}

	// 1. InitForInitRandom
	//---------------------
	_parameter->initForInitRANDOM();
	
	// 1rst RANDOM
	//-------------
	randomForInitRANDOMorUSER_PARTITION(
			tabIndividualCanBeUsedForInitRandom, tabClusterToInitialize);
	// Compute log-likelihood
	logLikelihood = getLogLikelihood(true); // true : to compute fik
	bestLogLikelihood = logLikelihood;
	bestParameter->recopy(_parameter);

	/*cout<<"initRandom"<<endl<<"1er essai : "<<endl<<"Parameter : "<<endl;
	_parameter->edit();
	cout<<"LL : "<< bestLogLikelihood<<endl;*/

	// Others RANDOM
	for (i = 1; i < nbTry; i++) {
		randomForInitRANDOMorUSER_PARTITION(
				tabIndividualCanBeUsedForInitRandom, tabClusterToInitialize);
		// Compute log-likelihood
		logLikelihood = getLogLikelihood(true); // true : to compute fik
		if (logLikelihood > bestLogLikelihood) {
			bestLogLikelihood = logLikelihood;
			bestParameter->recopy(_parameter);
		}
		/*cout<<endl<<"initRandom"<<endl<<i+1<<" eme essai : "<<endl<<"Parameter : "<<endl;
		  _parameter->edit();
		  cout<<"LL : "<< logLikelihood<<endl;*/
	}

	// set best parameter
	delete _parameter;
	_parameter = bestParameter;
	_parameter->setModel(this);
	/*cout<<endl<<"initRandom"<<endl<<"meilleur essai : "<<endl<<"Parameter : "<<endl;
	_parameter->edit();
	cout<<"LL : "<< bestLogLikelihood<<endl;*/

	//cout<<"fin de init RANDOM, nb d'essais effectues="<<i<<endl;
	delete [] tabIndividualCanBeUsedForInitRandom;
	delete [] tabClusterToInitialize;
}

//-----------------------------------------------
// random step for init RANDOM or USER_PARTITION
//----------------------------------------------
void Model::randomForInitRANDOMorUSER_PARTITION(
		bool * tabIndividualCanBeUsedForInitRandom, bool * tabClusterToInitialize)
{
	int64_t * tabIdxSampleForInit = new int64_t [_nbCluster];
	Sample ** tabSampleForInit = new Sample*[_nbCluster];
	double totalWeight = _data->_weightTotal;
	Sample ** tabSample = _data->_matrix;
	double * tabWeight = _data->_weight;
	int64_t k;

	for (k = 0; k < _nbCluster; k++) {
		if (tabClusterToInitialize[k]) {
			tabIdxSampleForInit[k] = generateRandomIndex(
					tabIndividualCanBeUsedForInitRandom, tabWeight, totalWeight);
			tabSampleForInit[k] = tabSample[tabIdxSampleForInit[k]];
		}
	} // end for k

	_parameter->updateForInitRANDOMorUSER_PARTITION(tabSampleForInit, tabClusterToInitialize);

	// update tabIndividualCanBeUsedForInitRandom for others runs
	for (k = 0; k < _nbCluster; k++) {
		if (tabClusterToInitialize[k]) {
			tabIndividualCanBeUsedForInitRandom[tabIdxSampleForInit[k]] = true;
		}
	}
	delete [] tabIdxSampleForInit;
	delete [] tabSampleForInit;
}

/*----------------------------------------
			initUSER
			--------
	updated in this method :
	- _parameter
-----------------------------------------*/
void Model::initUSER(Parameter * initParameter) {
	_algoName = UNKNOWN_ALGO_NAME;
	if (initParameter) {
		_parameter->initUSER(initParameter);
	}
	else {
		THROW(InputException, errorInitParameter);
	}
}

/*----------------------------------------
			initUSER_PARTITION
			------------------

    updated in this method :
        - _parameter

    Note : this method is only called in Classification context
           (not in Discriminant context). So, an Estep follows

-----------------------------------------*/
/*
Les partitions donnees servent a calculer
- les cik, les nk (dans fixKnownLabel) appele dans le constructeur
- les centres lorsque l'on a au moins un representant de la classe dans la initPartition
  (sinon on tire au hasard)
En revanche, on ne les utilise pas pour les dispersions.
On pourrait le faire si on a beaucoup d'information mais dans ce le cas ou l'on a peu d'information
(ex : un seul individu pour une des classes), on ne peut pas calculer de disperion.
On calcule donc la dispersion autour du centre (comme s'il y avait une seule classe)
dans le cas gaussien et on tire la dispersion au hasard dans le cas binaire.
 */
void Model::initUSER_PARTITION(Partition * initPartition, int64_t nbTryInInit) {

	_algoName = UNKNOWN_ALGO_NAME;
	int64_t nbInitializedCluster;
	bool * tabNotInitializedCluster = new bool[_nbCluster];

	// 1. InitForUSER_PARTITION
	//-------------------------
	_parameter->initForInitUSER_PARTITION(
			nbInitializedCluster, tabNotInitializedCluster, initPartition);

	// 2.init random if needed
	//------------------------
	if (nbInitializedCluster != _nbCluster) {
		// upadte tabIndividualCanBeUsedForInitRandom
		int64_t i, k;
		int64_t ** initLabelValue = initPartition->_tabValue;
		int64_t nbSampleCanBeUsedForInitRandom = _nbSample;
		bool * tabIndividualCanBeUsedForInitRandom = new bool[_nbSample];
		for (i = 0; i < _nbSample; i++) {
			tabIndividualCanBeUsedForInitRandom[i] = true;
			k = 0;
			while (k < _nbCluster && tabIndividualCanBeUsedForInitRandom[i]) {
				if (initLabelValue[i][k] == 1) {
					tabIndividualCanBeUsedForInitRandom[i] = false;
					nbSampleCanBeUsedForInitRandom--;
				}
				k++;
			}
		}
		if (nbSampleCanBeUsedForInitRandom < (_nbCluster - nbInitializedCluster)) {
			THROW(InputException,
					tooManySampleInInitPartitionAndTooManyClusterNotRepresented);
		}

		double logLikelihood, bestLogLikelihood;
		Parameter * bestParameter = _parameter->clone();

		// 1rst random
		//-------------
		//cout<<"1rst random"<<endl;
		randomForInitRANDOMorUSER_PARTITION(
				tabIndividualCanBeUsedForInitRandom, tabNotInitializedCluster);
		// Compute log-likelihood
		logLikelihood = getLogLikelihood(true); // true : to compute fik
		bestLogLikelihood = logLikelihood;
		bestParameter->recopy(_parameter);
		/*cout<<"initRandom"<<endl<<"1er essai : "<<endl<<"Parameter : "<<endl;
			_parameter->edit();
			cout<<"LL : "<< bestLogLikelihood<<endl;*/

		// Others RANDOM
		//-------------
		for (i = 1; i < nbTryInInit; i++) {
			//		cout<<i+1<<" random"<<endl;
			randomForInitRANDOMorUSER_PARTITION(
					tabIndividualCanBeUsedForInitRandom, tabNotInitializedCluster);
			// Compute log-likelihood
			logLikelihood = getLogLikelihood(true); // true : to compute fik
			if (logLikelihood > bestLogLikelihood) {
				bestLogLikelihood = logLikelihood;
				bestParameter->recopy(_parameter);
			}
			/*cout<<endl<<"initRandom"<<endl<<i+1<<" eme essai : "<<endl<<"Parameter : "<<endl;
				_parameter->edit();
				cout<<"LL : "<< logLikelihood<<endl;*/
		}

		// set best parameter
		delete _parameter;
		_parameter = bestParameter;
		_parameter->setModel(this);
		/*cout<<endl<<"initRandom"<<endl<<"meilleur essai : "<<endl<<"Parameter : "<<endl;
			_parameter->edit();
			cout<<"LL : "<< bestLogLikelihood<<endl;*/

		delete [] tabIndividualCanBeUsedForInitRandom;
	}

	delete [] tabNotInitializedCluster;
}

/*------------------------------------------------------
					initSMALL_EM
					------------

  updated in this method :
	- _parameter
	- _tabFik, _tabSumF, _tabCik, _tabTik, _tabNk
      (because an Estep is called to choose the bestParameter)
	Note : _tabFik, _tabSumF, _tabCik, _tabTik, _tabNk wil be 're'computed in the following EStep
	So only _parameter have to be updated in this method

-------------------------------------------------------*/
/*
void Model::initSMALL_EM(ClusteringStrategyInit * clusteringStrategyInit){
//  cout<<"init SMALL_EM, nbTryInInit="<<strategyInit->getNbTry()
//  <<", nbIteration = "<<strategyInit->getNbIteration()<<", epsilon = "<<strategyInit->getEpsilon()<<endl;
  _algoName = EM;
  double logLikelihood, bestLogLikelihood;
  Parameter * bestParameter = _parameter->clone();
  int64_t  i, nbRunOfSmallEMOk = 0;
  bestLogLikelihood = 0.0;
  for (i=0; i<clusteringStrategyInit->getNbTry(); i++){
	nbRunOfSmallEMOk++;
	try{
	  // one run of small EM
	  _parameter->reset();
	  oneRunOfSmallEM(clusteringStrategyInit, logLikelihood);
//cout<<"sortie de oneRunOfSmallEM, LL = "<<logLikelihood<<endl;
	  if ((nbRunOfSmallEMOk == 1) || (logLikelihood > bestLogLikelihood)){
		bestLogLikelihood = logLikelihood;
		bestParameter->recopy(_parameter);
//       cout<<"best LL dans SMALL_EM : " <<bestLogLikelihood<<endl;
	  }
	}
	catch (Exception&errorType){
	  nbRunOfSmallEMOk--;
	}
  }

  if (nbRunOfSmallEMOk == 0){
	THROW(InputException,SMALL_EM_error);
  }

  // set best parameter
  delete _parameter;
  _parameter = bestParameter;
  _parameter->setModel(this);
//  cout<<"fin de init SMALL_EM, nb d'essais effectues="<<i<<endl;
}
 */

//---------------------
// one run if small EM
//--------------------
/*
void Model::oneRunOfSmallEM(ClusteringStrategyInit * clusteringStrategyInit, double & logLikelihood){
  double lastLogLikelihood, eps;
  eps = 1000;
  initRANDOM(1);
  Estep();
  Mstep();
  logLikelihood = getLogLikelihood(true);  // true : to compute fik
  int64_t  nbIteration = 1;
  bool continueAgain = true;
  while (continueAgain){
//    cout<<"while de oneRunOfSmallEM, nbIteration = "<<nbIteration<<endl;
		 //(nbIteration < strategyInit->getNbIteration()) && (eps > strategyInit->getEpsilon())){
	lastLogLikelihood = logLikelihood;
	Estep();
	Mstep();
	nbIteration++;
	// update continueAgain
	switch (clusteringStrategyInit->getStopName()) {
	  case NBITERATION :
		continueAgain = (nbIteration < clusteringStrategyInit->getNbIteration());
		break;
	  case EPSILON :
		logLikelihood = getLogLikelihood(true);  // true : to compute fik
		eps = fabs(logLikelihood - lastLogLikelihood);
		//continueAgain = (eps > strategyInit->getEpsilon());
  		// on ajoute un test pour ne pas faire trop d'iterations quand meme ....
		continueAgain = (eps > clusteringStrategyInit->getEpsilon() && (nbIteration < maxNbIterationInInit));
		break;
	  case NBITERATION_EPSILON :
		logLikelihood = getLogLikelihood(true);  // true : to compute fi
		eps = fabs(logLikelihood - lastLogLikelihood);
		continueAgain = ((eps > clusteringStrategyInit->getEpsilon()) && (nbIteration < clusteringStrategyInit->getNbIteration()));
		break;
		default : THROW(OtherException,internalMixmodError);
	}
  }
  if (clusteringStrategyInit->getStopName() == NBITERATION){ // logLikelihood is an output
	logLikelihood = getLogLikelihood(true);  // true : to compute fi
  }
//cout<<"Fin de oneRunOfSmallEM, nb d'iterations effectuees = "<<nbIteration<<", logLikelihood = "<<logLikelihood<<endl;
}

 */

/*---------------------------------------------------
 initCEM_INIT
 ---------------------------------------------------
  updated in this method :
	- _parameter
	- _tabFik, _tabSumF, _tabCik, _tabTik, _tabNk (because an Estep and a CStep are called to choose the bestParameter)
	Note : _tabFik, _tabSumF, _tabCik, _tabTik, _tabNk wil be 're'computed in the following EStep
	So only _parameter have to be updated in this method

---------------------------------------------------*/
/*
void Model::initCEM_INIT(ClusteringStrategyInit * clusteringStrategyInit){
  //cout<<"init CEM, nbTryInInit="<<strategyInit->getNbTry()<<endl;
  _algoName = CEM;
  int64_t  i;
  double cLogLikelihood, oldLogLikelihood, bestCLogLikelihood;
  Parameter * bestParameter = _parameter->clone();
  int64_t  nbRunOfCEMOk = 0;
  bestCLogLikelihood = 0.0;

  for (i=0; i<clusteringStrategyInit->getNbTry(); i++){
	nbRunOfCEMOk++;
	try{
	  _parameter->reset(); // reset to default values
	  initRANDOM(1);
	  _algoName = CEM;
	  int64_t  nbIter = 0;
	  bool fin = false;
	  while (!fin && nbIter<=maxNbIterationInCEM_INIT){
		Estep();
		Cstep();
		Mstep();
		nbIter++;
		if (nbIter == 1){
		  oldLogLikelihood = getCompletedLogLikelihood();
		}
		else{
		  cLogLikelihood = getCompletedLogLikelihood();
		  if (cLogLikelihood == oldLogLikelihood){
			fin = true;
		  }
		  else{
			oldLogLikelihood = cLogLikelihood;
		  }
		}
	  }
	  //cout<<"dans init CEM, nb d'iterations effectuées : "<<nbIter<<endl;
	// Compute log-likelihood
	  cLogLikelihood = getCompletedLogLikelihood();
	// Comparaison of log-likelihood between step p and p-1
	  if ((nbRunOfCEMOk==1) || (cLogLikelihood > bestCLogLikelihood)){
		bestCLogLikelihood = cLogLikelihood;
		bestParameter->recopy(_parameter);
	  }
	  //cout<<"nbIter : "<<nbIter<<endl;
	}
	catch (Exception&errorType){
	  nbRunOfCEMOk--;
	}
  }

  if (nbRunOfCEMOk==0){
	delete _parameter;
	_parameter = bestParameter;
	_parameter->setModel(this);
	THROW(InputException,CEM_INIT_error);
  }

  //cout<<"fin de init CEM, nb d'essais effectues="<<i<<endl;
  // set Best parameter
  delete _parameter;
  _parameter = bestParameter;
  _parameter->setModel(this);

}
 */
/*---------------------------------------------------------
					Initialization by SEM
					---------------------

  updated in this method :
		- _parameter
		- _tabFik, _tabSumF, _tabCik, _tabTik, _tabNk (because an Estep and a SStep are called to choose the bestParameter)
	Note : _tabFik, _tabSumF, _tabCik, _tabTik, _tabNk wil be 're'computed in the following EStep
	So, only _parameter have to be updated in this method

-------------------------------------------------------*//*
void Model::initSEM_MAX(ClusteringStrategyInit * clusteringStrategyInit){
  //cout<<"init SEM_MAX, nbTryInInit="<<strategyInit->getNbIteration()<<endl;
  _algoName = SEM;
  int64_t  j;
  double logLikelihood, bestLogLikelihood;
  Parameter * bestParameter = _parameter->clone();
  int64_t  nbRunOfSEMMAXOk = 0;
  bestLogLikelihood = 0.0;
//  int64_t  bestIndex=0;

  for (j=0; j<clusteringStrategyInit->getNbIteration(); j++){
    nbRunOfSEMMAXOk++;
    try{
      _parameter->reset();
      initRANDOM(1);
      Estep();
      Sstep();
      Mstep();
      // Compute log-likelihood
      logLikelihood = getLogLikelihood(true);  // true : to compute fik
      if ((nbRunOfSEMMAXOk==1) || (logLikelihood > bestLogLikelihood)){
        bestLogLikelihood = logLikelihood;
        bestParameter->recopy(_parameter);
//        bestIndex = j;
      }
    }
    catch (Exception&errorType){
      nbRunOfSEMMAXOk--;
    }
  }

  if (nbRunOfSEMMAXOk==0){
    THROW(InputException,SEM_MAX_error);
  }

  //cout<<"fin de init SEM_MAX, nb d'iterations effectuees="<<j<<" meilleure solution : "<<bestIndex<<endl;
  // set best parameter
  delete _parameter;
  _parameter = bestParameter;
  _parameter->setModel(this);
}
*/


//-----------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------
// 																	Algorithms
//-----------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------

/*--------------------------------------------------------------------------------------------
	MAP step : Maximum a Posteriori procedure
	---------
  MAP procedure consists in assigning a point to the group maximing this conditional probabiliy

	Note : MAP follows an USER initialisation
	----

	already updated :
	- _parameter (by USER initialisation)

	updated in this method :
	- _tabFik, _tabCik, _tabTik, _tabNk (in Estep called in this method)
	- _tabCik, _tabNk (called by Cstep in this method)
----------------------------------------------------------------------------------------------*/
void Model::MAPstep() {
	Estep(); // to compute _tabFik, _tabTik, _tabCik, _tabNk
	Cstep(); // to compute _tabCik (by MAP procedure) and _tabNk
}

/*--------------------------------------------------------------------------------------------
	E step : Expectation
	------

	Note : Estep follows a Mstep or an initialsation (so _parameter is updated)
	----

	already updated :
	- _parameter

	updated in this method :
	- _tabFik, _tabCik, _tabTik, _tabNk

--------------------------------------------------------------------------------------------*/
void Model::Estep() {

	// 1. compute fik
	//---------------
	computeFik();

	//2. compute tik (conditional probabilities (posteriori probabilities))
	//---------------
	int64_t k;
	int64_t i;

	// Initialize timer info	
	ofstream progressFile;

	for (i = 0; i < _nbSample; i++) {
		if (_tabSumF[i] == 0.0) {
			_parameter->computeTikUnderflow(i, _tabTik);
		}
		else {
			for (k = 0; k < _nbCluster; k++) {
				_tabTik[i][k] = _tabFik[i][k] / _tabSumF[i];
			}
		}

		// 3. compute cik
		//---------------
		if (!_tabZiKnown[i]) {
			for (k = 0; k < _nbCluster; k++) {
				_tabCik[i][k] = _tabTik[i][k];
			}
    }
    //Write progress in file
    if (MASSICCC == 11) {
      progressFile.open ("progress.json");
      progressFile << "{ \"Progress\" :  " << (((double)i + 1)/((double)_nbSample) * 100.0)/2.0 << " }";
      progressFile.close();
    }

  }
	//4. compute nk
	//------------
	computeNk();
}

/*--------------------------------------------------------------------------------------------
	M step
	------

	Note : Mstep follows an Estep, Cstep, Step, or USER_PARTITION initialisation
	----

	already updated :
	- _tabFik, _tabCik, _tabTik, _tabNk

	updated in this method :
	- _parameter
--------------------------------------------------------------------------------------------*/
void Model::Mstep() {
	_parameter->MStep();
}

/*--------------------------------------------------------------------------------------------
	S step
	------
  // S Step : Stochastic Classification

	Note : Sstep follows an Estep
	----

	already updated :
	- _tabFik, _tabCik, _tabTik, _tabNk, _parameter

	updated in this method :
	- _tabCik, _tabNk
--------------------------------------------------------------------------------------------*/
void Model::Sstep() {
	int64_t i;
	int64_t k;
	double ** cumTabT = new double*[_nbSample];

	for (i = 0; i < _nbSample; i++) {
		cumTabT[i] = new double[_nbCluster];
		cumTabT[i][0] = _tabTik[i][0];
	}
	for (k = 1; k < _nbCluster; k++) {
		for (i = 0; i < _nbSample; i++) {
			cumTabT[i][k] = _tabTik[i][k] + cumTabT[i][k - 1];
		}
	}

	double * tabRnd = new double[_nbSample];
	for (i = 0; i < _nbSample; i++) {
		tabRnd[i] = rnd();
	}

	for (i = 0; i < _nbSample; i++) {
		if (!_tabZiKnown[i]) {
			for (k = 0; k < _nbCluster; k++) {
				_tabCik[i][k] = 0;
			}
			k = 0;
			while ((k < _nbCluster) && (tabRnd[i] > cumTabT[i][k])) {
				k++;
			}
			if (tabRnd[i] <= cumTabT[i][k]) {
				_tabCik[i][k] = 1;
			}
			else {
				THROW(OtherException, internalMixmodError);
			}
		}
	}

	for (i = 0; i < _nbSample; i++) {
		delete[] cumTabT[i];
	}
	delete[] cumTabT;
	delete[] tabRnd;

	//update _tabNk
	computeNk();
}

/*--------------------------------------------------------------------------------------------
	C step
	------
  Classification Step

	Note : Cstep follows an Estep
	----

	already updated :
	- _tabFik, _tabCik, _tabTik, _tabNk, _parameter

	updated in this method :
	- _tabCik, _tabNk
--------------------------------------------------------------------------------------------*/
void Model::Cstep() {
	int64_t k, kMax;
	int64_t i;
	double tikMax = 0;

	// Initialize timer info	
	ofstream progressFile;

	for (i = 0; i < _nbSample; i++) {
		if (!_tabZiKnown[i]) {
			kMax = 0;
			tikMax = _tabTik[i][0];
			for (k = 1; k < _nbCluster; k++) {
				if (_tabTik[i][k] > tikMax) {
					tikMax = _tabTik[i][k];
					kMax = k;
				}
			}
			for (k = 0; k < _nbCluster; k++) {
				_tabCik[i][k] = 0;
			}
			_tabCik[i][kMax] = 1;
		}
    //Write progress in file
    if (MASSICCC == 11) {
      progressFile.open ("progress.json");
      progressFile << "{ \"Progress\" :  " << 50 + (((double)i + 1)/((double)_nbSample) * 100.0)/2.0 << " }";
      progressFile.close();
    }
	}

	//update _tabNk
	if (_algoName == UNKNOWN_ALGO_NAME)
		throw;

	if (_algoName != MAP) {
		computeNk();
	}
}

//-----------------------
//  debug information
//-----------------------
void Model::editDebugInformation() {
	if (DEBUG > 0) {
		_parameter->edit();
		if (DEBUG > 1) {
			editFik();
			editTik();
		}
	}
}

//---------
// edit Fik (debug)
//---------
void Model::editFik() {
	int64_t i, k;
	for (i = 0; i < _nbSample; i++) {
		for (k = 0; k < _nbCluster; k++) {
			cout << "\tfik[" << i << "][" << k << "]=" << _tabFik[i][k];
		}
		cout << "\n";
	}
}

//---------
// edit Tik (debug)
//---------
void Model::editTik() {
	int64_t i, k;
	for (i = 0; i < _nbSample; i++) {
		for (k = 0; k < _nbCluster; k++) {
			cout << "\ttik[" << i << "][" << k << "]=" << _tabTik[i][k];
		}
		cout << "\n";
	}
}

//---------
// edit Cik (debug)
//---------
void Model::editCik() {
	int64_t i, k;
	for (i = 0; i < _nbSample; i++) {
		for (k = 0; k < _nbCluster; k++) {
			cout << "\tcik[" << i << "][" << k << "]=" << _tabCik[i][k];
		}
		cout << "\n";
	}
}

//-----------
// edit tabNk (debug)
//-----------
void Model::editNk() {
	int64_t k;
	for (k = 0; k < _nbCluster; k++) {
		cout << "\tnk[" << k << "]=" << _tabNk[k] << "\n";
	}
}

}
