/***************************************************************************
                             SRC/mixmod/Kernel/Parameter/BinaryParameter.cpp  description
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
#include "mixmod/Kernel/Parameter/BinaryParameter.h"
#include "mixmod/Kernel/Parameter/BinaryEkjhParameter.h"
#include "mixmod/Kernel/IO/BinaryData.h"
#include "mixmod/Kernel/IO/BinarySample.h"
#include "mixmod/Kernel/Model/Model.h"
#include "mixmod/Kernel/Model/ModelType.h"
#include "mixmod/Kernel/IO/Partition.h"
#include "mixmod/Utilities/Random.h"

namespace XEM {

//////////////////////////////////
//                              //
//    constructor/destructor    //
//                              //
//////////////////////////////////

//------------
// Constructor
//------------
// Default constructor
BinaryParameter::BinaryParameter() : Parameter() {
	THROW(OtherException, wrongConstructorType);
}

//------------
// Constructor
//------------
BinaryParameter::BinaryParameter(Model * iModel, ModelType * iModelType, int64_t * tabNbModality) 
: Parameter(iModel, iModelType) 
{
	int64_t j, k;

	_pbDimension = iModel->getBinaryData()->getPbDimension();
	// complete scatter
	// tab modality //
	_tabNbModality = copyTab(tabNbModality, _pbDimension);

	// total number of modality //
	_totalNbModality = 0;
	for (j = 0; j < _pbDimension; j++) {
		_totalNbModality += _tabNbModality[j];
	}

	// Centers 
	_tabCenter = new int64_t*[_nbCluster];
	for (k = 0; k < _nbCluster; k++) {
		_tabCenter[k] = new int64_t[_pbDimension];
		for (j = 0; j < _pbDimension; j++) {
			_tabCenter[k][j] = 0;
		}
	}

	if (hasFreeProportion(iModelType->_nameModel)) {
		_freeProportion = true;
	}
	else {
		_freeProportion = false;
	}
}

//------------
// Constructor (if USER initialisation)
//------------
BinaryParameter::BinaryParameter(int64_t iNbCluster, int64_t iPbDimension, 
		ModelType * iModelType, int64_t * tabNbModality) 
: Parameter(iNbCluster, iPbDimension, iModelType) 
{
	int64_t j, k;

	// tab modality //
	_tabNbModality = copyTab(tabNbModality, _pbDimension);

	// total number of modality //
	_totalNbModality = 0;
	for (j = 0; j < _pbDimension; j++) {
		_totalNbModality += _tabNbModality[j];
	}

	// Centers
	_tabCenter = new int64_t*[_nbCluster];
	for (k = 0; k < _nbCluster; k++) {
		_tabCenter[k] = new int64_t[_pbDimension];
		for (j = 0; j < _pbDimension; j++) {
			_tabCenter[k][j] = 0;
		}
	}

	if (hasFreeProportion(iModelType->_nameModel)) {
		_freeProportion = true;
	}
	else {
		_freeProportion = false;
	}
}

//------------
// Constructor
//------------
// copy constructor
BinaryParameter::BinaryParameter(const BinaryParameter * iParameter) : Parameter(iParameter) {

	// tab modality //
	_tabNbModality = copyTab(iParameter->getTabNbModality(), _pbDimension);

	// total number of modality //
	_totalNbModality = iParameter->getTotalNbModality();

	// Centers
	_tabCenter = copyTab(iParameter->getTabCenter(), _nbCluster, _pbDimension);
}

//-----------
// Destructor
//-----------
BinaryParameter::~BinaryParameter() {
	int64_t k;

	if (_tabCenter) {
		for (k = 0; k < _nbCluster; k++) {
			delete[] _tabCenter[k];
			_tabCenter[k] = NULL;
		}
		delete[] _tabCenter;
		_tabCenter = NULL;
	}

	if (_tabNbModality) {
		delete[] _tabNbModality;
		_tabNbModality = NULL;
	}
}

//---------------------
/// Comparison operator
//---------------------
bool BinaryParameter::operator ==(const BinaryParameter & param) const {
	if (!Parameter::operator==(param)) return false;
	if (_totalNbModality != param.getTotalNbModality()) return false;
	for (int64_t k = 0; k < _nbCluster; k++) {
		for (int64_t j = 0; j < _pbDimension; j++) {
			if (_tabCenter[k][j] != param.getTabCenter()[k][j]) return false;
		}
	}
	for (int64_t j = 0; j < _pbDimension; j++) {
		if (_tabNbModality[j] != param.getTabNbModality()[j]) return false;
	}
	return true;
}

//-------------------------
// reset to default values
//------------------------
void BinaryParameter::reset() {
	int64_t k, j;
	for (k = 0; k < _nbCluster; k++) {
		for (j = 0; j < _pbDimension; j++) {
			_tabCenter[k][j] = 0;
		}
	}
	Parameter::reset();
}

//------------------------------------------------------
// getTabCenterIfOneCluster
//------------------------------------------------------
// Compute Center for One cluster 
// _tabCenter will not be changed
void BinaryParameter::getTabCenterIfOneCluster(int64_t * tabCenter, 
		double * tabNbSampleInMajorModality, double ** tabNbSamplePerModality) const 
{
	int64_t i;
	int64_t j, h;
	double max, nbSamplePerModality;
	int64_t nbSample = _model->getNbSample();

	Data * data = _model->getData();
	Sample ** dataMatrix = data->_matrix;
	BinarySample * curSample;

	for (j = 0; j < _pbDimension; j++) {
		max = 0.0;
		for (h = 1; h <= _tabNbModality[j]; h++) {
			nbSamplePerModality = 0.0;
			for (i = 0; i < nbSample; i++) {
				curSample = dataMatrix[i]->getBinarySample();
				if (curSample->getDataValue(j) == h)
					nbSamplePerModality += data->_weight[i];
			}
			// test if it is majority modality //
			if (nbSamplePerModality > max) {
				max = nbSamplePerModality;
				tabCenter[j] = h;
			}
			if (tabNbSamplePerModality != NULL) {
				tabNbSamplePerModality[j][h - 1] = nbSamplePerModality;
			}
		}
		tabNbSampleInMajorModality[j] = max;
	}
}

//////////////////////////////////
//                              //
//      compute method          //
//                              //
//////////////////////////////////
void BinaryParameter::getAllPdf(double** tabFik, double* tabProportion) const {
	int64_t nbSample = _model->getNbSample();
	int64_t i;
	int64_t k;

	for (i = 0; i < nbSample; i++) {
		for (k = 0; k < _nbCluster; k++) {
			tabFik[i][k] = getPdf(i, k) * tabProportion[k];
		}
	}
}

//-----------------------------------------
// computeTik when underflow
// -model->_tabSumF[i] pour ith sample = 0
// i : 0 ->_nbSample-1
//-----------------------------------------
void BinaryParameter::computeTikUnderflow(int64_t i, double ** tabTik) {
	THROW(NumericException, sumFiNullInMultinomialCase);
	/*
	long double * lnFk      = new long double[_nbCluster];
	long double * lnFkPrim  = new long double[_nbCluster];
	long double * fkPrim    = new long double[_nbCluster];
	double *  tabTik_i = tabTik[i];  
	int64_t k,k0;
	long double lnFkMax, fkTPrim;
    
	for (k=0; k<_nbCluster; k++){
	  lnFk[k]       =  getLogPdf(i,k);
	} 

	lnFkMax = lnFk[0];
	for (k=1; k<_nbCluster; k++){
	  if (lnFk[k] > lnFkMax){
		lnFkMax = lnFk[k];
	  }
	}

	fkTPrim = 0.0;
	for (k=0; k<_nbCluster; k++){
	  lnFkPrim[k] = lnFk[k] - lnFkMax;
	  fkPrim[k]   = exp(lnFkPrim[k]);
	  fkTPrim    += fkPrim[k];
	}

	// compute tabTik
	if (fkTPrim == 0){
	  THROW(XEMInputException,sumFiNullAndfkTPrimNull);
	  // reset tabTik
  //     initToZero(tabTik_i, _nbCluster);
  //     k0           = XEMGaussianParameter::computeClassAssigment(i);
  //     tabTik_i[k0] = 1.0;
	}
	else{
	  for (k=0; k<_nbCluster; k++){
		tabTik_i[k] = fkPrim[k] / fkTPrim;
	  }
	}
  
	delete [] lnFk;
	delete [] lnFkPrim;
	delete [] fkPrim;*/
}

//-----------------------------------------------------
// Compute table of centers of samples for each cluster
// outputs :
// -  nbInitializedCluster
// - tabInitializedCluster (array of size _nbCluster)
//-----------------------------------------------------
void BinaryParameter::computeTabCenterInitUSER_PARTITION(int64_t & nbInitializedCluster, 
		bool * tabInitializedCluster, Partition * initPartition) 
{
	int64_t j, k, h;
	int64_t i;
	double argMax, bestArgMax;

	double ** tabTik = _model->getTabTik();
	int64_t ** initPartitionValue = initPartition->_tabValue;
	BinaryData * data = _model->getBinaryData();
	Sample ** dataMatrix = data->getDataMatrix();
	BinarySample * curSample;
	int64_t * tabNbModality = data->getTabNbModality();
	int64_t nbSample = _model->getNbSample();

	for (k = 0; k < _nbCluster; k++) {
		//cout<<"k = "<<k<<endl;
		for (j = 0; j < _pbDimension; j++) {
			//cout<<"j = "<<j<<endl;
			bestArgMax = 0.0;
			_tabCenter[k][j] = 0;
			for (h = 1; h <= tabNbModality[j]; h++) {
				//cout<<"h = "<<h<<endl;
				argMax = 0.0;
				for (i = 0; i < nbSample; i++) {
					curSample = dataMatrix[i]->getBinarySample();
					if (curSample->getDataValue(j) == h) {
						if (initPartitionValue[i][k] == 1) {
							argMax += data->_weight[i];
						}
						else {
							argMax += tabTik[i][k] * data->_weight[i];
						}
					}
					//cout<<"argMax = "<<argMax<<endl;
				} // end for i

				if (argMax > bestArgMax) {
					bestArgMax = argMax;
					_tabCenter[k][j] = h;
				}
			} // end for h
			//cout<<"bestArgMax = "<<bestArgMax<<endl;
		} // end for j
	} // end for k

	nbInitializedCluster = _nbCluster;
	for (k = 0; k < _nbCluster; k++) {
		tabInitializedCluster[k] = true;
	}
}

void BinaryParameter::computeTabCenter() {
	int64_t j, k, h;
	int64_t i;
	double argMax, bestArgMax;

	double * tabNk = _model->getTabNk();
	double ** tabCik = _model->getTabCik();

	BinaryData * data = _model->getBinaryData();
	Sample ** dataMatrix = data->getDataMatrix();
	BinarySample * curSample;
	int64_t * tabNbModality = data->getTabNbModality();
	int64_t nbSample = _model->getNbSample();

	for (k = 0; k < _nbCluster; k++) {
		for (j = 0; j < _pbDimension; j++) {

			bestArgMax = 0.0;
			for (h = 1; h <= tabNbModality[j]; h++) {
				argMax = 0.0;
				for (i = 0; i < nbSample; i++) {
					curSample = dataMatrix[i]->getBinarySample();
					if (curSample->getDataValue(j) == h) {
						argMax += tabCik[i][k] * data->_weight[i];
					}
				} // end for i

				argMax /= tabNk[k];

				if (argMax > bestArgMax) {
					bestArgMax = argMax;
					_tabCenter[k][j] = h;
				}

			} // end for h

		} // end for j
	} // end for k
}

//-----------------------------
// updateForOneRunOfInitRANDOM
//-----------------------------
void BinaryParameter::updateForInitRANDOMorUSER_PARTITION(
		Sample ** tabSampleForInit, bool * tabClusterToInitialize) 
{
	// set _tabCenter value
	int64_t * sampleValue = NULL;
	Sample * curSample = NULL;
	for (int64_t k = 0; k < _nbCluster; k++) {
		if (tabClusterToInitialize[k]) {
			curSample = tabSampleForInit[k];
			sampleValue = (curSample->getBinarySample())->getTabValue();
			recopyTab(sampleValue, _tabCenter[k], _pbDimension);
		}
	}

	// set _tabScatter value
	computeRandomScatter(); // virtual
}


//////////////////////////////////
//                              //
//    initialization method     //
//                              //
//////////////////////////////////

//----------
// init USER
//----------
// note : iParam est de type XEMBinaryEkjhParameter (et proportion libre) 
void BinaryParameter::initUSER(Parameter * iParam) {
	/*
	 updated : _tabProportion, _tabCenter, _scatter
	 */
	int64_t j, k;
	BinaryParameter * param = iParam->getBinaryParameter();


	double * iTabProportion = param->getTabProportion();
	int64_t ** iTabCenter = param->getTabCenter();
	int64_t * iTabNbModality = param->getTabNbModality();

	_totalNbModality = param->getTotalNbModality();
	for (j = 0; j < _pbDimension; j++)
		_tabNbModality[j] = iTabNbModality[j];

	for (k = 0; k < _nbCluster; k++) {
		// proportion (no respecting model type)
		if (!hasFreeProportion(_modelType->_nameModel)) {
			_tabProportion[k] = 1.0 / _nbCluster;
		}
		else {
			_tabProportion[k] = iTabProportion[k];
		}

		// center
		for (j = 0; j < _pbDimension; j++) {
			_tabCenter[k][j] = iTabCenter[k][j];
		}
	}

	// scatters
	if (typeid (*param) == typeid (*this)) {
		recopyScatter(param); // virtual
	}
	else {
		createScatter(param->scatterToArray());
	}
}

//-------------------------------------------
// initialize attributes before an InitRandom 
//-------------------------------------------
void BinaryParameter::initForInitRANDOM() {
	// none !
}

//---------------------------
// initForInitUSER_PARTITION
// outputs :
// -  nbInitializedCluster
// - tabInitializedCluster (array of size _nbCluster)
//-------------------------
void BinaryParameter::initForInitUSER_PARTITION(int64_t & nbInitializedCluster, 
		bool * tabInitializedCluster, Partition * initPartition) 
{
	// Set Centers 
	computeTabCenterInitUSER_PARTITION(nbInitializedCluster, tabInitializedCluster, initPartition);

	if (nbInitializedCluster == _nbCluster) {
		computeRandomScatter();
	}
	// sinon c'est fait plus tard 
}


//////////////////////////////////
//                              //
//         algorithms           //
//                              //
//////////////////////////////////

/*-------------------------------------------------------------------------------------------
	M step
	------
	
	already updated :
	- _model (_tabFik, _tabTik, _tabCik, _tabNk) if Estep is done before
	- _model->_tabCik, _model->_tabNk if USER_PARTITION is done before (Disciminant analysis)
	In all cases, only  _tabCik and _tabNk are needed 
	
	updated in this method :
	- _tabProportion
	- _tabCenter
	- _tabScatter
--------------------------------------------------------------------------------------------*/
void BinaryParameter::MStep() {
	/* Proportion estimator (if proportions free) */
	computeTabProportion();

	/* Centers centers */
	computeTabCenter();

	/* Compute Scatter */
	computeScatter(); // virtual
}


//////////////////////////////////
//                              //
//        input/output          //
//                              //
//////////////////////////////////

//--------------------------
// Edit  : debug information
//--------------------------
void BinaryParameter::edit() {
	int64_t k, i;
	for (k = 0; k < _nbCluster; k++) {
		cout << "\tcomponent : " << k << endl;
		cout << "\tproportion : " << _tabProportion[k] << endl;
		cout << "\tcenter : ";
		for (i = 0; i < _pbDimension; i++) {
			cout << "\t" << _tabCenter[k][i];
		}
		cout << endl;
		cout << "\tscatter : ";
		editScatter(k); // virtual
		cout << endl;
	}
}

//-----
// Edit
//-----
void BinaryParameter::edit(std::ofstream & oFile, bool text) {
	int64_t j, k;

	// Write in output text files
	if (text) {
		for (k = 0; k < _nbCluster; k++) {

			oFile << "\t\t\tComponent " << k + 1 << endl;
			oFile << "\t\t\t---------" << endl;

			oFile << "\t\t\tMixing proportion : " << _tabProportion[k] << endl;

			oFile << "\t\t\tCenter : ";
			for (j = 0; j < _pbDimension; j++)
				oFile << _tabCenter[k][j] << "  ";
			oFile << endl;

			editScatter(oFile, k, text); // virtual
			oFile << endl << endl;

		}// end for k
		oFile << endl;
	}

	// Write in output numeric files
  else {
    for (k = 0; k < _nbCluster; k++) {
			putDoubleInStream(oFile, _tabProportion[k]);

      for (j = 0; j < _pbDimension; j++)
        oFile << _tabCenter[k][j] << "  ";
      oFile << endl;
      editScatter(oFile, k, text); // virtual
      oFile << endl << endl;

    }// end for k
    oFile << endl;
  }
}

//------------------------------
// Read Parameters in input file
//------------------------------
void BinaryParameter::input(std::ifstream & fi) {
	int64_t j, k;

  for (k = 0; k < _nbCluster; k++) {

    // Proportions //
		_tabProportion[k] = getDoubleFromStream(fi);

    // Centers //
    for (j = 0; j < _pbDimension; j++)
			fi >> _tabCenter[k][j];

		// Scatters //
		inputScatter(fi, k); // virtual
	}
}

//------------------------------
// Read Parameters in input containers
//------------------------------
void BinaryParameter::input(
		double * proportions, 
		double ** centers, 
		double *** scatters) 
{
	int64_t j, k;

	for (k = 0; k < _nbCluster; k++) {

		// Proportions //
		_tabProportion[k] = proportions[k];

		// Centers //
		for (j = 0; j < _pbDimension; j++) {
			_tabCenter[k][j] = centers[k][j];
		}

		// Scatters //
		inputScatter(scatters); // virtual
	}
}

void BinaryParameter::updateForCV(Model * originalModel, CVBlock & CVBlock) {
	//updates tabProportion, tabCenter and tabScatter
	// Mstep could be could to do that
	// even if this is a few slower function versus a real update
	MStep();
}

//------------------//
// Others functions //
//------------------//
// compute Pdf in case nbCluster=1 (Scatter is a scalar)
double computePdfOneCluster(Sample * x, int64_t * Center, 
		double Scatter, int64_t * tabNbModality) 
{
	double bernPdf = 1.0;
	int64_t j;
	int64_t pbDimension = x->getPbDimension();
	for (j = 0; j < pbDimension; j++) {
		if ((x->getBinarySample())->getDataValue(j) == Center[j])
			bernPdf *= 1.0 - Scatter;
		else
			bernPdf *= Scatter / (tabNbModality[j] - 1);
	}
	return bernPdf;
}

// compute Pdf in case nbCluster=1 (Scatter is a array of double, depends on variables)
double computePdfOneCluster(Sample * x, int64_t * Center, 
		double * Scatter, int64_t * tabNbModality) 
{
	double bernPdf = 1.0;
	int64_t j;
	int64_t pbDimension = x->getPbDimension();
	for (j = 0; j < pbDimension; j++) {
		if ((x->getBinarySample())->getDataValue(j) == Center[j])
			bernPdf *= 1.0 - Scatter[j];
		else
			bernPdf *= Scatter[j] / (tabNbModality[j] - 1);
	}
	return bernPdf;
}

// compute Pdf in case nbCluster=1 
// (Scatter is a array of double*double, depends on variables and modalities)
double computePdfOneCluster(Sample * x, int64_t * Center, 
		double ** Scatter, int64_t * tabNbModality) 
{
	double bernPdf = 1.0;
	int64_t j;
	int64_t pbDimension = x->getPbDimension();
	int64_t value;
	for (j = 0; j < pbDimension; j++) {
		value = (x->getBinarySample())->getDataValue(j);
		if (value == Center[j])
			bernPdf *= 1.0 - Scatter[j][value - 1];
		else
			bernPdf *= Scatter[j][value - 1];
	}
	return bernPdf;
}

void BinaryParameter::recopy(Parameter * otherParameter) {

	BinaryParameter * iParameter = otherParameter->getBinaryParameter();

	// tab modality //
	recopyTab(iParameter->getTabNbModality(), _tabNbModality, _pbDimension);

	// total number of modality //
	_totalNbModality = iParameter->getTotalNbModality();

	// Centers 
	recopyTab(iParameter->getTabCenter(), _tabCenter, _nbCluster, _pbDimension);

	// Scatters
	recopyScatter(otherParameter->getBinaryParameter()); // virtual
}

}
