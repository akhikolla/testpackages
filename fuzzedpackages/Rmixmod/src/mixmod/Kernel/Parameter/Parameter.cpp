/***************************************************************************
                             SRC/mixmod/Kernel/Parameter/Parameter.cpp  description
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
#include "mixmod/Kernel/Parameter/Parameter.h"
#include "mixmod/Utilities/Random.h"
#include "mixmod/Kernel/Model/Model.h"
#include "mixmod/Kernel/IO/Data.h"

namespace XEM {

//------------
// Constructor
//------------
// Default constructor
Parameter::Parameter() {
	THROW(OtherException, wrongConstructorType);
}

//------------
// Constructor
// called if USER initialisation
//------------
Parameter::Parameter(int64_t iNbCluster, int64_t iPbDimension, ModelType * iModelType) {
	_modelType = iModelType;
	_nbCluster = iNbCluster;
	_pbDimension = iPbDimension;
	_tabProportion = new double[_nbCluster];
	for (int64_t k = 0; k < _nbCluster; k++) {
		_tabProportion[k] = 1.0 / _nbCluster;
	}
	_model = NULL;
	_filename = "";
	_format = FormatNumeric::defaultFormatNumericFile;
	_freeProportion = true; //[bauder] TODO: or false? true seems more logical but...
}

//------------
// Constructor
//------------
// called by XEMModel::XEMModel
// ! iModel isn't updated (iModel._param is this)
Parameter::Parameter(Model * iModel, ModelType * iModelType) {
	_modelType = iModelType;
	_model = iModel;
	_nbCluster = _model->getNbCluster();
	//Now pbDimention gets updated in derived classes.
	//_pbDimension   = _model->getData()->_pbDimension ;
	_tabProportion = new double[_nbCluster];

	// _tabProportion isn't computed because model->_tabNk isn't updated (=0)
	for (int64_t k = 0; k < _nbCluster; k++) {
		_tabProportion[k] = 1.0 / _nbCluster;
	}
	_filename = "";
	_format = FormatNumeric::defaultFormatNumericFile;
	_freeProportion = true; //[bauder] TODO: or false? true seems more logical but...
}

//------------
// Constructor
//------------
Parameter::Parameter(const Parameter * iParameter) {
	_nbCluster = iParameter->getNbCluster();
	_pbDimension = iParameter->getPbDimension();
	_tabProportion = copyTab(iParameter->getTabProportion(), _nbCluster);
	if (iParameter->getModel() != NULL)
		_model = iParameter->getModel();
	else
		_model = NULL;
	if (iParameter->getModelType() != NULL)
		_modelType = iParameter->getModelType();
	else
		_modelType = NULL;
	_freeProportion = iParameter->getFreeProportion();
	_filename = iParameter->getFilename();
	_format = iParameter->getFormat();
	_freeProportion = iParameter->getFreeProportion();
}

//-----------
// Destructor
//-----------
Parameter::~Parameter() {

	if (_tabProportion) {
		delete[] _tabProportion;
		_tabProportion = NULL;
	}
}

//---------------------
/// Comparison operator
//---------------------
bool Parameter::operator ==(const Parameter & param) const {
	if (_pbDimension != param.getPbDimension()) return false;
	if (_nbCluster != param.getNbCluster()) return false;
	if (_freeProportion != param.getFreeProportion()) return false;
	for (int64_t k = 0; k < _nbCluster; k++) {
		if (_tabProportion[k] != param.getTabProportion()[k]) return false;
	}
	return true;
}

//------------------------
// reset to default values
//------------------------
void Parameter::reset() {
	for (int64_t k = 0; k < _nbCluster; k++) {
		_tabProportion[k] = 1.0 / _nbCluster;
	}
}

//-------------------
//generateRandomIndex
//-------------------
int64_t Parameter::generateRandomIndex(bool * tabIndividualCanBeUsedForInitRandom, 
		double * weight, double totalWeight) 
{
	double rndWeight, sumWeight;
	int64_t idxSample;

	/* Generate a random integer between 0 and _nbSample-1 */
	bool IdxSampleCanBeUsed = false; // idxSample can be used
	while (!IdxSampleCanBeUsed) {
		// get index of sample with weight //
		rndWeight = (int64_t) (totalWeight * rnd() + 1);
		sumWeight = 0.0;
		idxSample = -1;
		while (sumWeight < rndWeight) {
			idxSample++;
			sumWeight += weight[idxSample];
		}
		//cout<<"index tire au hasard :"<<idxSample<<endl;
		IdxSampleCanBeUsed = tabIndividualCanBeUsedForInitRandom[idxSample];
	}
	// on indique que cet individu ne pourra pas etre tire au  hasard pour une autre classe
	tabIndividualCanBeUsedForInitRandom[idxSample] = false;
	//cout<<"choisi"<<endl;
	return idxSample;
}


//-----------------
//-----------------
// compute methods
//-----------------
//-----------------

//-----------------------
// compute tab proportion
//-----------------------
void Parameter::computeTabProportion() {
	int64_t k;
	double * tabNk = _model->getTabNk();
	double weightTotal = (_model->getData())->_weightTotal;

	if (_freeProportion) {
		for (k = 0; k < _nbCluster; k++)
			_tabProportion[k] = tabNk[k] / weightTotal;
	}
	else {
		for (k = 0; k < _nbCluster; k++)
			_tabProportion[k] = 1.0 / _nbCluster;
	}
}

//-----------------------------------------
// computeTik when underflow
// -model->_tabSumF[i] pour ith sample = 0
// i : 0 ->_nbSample-1
//-----------------------------------------
void Parameter::computeTikUnderflow(int64_t i, double ** tabTik) {
	THROW(OtherException, nonImplementedMethod);
}

}
