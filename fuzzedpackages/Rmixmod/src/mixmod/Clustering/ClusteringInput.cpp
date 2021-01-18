/***************************************************************************
                             SRC/mixmod/Clustering/ClusteringInput.cpp  description
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

#include "mixmod/Clustering/ClusteringInput.h"
#include "mixmod/Clustering/ClusteringStrategy.h"
#include "mixmod/Kernel/Model/ModelType.h"

namespace XEM {

//--------------------
// Default Constructor
//--------------------
ClusteringInput::ClusteringInput() {
	_strategy = new ClusteringStrategy();
}

//-----------------
//  Copy constructor
//-----------------
ClusteringInput::ClusteringInput( const ClusteringInput & cInput )
: Input(cInput)
{
	_strategy = new ClusteringStrategy(*cInput.getStrategy());
}

//---------------------------
// Initialisation Constructor
//---------------------------
ClusteringInput::ClusteringInput( const std::vector<int64_t> & iNbCluster,
		const DataDescription & iDataDescription)
: Input(iNbCluster, iDataDescription)
{
	_strategy = new ClusteringStrategy();
}

//-----------
// Destructor
//-----------
ClusteringInput::~ClusteringInput() {
	if ( _strategy ) {
		delete _strategy;
		_strategy = NULL;
	}
}

//-----------
// setStrategy
//-----------
void ClusteringInput::setStrategy(ClusteringStrategy * strat) {
	_strategy = new ClusteringStrategy(*strat);
}

//setCriterion
//----------------
void ClusteringInput::setCriterion(const CriterionName criterionName, unsigned int index) {
	if (index < _criterionName.size()) {
		switch (criterionName) {
		case BIC: _criterionName[index] = BIC;
			break;
		case CV:  THROW(InputException, DAInput);
		case ICL: _criterionName[index] = ICL;
			break;
		case NEC: _criterionName[index] = NEC;
			break;
		case UNKNOWN_CRITERION_NAME: THROW(OtherException, internalMixmodError);
			break;
		default: THROW(OtherException, internalMixmodError);
		}
	}
	else {

		THROW(InputException, wrongCriterionPositionInSet);
	}
	_finalized = false;
}

//setCriterion
//----------------
void ClusteringInput::setCriterion(std::vector<CriterionName> const & criterionName) {

	// copy vector contents
	_criterionName = criterionName;

	// check vector contents
	for ( unsigned int iCriterion = 0; iCriterion < _criterionName.size(); iCriterion++ ) {
		switch (_criterionName[iCriterion]) {
		case BIC: break;
		case CV:  THROW(InputException, DAInput);
		case ICL: break;
		case NEC: break;
		case UNKNOWN_CRITERION_NAME: THROW(OtherException, internalMixmodError);
			break;
		default: THROW(OtherException, internalMixmodError);
		}
	}
	_finalized = false;
}

// insertCriterion
//-----------------

/* Note Bug 15361
When adding a new criterion via mixmodGUI, we do :
1. insertCriterion(UNKNOWN_CRITERION_NAME, idnex)
2. setCriterion(<theGoodCriterionName>, index)
so we have to remove criterionName tests in this method
A new task 5290 has been created to try to do only a insert
 */
void ClusteringInput::insertCriterion(const CriterionName criterionName, unsigned int index) {
	if (index >= 0 && index <= _criterionName.size()) {
		switch (criterionName) {
		case BIC:
			_criterionName.insert(_criterionName.begin() + index , BIC);
			break;
		case CV:
			THROW(InputException, DAInput);
		case ICL:
			_criterionName.insert(_criterionName.begin() + index , ICL);
			break;
		case NEC:
			_criterionName.insert(_criterionName.begin() + index , NEC);
			break;
			/*Correction bug 15361
			case UNKNOWN_CRITERION_NAME : THROW(XEMOtherException,internalMixmodError);break;*/
		case UNKNOWN_CRITERION_NAME:
			_criterionName.insert(_criterionName.begin() + index , UNKNOWN_CRITERION_NAME);
			break;
		default:
			THROW(OtherException, internalMixmodError);
		}
	}
	else {
		THROW(InputException, wrongCriterionPositionInInsert);
	}
	_finalized = false;
}

// add Criterion
//-----------------
void ClusteringInput::addCriterion(const CriterionName criterionName) {

	bool found = false;
	for (unsigned int iCriterion = 0; iCriterion < _criterionName.size(); iCriterion++ ) {
		if ( _criterionName[iCriterion] == criterionName ) found = true;
	}
	if (!found) {
		switch (criterionName) {
		case BIC: _criterionName.push_back(BIC);
			break;
		case CV:  THROW(InputException, DAInput);
		case ICL: _criterionName.push_back(ICL);
			break;
		case NEC: _criterionName.push_back(NEC);
			break;
		case UNKNOWN_CRITERION_NAME: THROW(OtherException, internalMixmodError);
			break;
		default: THROW(OtherException, internalMixmodError);
		}
	}
	_finalized = false;
}

//------------
// print input
//------------
void ClusteringInput::edit(std::ostream & out ) const {
	Input::edit(out);
	_strategy->edit(out);
}


// ----------------
// Verif
//-----------------
bool ClusteringInput::verif() {
	bool res = Input::verif();
	if (res) {
		res = _strategy->verify();
	}
	return res;
}

/// setModelType
void ClusteringInput::setModelType(const ModelType * modelType, unsigned int index){
	if (isHD(modelType->getModelName())){
		THROW(InputException, HDModelsAreNotAvailableInClusteringContext);
	}
	else
		Input::setModelType(modelType, index);
}


/// insertModelType
void ClusteringInput::insertModelType(const ModelType * modelType, unsigned int index){
	if (isHD(modelType->getModelName())){
		THROW(InputException, HDModelsAreNotAvailableInClusteringContext);
	}
	else
		Input::insertModelType(modelType, index);
}


/// add new model type (at the end)
void ClusteringInput::addModelType(const ModelType * modelType){
	if (isHD(modelType->getModelName())){
		THROW(InputException, HDModelsAreNotAvailableInClusteringContext);
	}
	else
		Input::addModelType(modelType);
}

/// add new model (modelName -> modelType)
void ClusteringInput::addModel(ModelName const modelName){
	if (isHD(modelName)){
		THROW(InputException, HDModelsAreNotAvailableInClusteringContext);
	}
	else
		Input::addModel(modelName);
}

/// setModel (modelName -> modelType)
void ClusteringInput::setModel(std::vector<ModelName> const & modelName){
	for (unsigned int iModel = 0; iModel < modelName.size(); iModel++) {
		if (isHD(modelName[iModel])){
			THROW(InputException, HDModelsAreNotAvailableInClusteringContext);
		}
	}
	Input::setModel(modelName);
}



}
