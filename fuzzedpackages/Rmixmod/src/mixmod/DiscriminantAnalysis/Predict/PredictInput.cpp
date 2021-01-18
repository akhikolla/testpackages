/***************************************************************************
                             SRC/mixmod/DiscriminantAnalysis/Predict/PredictInput.cpp  description
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
#include "mixmod/DiscriminantAnalysis/Predict/PredictInput.h"
#include "mixmod/DiscriminantAnalysis/Predict/PredictStrategy.h"
#include "mixmod/Kernel/Parameter/Parameter.h"
#include "mixmod/Kernel/IO/ParameterDescription.h"
#include "mixmod/Kernel/Model/ModelType.h"

namespace XEM {

//--------------------
// Default Constructor
//--------------------
PredictInput::PredictInput() : Input() {
}

//-----------------
//  Copy constructor
//-----------------
PredictInput::PredictInput(const PredictInput & cInput)
: Input(cInput), _classificationRule(cInput.getClassificationRule()) {
}

//---------------------------
// Initialisation Constructor
//---------------------------
PredictInput::PredictInput(DataDescription * predictData, 
		ParameterDescription * classificationRule)
: Input(std::vector<int64_t>(1, classificationRule->getNbCluster()), *predictData)
, _classificationRule(classificationRule->getParameter()) 
{
	// replace default model type by the input model type
	delete _modelType[0];
	_modelType[0] = new ModelType(*classificationRule->getModelType());
}

//-----------
// Destructor
//-----------
PredictInput::~PredictInput() {
}

//------------
//------------
// Modifiers
//------------
//------------

//------ Criterion  ----//

//getCriterion[i]
//-------------------
CriterionName PredictInput::getCriterionName(unsigned int index) const {
	THROW(InputException, notAvailableForPrediction);
}

// removeCriterionName
//--------------------
void PredictInput::removeCriterion(unsigned int index) {
	THROW(InputException, notAvailableForPrediction);
}

//setCriterionName
//----------------
void PredictInput::setCriterion(std::vector<CriterionName> const & criterionName) {
	THROW(InputException, notAvailableForPrediction);
}

//setCriterionName
//----------------
void PredictInput::setCriterion(const CriterionName criterionName, unsigned int index) {
	THROW(InputException, notAvailableForPrediction);
}

// insertCriterionName
//-----------------
void PredictInput::insertCriterion(const CriterionName criterionName, unsigned int index) {
	THROW(InputException, notAvailableForPrediction);
}

// add Criterion
//-----------------
void PredictInput::addCriterion(const CriterionName criterionName) {
	THROW(InputException, notAvailableForPrediction);
}

//------ modelType  ----//

//setModelType
//----------------
void PredictInput::setModelType(const ModelType * modelType, unsigned int index) {
	THROW(InputException, notAvailableForPrediction);
}

// insertModelType
//-----------------
void PredictInput::insertModelType(const ModelType * modelType, unsigned int index) {
	THROW(InputException, notAvailableForPrediction);
}

// add new model type
void PredictInput::addModel(ModelName const modelName) {
	THROW(InputException, notAvailableForPrediction);
}

// removeModelType
//--------------------
void PredictInput::removeModelType(unsigned int index) {
	THROW(InputException, notAvailableForPrediction);
}

// ----------------
// Verif
//-----------------
bool PredictInput::verif() {
	bool res = Input::verif();

	return res;
}

}
