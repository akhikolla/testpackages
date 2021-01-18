/***************************************************************************
                             SRC/mixmod/Kernel/IO/ModelOutput.cpp  description
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

#include "mixmod/Kernel/IO/ModelOutput.h"
#include "mixmod/Kernel/Model/Model.h"
#include "mixmod/Kernel/Model/ModelType.h"
#include "mixmod/Kernel/IO/ProbaDescription.h"
#include "mixmod/Kernel/IO/LabelDescription.h"

namespace XEM {

//--------------------
// Default Constructor
//--------------------
ModelOutput::ModelOutput() {
}

//-----------------
//  Copy constructor
//-----------------
ModelOutput::ModelOutput(const ModelOutput & modelOutput) {
	THROW(OtherException, internalMixmodError);
}

//------------------------------
//  Initialization constructor 1
//------------------------------
ModelOutput::ModelOutput(Model * estimation) {
	if (!estimation) {
		THROW(OtherException, nullPointerError);
	}
	_modelType = *(estimation->getModelType());
	_nbCluster = estimation->getNbCluster();
	//_strategyRunError = &(estimation->getErrorType());
	//changed to clone instead of reference
	_strategyRunError = ((estimation->getErrorType())).clone();
	if (*(dynamic_cast<Exception*> (_strategyRunError)) == NOERROR) {
		_probaDescription = new ProbaDescription(estimation);
		_labelDescription = new LabelDescription(estimation);
		_parameterDescription = new ParameterDescription(estimation);
	}
	else {
		_probaDescription = NULL;
		_labelDescription = NULL;
		_parameterDescription = NULL;
	}

	_likelihood = estimation->getLogLikelihood(false);
}

//------------------------------
//  Initialization constructor 2
//------------------------------
ModelOutput::ModelOutput(
		ModelType & modelType, 
		int64_t nbCluster, 
		std::vector<CriterionOutput*> & criterionOutput, 
		double likelihood, 
		ParameterDescription & parameterDescription, 
		LabelDescription & labelDescription, 
		ProbaDescription & probaDescription) 
{
  _modelType = modelType;
  _nbCluster = nbCluster;
  _strategyRunError = NOERROR.clone(); // TODO ??
  if (*(dynamic_cast<Exception*> (_strategyRunError)) == NOERROR) {
    _probaDescription = new ProbaDescription(probaDescription);
    _labelDescription = new LabelDescription(labelDescription);
    _parameterDescription = new ParameterDescription(parameterDescription);
    for (int64_t i = 0; i < criterionOutput.size(); i++) {
      _criterionOutput[i].setValue(criterionOutput[i]->getValue());
      _criterionOutput[i].setCriterionName(criterionOutput[i]->getCriterionName());
      _criterionOutput[i].setError(criterionOutput[i]->getError());
    }
  }
  else {
    _probaDescription = NULL;
    _labelDescription = NULL;
    _parameterDescription = NULL;
  }

  _likelihood = likelihood;
}

//------------------------------
//  Initialization constructor 3
//------------------------------
ModelOutput::ModelOutput(ModelType & modelType, int64_t nbCluster, Exception& error) {
	_modelType = modelType;
	_nbCluster = nbCluster;
	_strategyRunError = error.clone();
	_probaDescription = NULL;
	_labelDescription = NULL;
	_parameterDescription = NULL;
	_likelihood = 0;
}

//-----------
// Destructor
//-----------
ModelOutput::~ModelOutput() {
	if (_labelDescription != NULL) delete _labelDescription;
	if (_parameterDescription != NULL) delete _parameterDescription;
	if (_probaDescription != NULL) delete _probaDescription;
	if (_strategyRunError != NULL) delete _strategyRunError;
}

//----------------------
/// Comparison operator
//----------------------
bool ModelOutput::operator ==(const ModelOutput & modelOutput) const {

	if (_nbCluster != modelOutput.getNbCluster()) return false;
	if (!(_modelType == modelOutput.getModelType())) return false;
	for (int iCriterion = 0; iCriterion < maxNbCriterion; iCriterion++) {
		if (!(_criterionOutput[iCriterion] == modelOutput.getCriterionOutput(iCriterion))) 
			return false;
	}
	if (!(_parameterDescription == modelOutput.getParameterDescription())) return false;
	if (!(_labelDescription == modelOutput.getLabelDescription())) return false;
	if (!(_probaDescription == modelOutput.getProbaDescription())) return false;
	return true;
}

// set criterion output
void ModelOutput::setCriterionOutput(CriterionOutput const & criterionOutput) {
	// get criterion name
	CriterionName criterionName = criterionOutput.getCriterionName();
	_criterionOutput[criterionName].setCriterionName(criterionName);
	_criterionOutput[criterionName].setValue(criterionOutput.getValue());
	_criterionOutput[criterionName].setError(criterionOutput.getError());
}

}
