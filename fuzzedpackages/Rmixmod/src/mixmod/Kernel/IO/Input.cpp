/***************************************************************************
                             SRC/mixmod/Kernel/IO/Input.cpp  description
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

#include "mixmod/Kernel/IO/Input.h"
#include "mixmod/Kernel/IO/Data.h"
#include "mixmod/Kernel/IO/Partition.h"
#include "mixmod/Kernel/IO/LabelDescription.h"
#include "mixmod/Kernel/Model/ModelType.h"

namespace XEM {

//--------------------
// Default Constructor
//--------------------
Input::Input() {
	_nbSample = 0;
	_pbDimension = 0;
	_nbCluster.clear();
	//_knownPartition   = NULL;
	_knownLabelDescription = NULL;
	_finalized = false;
	//_dataDescription will be created by default DataDescription constructor
}

//-----------------
// Copy Constructor
//-----------------
Input::Input(const Input & input) {
	_finalized = input._finalized;

	_nbSample = input._nbSample;
	_pbDimension = input._pbDimension;
	_nbCluster = input._nbCluster;

	_dataDescription = (input._dataDescription); // op =


	_knownPartition = NULL;
	if (input._knownPartition != NULL) {
		_knownPartition = new Partition(input._knownPartition);
	}

//if (!_knownLabelDescription) ??
//  delete _knownLabelDescription;

	_knownLabelDescription = NULL;
	if (input._knownLabelDescription != NULL) {
		_knownLabelDescription = new LabelDescription(*input._knownLabelDescription);
	}

	_criterionName = input._criterionName;
	_modelType = input._modelType;

}

//---------------------------
// Initialisation Constructor
//---------------------------
Input::Input(const std::vector<int64_t> & iNbCluster, const DataDescription & iDataDescription) {
	cloneInitialisation(iNbCluster, iDataDescription);
}

//---------------------
// Clone Initialisation 
//---------------------
void Input::cloneInitialisation(const std::vector<int64_t> & iNbCluster, 
		const DataDescription & iDataDescription) 
{
	_finalized = false;

	_nbSample = iDataDescription.getNbSample();
	_pbDimension = iDataDescription.getPbDimension();
	_nbCluster = iNbCluster; // copy Constructor

	//_data
	_dataDescription = iDataDescription;

	_knownPartition = NULL;
	_knownLabelDescription = NULL;
	_criterionName.push_back(defaultCriterionName);

	/*if ( _dataDescription.getDataType()==false ){
	  _modelType.push_back(new ModelType());
	}
	else{
	  _modelType.push_back(new ModelType(defaultBinaryModelName));
	}*/
	//cout<<_dataDescription.getDataType()<<endl;
	switch (_dataDescription.getDataType()) {
	case QualitativeData:
		_modelType.push_back(new ModelType(defaultBinaryModelName));
		break;
	case QuantitativeData:
		_modelType.push_back(new ModelType());
		break;
	case HeterogeneousData:
		_modelType.push_back(new ModelType(defaultHeterogeneousModelName));
		break;
	default:
		break;
	}
}

//---------------------------------
// Constructor used in DCV context
//--------------------------------
Input::Input(Input * originalInput, CVBlock & learningBlock) {
	THROW(OtherException, internalMixmodError);
}

//-----------
// Destructor
//-----------
Input::~Input() {
	if (!_knownPartition) {
		delete _knownPartition;
	}
	if (!_knownLabelDescription) {
		delete _knownLabelDescription;
	}
	for (unsigned int i = 0; i < _modelType.size(); i++) {
		delete _modelType[i];
		_modelType[i] = NULL;
	}
}


//------------
//------------
// Modifiers
//------------
//------------


//------ Criterion  ----//

//getCriterion[i]
//-------------------
CriterionName Input::getCriterionName(unsigned int index) const {
	if (index < _criterionName.size()) {
		return _criterionName[index];
	}
	else {
		THROW(InputException, wrongCriterionPositionInGet);
	}
}

// removeCriterionName
//--------------------
void Input::removeCriterion(unsigned int index) {
	if (index < _criterionName.size()) {
		_criterionName.erase(_criterionName.begin() + index);
	}
	else {
		THROW(InputException, wrongCriterionPositionInRemove);
	}
	_finalized = false;
}


//------ modelType  ----//

/// setModel
void Input::setModel(std::vector<ModelName> const & modelName) {
	// resize vector
	_modelType.resize(modelName.size());
	for (unsigned int iModel = 0; iModel < _modelType.size(); iModel++) {
		// copy vector contents
      if(_modelType[iModel]) delete _modelType[iModel];
      _modelType[iModel] = new ModelType(modelName[iModel]);
	}
}

//setModelType
//----------------
void Input::setModelType(const ModelType * modelType, unsigned int index) {
	if (index < _modelType.size()) {
		if (_modelType[index]) delete _modelType[index];
		_modelType[index] = new ModelType(*modelType);
	}
	else {  
		THROW(InputException, wrongModelPositionInSet);
	}
	_finalized = false;
}

// insertModelType
//-----------------
void Input::insertModelType(const ModelType * modelType, unsigned int index) {
	if (index >= 0 && index <= _modelType.size()) {
		_modelType.insert(_modelType.begin() + index, new ModelType(*modelType));
	}
	else {
		THROW(InputException, wrongModelPositionInInsert);
	}
	_finalized = false;
}


// add new model type 
void Input::addModelType(const ModelType * modelType) {
	if (getDataType() == QualitativeData)
  	if (getModelGenre(modelType->getModelName()) != QualitativeModel) return;


	if (getDataType() == QuantitativeData)
		if (getModelGenre(modelType->getModelName()) != QuantitativeModel) return;

	if (getDataType() == HeterogeneousData)
		if (getModelGenre(modelType->getModelName()) != HeterogeneousModel) return;

	bool found = false;
	for (unsigned int iModel = 0; iModel < _modelType.size(); iModel++) {
		if (_modelType[iModel]->getModelName() == modelType->getModelName()) found = true;
	}
	if (!found) _modelType.push_back(new ModelType(*modelType));
}

// add new model (modelName -> modelType)
void Input::addModel(ModelName const modelName) {

	if (getDataType() == QualitativeData)
		if (getModelGenre(modelName) != QualitativeModel) return;

	if (getDataType() == QuantitativeData)
		if (getModelGenre(modelName) != QuantitativeModel) return;

	if (getDataType() == HeterogeneousData)
		if (getModelGenre(modelName) != HeterogeneousModel) return;

	bool found = false;
	for (unsigned int iModel = 0; iModel < _modelType.size(); iModel++) {
		if (_modelType[iModel]->getModelName() == modelName) found = true;
	}
	if (!found) _modelType.push_back(new ModelType(modelName));
}

// removeModelType
//--------------------
void Input::removeModelType(unsigned int index) {
	if (index < _modelType.size()) {
		delete _modelType[index];
		_modelType.erase(_modelType.begin() + index);
	}
	else {
		THROW(InputException, wrongModelPositionInRemove);
	}
	_finalized = false;
}


// ---- Weight ----//

//setWeight()
//-----------
//TODO a enlever
void Input::setWeight(std::string weightFileName) {
	//_data->setWeight(weightFileName);
	_finalized = false;
}

void Input::setWeight(double* weight) {
	// get pointer to Data
	Data * data = getData();
	data->setWeight(weight);
	_finalized = false;
}

//TODO a enlever
void Input::insertWeight(std::string weightFileName) {
	//_data->setWeight(weightFileName);
	_finalized = false;
}

//TODO a enlever
void Input::removeWeight() {
	//_data->setWeightDefault();
	_finalized = false;
}


// ----- KnownPartition ----//

// setKnownPartition
//------------------
void Input::setKnownPartition(std::string iFileName) {
	if (_nbCluster.size() != 1) {
		THROW(InputException, badSetKnownPartition);
	}
	else {
		if (_knownPartition) {
			delete _knownPartition;
		}
		NumericPartitionFile partitionFile;
		partitionFile._fileName = iFileName;
		partitionFile._format = FormatNumeric::defaultFormatNumericFile;
		partitionFile._type = TypePartition::defaultTypePartition;
		_knownPartition = new Partition(_nbSample, _nbCluster[0], partitionFile);
	}
	_finalized = false;
}

// insertKnownPartition
//---------------------
void Input::insertKnownPartition(NumericPartitionFile partitionFile) {
	if (_nbCluster.size() != 1) {
		THROW(InputException, badSetKnownPartition);
	}
	else {
		if (_knownPartition) {
			delete _knownPartition;
		}
		_knownPartition = new Partition(_nbSample, _nbCluster[0], partitionFile);
	}
	_finalized = false;
}

// removeKnownPartition
//---------------------
void Input::removeKnownPartition() {
	if (_knownPartition) {
		delete _knownPartition;
	}
	_finalized = false;
}

void Input::setKnownLabelDescription(LabelDescription & labeldescription) {
	removeKnownLabelDescription();

	_knownLabelDescription = new LabelDescription(labeldescription);
}
void Input::setKnownLabelDescription(LabelDescription *labeldescription) {
	removeKnownLabelDescription();

	_knownLabelDescription = labeldescription;
}

void Input::removeKnownLabelDescription() {
	if (!_knownLabelDescription)
		delete _knownLabelDescription;

	_knownLabelDescription = NULL;
}

// --------
// finalize
// --------
void Input::finalize() {
	if (!_finalized) {
		_finalized = verif();
	}
}

// ----------------
// Verif
//-----------------
bool Input::verif() {
	bool res = true;

	// verif 1 : required inputs
	if (_nbSample == 0
		|| _pbDimension == 0
		|| _nbCluster.empty())
	{
		res = false;
		THROW(InputException, missingRequiredInputs);
	}

	// verif 2 : _data->verify
	res = _dataDescription.verifyData();

	return res;
}

// ----------------
// print out input
// ----------------
void Input::edit(std::ostream & out) const {
	out << "Models : ";
	for (unsigned int iModel = 0; iModel < _modelType.size(); iModel++)
		//out << endl << "  " << ModelNameToString(_modelType[iModel]->getModelName());
		out << endl <<*(_modelType[iModel]);
	out << std::endl;

	out << "Criterions : ";
	for (unsigned int iCriterion = 0; iCriterion < _criterionName.size(); iCriterion++)
		out << endl << "  " << CriterionNameToString(_criterionName[iCriterion]);
	out << std::endl;
}

}
