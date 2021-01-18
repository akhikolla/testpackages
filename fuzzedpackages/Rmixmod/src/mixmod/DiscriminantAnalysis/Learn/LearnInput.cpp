/***************************************************************************
                             SRC/mixmod/DiscriminantAnalysis/Learn/LearnInput.cpp  description
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
#include "mixmod/DiscriminantAnalysis/Learn/LearnInput.h"
#include "mixmod/Kernel/IO/DataDescription.h"
#include "mixmod/Kernel/IO/LabelDescription.h"

namespace XEM {

//--------------------
// Default Constructor
//--------------------
LearnInput::LearnInput() : _nbCVBlock(defaultCVnumberOfBlocks) {
}

//-----------------
//  Copy constructor
//-----------------
LearnInput::LearnInput(const LearnInput & cInput)
: Input(cInput), _nbCVBlock(cInput.getNbCVBlock()) {
	// TODO clone strategy
}

//---------------------------
// Initialisation Constructor
//---------------------------
LearnInput::LearnInput(DataDescription * learnData, LabelDescription * knownLabelDescription)
: Input(std::vector<int64_t>(1, knownLabelDescription->getNbCluster()), *learnData) 
{
	// set partition
	setKnownLabelDescription(*knownLabelDescription);
	// set CV as default criterion 
	setCriterion(defaultLearnCriterionName, 0);
	// set the number of CV blocks
	_nbCVBlock = defaultCVnumberOfBlocks;
}

//-----------
// Destructor
//-----------
LearnInput::~LearnInput() {
}

//setCriterion
//----------------
void LearnInput::setCriterion(std::vector<CriterionName> const & criterionName) {

	// copy vector contents
	_criterionName = criterionName;

	// check vector contents
	for (unsigned int iCriterion = 0; iCriterion < _criterionName.size(); iCriterion++) {
		switch (_criterionName[iCriterion]) {
		case BIC: break;
		case CV: break;
		case ICL: THROW(InputException, badCriterion);
		case NEC: THROW(InputException, badCriterion);
		case UNKNOWN_CRITERION_NAME: THROW(OtherException, internalMixmodError);
			break;
		default: THROW(OtherException, internalMixmodError);
		}
	}
	_finalized = false;
}

//setCriterionName
//----------------
void LearnInput::setCriterion(const CriterionName criterionName, unsigned int index) {
	if (index < _criterionName.size()) {
		switch (criterionName) {
		case BIC: _criterionName[index] = BIC;
			break;
		case CV: _criterionName[index] = CV;
			break;
		case ICL: THROW(InputException, badCriterion);
		case NEC: THROW(InputException, badCriterion);
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

// insertCriterionName
//-----------------
void LearnInput::insertCriterion(const CriterionName criterionName, unsigned int index) {
	if (index <= _criterionName.size()) {
		switch (criterionName) {
		case BIC: _criterionName.insert(_criterionName.begin() + index, BIC);
			break;
		case CV: _criterionName.insert(_criterionName.begin() + index, CV);
			break;
		case ICL: THROW(InputException, badCriterion);
		case NEC: THROW(InputException, badCriterion);
		case UNKNOWN_CRITERION_NAME: THROW(OtherException, internalMixmodError);
			break;
		default: THROW(OtherException, internalMixmodError);
		}
	}
	else {
		THROW(InputException, wrongCriterionPositionInInsert);
	}
	_finalized = false;
}

// add Criterion
//-----------------
void LearnInput::addCriterion(const CriterionName criterionName) {

	bool found = false;
	for (unsigned int iCriterion = 0; iCriterion < _criterionName.size(); iCriterion++) {
		if (_criterionName[iCriterion] == criterionName) found = true;
	}
	if (!found) {
		switch (criterionName) {
		case BIC: _criterionName.push_back(BIC);
			break;
		case CV: _criterionName.push_back(CV);
			break;
		case ICL: THROW(InputException, badCriterion);
		case NEC: THROW(InputException, badCriterion);
		case UNKNOWN_CRITERION_NAME: THROW(OtherException, internalMixmodError);
			break;
		default: THROW(OtherException, internalMixmodError);
		}
	}
	_finalized = false;
}

// ---------------------------
// set the number of CV blocks
// ---------------------------
void LearnInput::setNbCVBlock(int64_t nbCVBlock) {
	_nbCVBlock = nbCVBlock;
}

// ----------------
// Verif
//-----------------
bool LearnInput::verif() {
	bool res = Input::verif();

	return res;
}

}
