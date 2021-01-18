/***************************************************************************
                             SRC/mixmod/Kernel/Criterion/CriterionOutput.cpp  description
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

#include "mixmod/Kernel/Criterion/CriterionOutput.h"

namespace XEM {

//------------
// Constructor
//------------
CriterionOutput::CriterionOutput() {
	_criterionName = UNKNOWN_CRITERION_NAME;
	_value = 0;
	_error = NOERROR.clone();
}

//------------
// Copy Constructor
//------------
CriterionOutput::CriterionOutput(const CriterionOutput& criterionOutput) {
	_value = criterionOutput._value;
	_criterionName = criterionOutput._criterionName;
	_error = criterionOutput._error->clone();
}

//------------
// Constructor
//------------
CriterionOutput::CriterionOutput(CriterionName criterionName) {
	_value = 0;
	_criterionName = criterionName;
	_error = NOERROR.clone();
}

//------------
// Constructor
//------------
CriterionOutput::CriterionOutput(CriterionName criterionName, 
		double criterionValue, Exception& criterionErrorType) 
{
	_value = criterionValue;
	_criterionName = criterionName;
	//_error         = &criterionErrorType;
	//changed to clone instead of reference
	_error = (criterionErrorType).clone();
}

//-----------
// Destructor
//-----------
CriterionOutput::~CriterionOutput() {
	if (_error) {
		delete _error;
		_error = NULL;
	}
}

/// Comparison operator
bool CriterionOutput::operator ==(const CriterionOutput & criterionOutput) const {
	if (_value != criterionOutput.getValue()) return false;
	if (_criterionName != criterionOutput.getCriterionName()) return false;
	if (*(dynamic_cast<Exception*> (_error)) != (criterionOutput.getError())) return false;
	return true;
}

//----------
// edit Type
//----------
void CriterionOutput::editType(std::ostream & oFile) const {
	oFile << "Criterion Name : ";
	if (_criterionName == BIC) {
		oFile << "BIC";
	}
	else if (_criterionName == CV) {
		oFile << "CV";
	}
	else if (_criterionName == DCV) {
		oFile << "DCV";
	}
	else if (_criterionName == NEC) {
		oFile << "NEC";
	}
	else if (_criterionName == ICL) {
		oFile << "ICL";
	}
	oFile << endl << "---------------" << endl << endl;
}

//-----------
// edit Value
//-----------
void CriterionOutput::editValue(std::ostream & oFile, bool text) const {
	if (text) {
		oFile << "\t\t\tCriterion Value : ";
		if (*(dynamic_cast<Exception*> (_error)) == NOERROR) {
			oFile << _value << endl << endl;
		}
		else {
			oFile << "numeric Error" << endl << endl;
		}
	}
	else {
		if (*(dynamic_cast<Exception*> (_error)) == NOERROR) {
			oFile << _value << endl << endl;
		}
	}
}

//--------------------
// edit Type And Value
//--------------------
void CriterionOutput::editTypeAndValue(std::ostream & oFile) const {
	if (_criterionName == BIC) {
		oFile << "\t\t\tBIC ";
	}
	else if (_criterionName == CV) {
		oFile << "\t\t\tCV ";
	}
	else if (_criterionName == DCV) {
		oFile << "\t\t\tDCV ";
	}
	else if (_criterionName == NEC) {
		oFile << "\t\t\tNEC ";
	}
	else if (_criterionName == ICL) {
		oFile << "\t\t\tICL ";
	}

	oFile << "Criterion Value : ";
	if (*(dynamic_cast<Exception*> (_error)) == NOERROR) {
		oFile << _value << endl << endl;
	}
	else {
		//cout << "----<<<<<<<>>>>>>>>>>" << << endl;
		oFile << "numeric Error" << endl << endl;
	}
}

}
