/***************************************************************************
                             SRC/mixmod/Kernel/Criterion/CriterionOutput.h  description
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
#ifndef XEMCriterionOutput_H
#define XEMCriterionOutput_H

#include "mixmod/Utilities/Util.h"

namespace XEM {

/** @brief Base class for Label(s)
	@author F Langrognet
 */

class CriterionOutput {

public:

	/// Default constructor
	CriterionOutput();

	/// Copy constructor
	CriterionOutput(const CriterionOutput& criterionOutput);

	/// constructor
	CriterionOutput(CriterionName criterionName);

	/// Constructor
	CriterionOutput(CriterionName criterionName, double criterionValue, Exception& criterionErrorType);

	/// Destructor
	virtual ~CriterionOutput();

	/// Comparison operator
	bool operator ==(const CriterionOutput & criterionOutput) const;

	///editType
	void editType(std::ostream & oFile) const;

	///editValue
	void editValue(std::ostream & oFile, bool text = false) const;

	/// editTypeAndValue
	void editTypeAndValue(std::ostream & oFile) const;

	//--------------
	// get
	//--------------
	/// getCriterionName
	CriterionName const getCriterionName() const;

	/// getValue
	double const getValue() const;

	/// getError
	Exception & getError() const;

	//--------------
	// set
	//--------------
	/// setCriterionName
	void setCriterionName(CriterionName criterionName);

	/// setValue
	void setValue(double value);

	/// setError
	void setError(Exception& e);

private:
	
	/// Criterion value
	double _value;

	/// Error type in calculation of criterion value
	Exception * _error;

	/// criterion name
	CriterionName _criterionName;
};

inline CriterionName const CriterionOutput::getCriterionName() const {
	return _criterionName;
}

inline double const CriterionOutput::getValue() const {
	return _value;
}

inline Exception & CriterionOutput::getError() const {
	return *_error;
}

inline void CriterionOutput::setCriterionName(CriterionName criterionName) {
	_criterionName = criterionName;
}

inline void CriterionOutput::setValue(double value) {
	_value = value;
}

inline void CriterionOutput::setError(Exception& e) {
	if (_error) delete _error;
	_error = e.clone();
}

}

#endif
