/***************************************************************************
                             SRC/mixmod/Utilities/exceptions/DCVonlyInGaussianCaseException.cpp  description
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
#include "mixmod/Utilities/exceptions/DCVonlyInGaussianCaseException.h"

namespace XEM {

// mapping: error type --> message translations
std::map<DCVonlyInGaussianCaseError, const char*> DCVonlyInGaussianCaseException::mapErrorMsg = DCVonlyInGaussianCaseException::create_map();
//DCVonlyInGaussianCaseException::mapErrorMsg.insert(std::make_pair(allCVCriterionErrorForAnEstimationInDCVContext, "All CV Criterion have error for an estimation in DCV context"));
//DCVonlyInGaussianCaseException::mapErrorMsg.insert(std::make_pair(NbDCVBlocksTooSmall, "Number of DCV Blocks is too small (or NbSample is too small)"));

DCVonlyInGaussianCaseException::DCVonlyInGaussianCaseException(std::string file, int line, DCVonlyInGaussianCaseError error) throw () {
	_errorType = error;
	_filename = file;
	_linenumber = line;
}

Exception * DCVonlyInGaussianCaseException::clone() throw () {
	return new DCVonlyInGaussianCaseException(*this);
}

DCVonlyInGaussianCaseException::DCVonlyInGaussianCaseException(DCVonlyInGaussianCaseError error) throw () {
	_errorType = error;
	_filename = "Defaulter";
	_linenumber = 0;
}

const char* DCVonlyInGaussianCaseException::what() const throw () {
	return (mapErrorMsg.find(_errorType))->second;
}

bool DCVonlyInGaussianCaseException::operator ==(const Exception& other) const throw () {
	if (typeid (*this) != typeid (other)) return false;

	return (*this)._errorType == dynamic_cast<const DCVonlyInGaussianCaseException&> (other)._errorType;
}

void DCVonlyInGaussianCaseException::run(std::ostream & flux) const throw () {
	flux << "In file: " << _filename << " at line: " << _linenumber << "\n";
	flux << "MIXMOD ERROR (DCVonlyInGaussianCase Error type " << _errorType << ") :" << "\n"; //new line will flush automatically
	flux << (*this).what() << "\n\n";
}

}
