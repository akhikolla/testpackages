/***************************************************************************
                             SRC/mixmod/Utilities/exceptions/DCVException.cpp  description
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
#include "mixmod/Utilities/exceptions/DCVException.h"

namespace XEM {

// mapping: error type --> message translations
std::map<DCVError, const char*> DCVException::mapErrorMsg = DCVException::create_map();
//DCVException::mapErrorMsg.insert(make_pair(wrongDCVinitBlocks, "DCV error : wrong init block specification, must be either RANDOM or DIAG"));
//DCVException::mapErrorMsg.insert(make_pair(wrongDCVnumberOfBlocks, "DCV error : wrong number of blocks, must be between 2 and the number of samples"));
//DCVException::mapErrorMsg.insert(make_pair(DCVmustBeDIAG, "DCV error : in this situation DCV init block specification must be DIAG"));
//DCVException::mapErrorMsg.insert(make_pair(forbiddenCallToGetBestCVModel, "DCV error : call to getBestCVModel is forbidden in the current context"));

DCVException::DCVException(std::string file, int line, DCVError error) throw () {
	_errorType = error;
	_filename = file;
	_linenumber = line;
}

DCVException::DCVException(DCVError error) throw () {
	_errorType = error;
	_filename = "Defaulter";
	_linenumber = 0;
}

const char* DCVException::what() const throw () {
	return (mapErrorMsg.find(_errorType))->second;
}

Exception * DCVException::clone() throw () {
	return new DCVException(*this);
}

bool DCVException::operator ==(const Exception& other) const throw () {
	if (typeid (*this) != typeid (other)) return false;

	return (*this)._errorType == dynamic_cast<const DCVException&> (other)._errorType;
}

void DCVException::run(std::ostream & flux) const throw () {
	flux << "In file: " << _filename << " at line: " << _linenumber << "\n";
	flux << "MIXMOD ERROR (DCV Error type " << _errorType << ") :" << "\n"; //new line will flush automatically
	flux << (*this).what() << "\n\n";
}

}
