/***************************************************************************
                             SRC/mixmod/Utilities/exceptions/InputException.cpp  description
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
#include "mixmod/Utilities/exceptions/InputException.h"

namespace XEM {

// mapping: error type --> message translations
std::map<InputError, const char*> InputException::mapErrorMsg = InputException::create_map();

InputException::InputException(std::string file, int line, InputError error) throw () {
	_errorType = error;
	_filename = file;
	_linenumber = line;
}

InputException::InputException(InputError error) throw () {
	_errorType = error;
	_filename = "Defaulter";
	_linenumber = 0;
}

Exception * InputException::clone() throw () {
	return new InputException(*this);
}

InputException::InputException(const InputException & inputException) : Exception(inputException) {
	_errorType = inputException._errorType;
}

const char* InputException::what() const throw () {
	return (mapErrorMsg.find(_errorType))->second;
}

bool InputException::operator ==(const Exception& other) const throw () {
	if (typeid (*this) != typeid (other)) return false;

	return (*this)._errorType == dynamic_cast<const InputException&> (other)._errorType;
}

void InputException::run(std::ostream & flux) const throw () {
	flux << "In file: " << _filename << " at line: " << _linenumber << "\n";
	flux << "MIXMOD ERROR (Input Error type " << _errorType << ") :" << "\n"; //new line will flush automatically
	flux << (*this).what() << "\n\n";
}

}
