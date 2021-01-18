/***************************************************************************
                             SRC/mixmod/Utilities/exceptions/Exception.cpp  description
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
#include "mixmod/Utilities/exceptions/Exception.h"

namespace XEM {

Exception::Exception(const std::string& msg, const std::string& file, int line) throw () {
	_errorMsg = msg;
	_filename = file;
	_linenumber = line;
}

Exception::Exception(const std::string& msg) throw () {
	_errorMsg = msg;
	_filename = "Defaulter";
	_linenumber = 0;
}

Exception::Exception(const Exception & exception) {
	_errorMsg = exception._errorMsg;
	_filename = exception._filename;
	_linenumber = exception._linenumber;
}

const char* Exception::what() const throw () {
	return _errorMsg.c_str();
}

Exception * Exception::clone() throw () {
	return new Exception(*this);
}

bool Exception::operator ==(const Exception& other) const throw () {
	if (typeid (*this) != typeid (other)) return false;
	if (strcmp(((*this)._errorMsg).c_str(), (other._errorMsg).c_str()) == 0)
		return true;
	else
		return false;
}

void Exception::run(std::ostream & flux) const throw () {
	flux << (*this).what() << "\n\n";
}

}
