/***************************************************************************
                             SRC/mixmod/Utilities/exceptions/Exception.h  description
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
/** @file XEMException.h
 *  @brief Base class for MIXMOD Exceptions. All XEMException classes must inherit this class.
 *  @author Parmeet Bhatia
 **/

#ifndef XEM_EXCEPTION_H
#define XEM_EXCEPTION_H

#include <exception>
#include <string.h>
#include <map>
#include <iostream>
#include <typeinfo>
#include "mixmod/Utilities/exceptions/ErrorEnumerations.h"

namespace XEM {

class Exception : public std::exception {

public:

	Exception() {
	}
	Exception(const std::string& what_arg) throw ();
	Exception(const std::string& msg, const std::string& file, int line) throw ();
	Exception(const Exception & exception);
	virtual Exception* clone() throw ();

	/** Same as standard what() function.*/
	virtual const char* what() const throw ();

	/**Getter function for error.*/
	virtual Exception const& getError() const throw () {
		return *this;
	}
	/** Equality check : Return true if both the Exception class as well as it's ErrorType_ is same.*/
	virtual bool operator==(const Exception&) const throw ();

	/**Inequality check: Return true if equality check fails.*/
	virtual bool operator!=(const Exception& other) const throw () {
		return !(*this == other);
	}
	/**Interface for runner method. It will print the error and its content to the stream passed as
	 * argument.*/
	virtual void run(std::ostream & flux = std::cout) const throw ();

	virtual ~Exception() throw () {
	}

protected:
	
	std::string _errorMsg;
	std::string _filename;
	int _linenumber;
};

}

#endif /* XEMEXCEPTION_H_ */
