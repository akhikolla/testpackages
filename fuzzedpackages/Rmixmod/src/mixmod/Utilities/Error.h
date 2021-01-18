/***************************************************************************
                             SRC/mixmod/Utilities/Error.h  description
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
/**
  @file Error.h
  @brief Base class for Error(s)
  @author F Langrognet
 */

#ifndef XEMERROR_H
#define XEMERROR_H

#include "mixmod/Utilities/Util.h"

namespace XEM {

/**
  @brief Base class for Error(s)
  @author F Langrognet
 */
class Error {

public:

	/// Default constructor
	Error();

	/// Constructor
	Error(Exception & errorType);

	/// Destructor
	~Error();

	// setter
	inline void setError(Exception & errorType) {
		delete _errorType;
		_errorType = errorType.clone();
	}

	// getter
	inline Exception& getError() const {
		return *_errorType;
	}

	/// Run method (for debug)
	void run();

private:

	/// Type of error	
	Exception * _errorType;
};

}

#endif
