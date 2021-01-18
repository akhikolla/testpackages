/***************************************************************************
                             SRC/mixmod/Utilities/Error.cpp  description
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
#include "mixmod/Utilities/Error.h"

namespace XEM {

//-----------
//Constructor
//-----------
Error::Error() : _errorType(NOERROR.clone()) {
}

//-----------
//Constructor
//-----------
Error::Error(Exception& errorType) : _errorType((((errorType)).clone())) {
}

//----------
//Destructor
//----------
Error::~Error() {
	if (_errorType) delete _errorType;
}

//---
//Run
//---
void Error::run() {
	((Exception*) _errorType)->run();
}

}
