/***************************************************************************
                             SRC/mixmod/Utilities/exceptions/DCVException.h  description
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
/** @file XEMDCVException.h
 *  @brief Exception class for DCV errors.
 *  @author Parmeet Bhatia
 **/

#ifndef XEM_DCVEXCEPTION_H
#define XEM_DCVEXCEPTION_H

#include "mixmod/Utilities/exceptions/Exception.h"

namespace XEM {

class DCVException : public Exception {

public:
	
	DCVException(std::string, int, DCVError) throw ();
	Exception * clone() throw ();
	DCVException(DCVError) throw ();
	virtual const char* what() const throw ();
	virtual bool operator==(const Exception&) const throw ();
	virtual void run(std::ostream & flux = std::cout) const throw ();

	virtual ~DCVException() throw () {
	}

  static std::map<DCVError, const char*> create_map()
  {
    std::map<DCVError, const char*> m;
    m.insert(std::make_pair(wrongDCVinitBlocks, "DCV error : wrong init block specification, must be either RANDOM or DIAG"));
    m.insert(std::make_pair(wrongDCVnumberOfBlocks, "DCV error : wrong number of blocks, must be between 2 and the number of samples"));
    m.insert(std::make_pair(DCVmustBeDIAG, "DCV error : in this situation DCV init block specification must be DIAG"));
    m.insert(std::make_pair(forbiddenCallToGetBestCVModel, "DCV error : call to getBestCVModel is forbidden in the current context"));

    return m;
  }

	static std::map<DCVError, const char*> mapErrorMsg;

protected:
	
	DCVError _errorType;
};

}

#endif /* XEMDCVEXCEPTION_H_ */
