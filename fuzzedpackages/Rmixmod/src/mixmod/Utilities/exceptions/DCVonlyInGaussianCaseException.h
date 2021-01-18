/***************************************************************************
                             SRC/mixmod/Utilities/exceptions/DCVonlyInGaussianCaseException.h  description
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
/** @file XEMDCVonlyInGaussianCaseException.h
 *  @brief Exception class for DCV Errors in Gaussian case only.
 *  @author Parmeet Bhatia
 **/

#ifndef XEM_DCVONLYINGAUSSIANCASEEXCEPTION_H
#define XEM_DCVONLYINGAUSSIANCASEEXCEPTION_H

#include "mixmod/Utilities/exceptions/Exception.h"

namespace XEM {

class DCVonlyInGaussianCaseException : public Exception {

public:
	
	Exception * clone() throw ();
	DCVonlyInGaussianCaseException(std::string, int, DCVonlyInGaussianCaseError) throw ();
	DCVonlyInGaussianCaseException(DCVonlyInGaussianCaseError) throw ();
	virtual const char* what() const throw ();
	virtual bool operator==(const Exception&) const throw ();
	virtual void run(std::ostream & flux = std::cout) const throw ();

	virtual ~DCVonlyInGaussianCaseException() throw () {
	}

	static std::map<DCVonlyInGaussianCaseError, const char*> create_map()
	{
		std::map<DCVonlyInGaussianCaseError, const char*> m;
		m.insert(std::make_pair(allCVCriterionErrorForAnEstimationInDCVContext, "All CV Criterion have error for an estimation in DCV context"));
		m.insert(std::make_pair(NbDCVBlocksTooSmall, "Number of DCV Blocks is too small (or NbSample is too small)"));

		return m;
	}
	static std::map<DCVonlyInGaussianCaseError, const char*> mapErrorMsg;

protected:
	
	DCVonlyInGaussianCaseError _errorType;
};

}

#endif /* XEMDCVONLYINGAUSSIANCASEEXCEPTION_H_ */
