/***************************************************************************
                             SRC/mixmod/Utilities/exceptions/OtherException.h  description
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
/** @file XEMOtherException.h
 *  @brief Exception class for Other types of error handling.
 *  @author Parmeet Bhatia
 **/

#ifndef XEM_OTHEREXCEPTION_H
#define XEM_OTHEREXCEPTION_H

#include "mixmod/Utilities/exceptions/Exception.h"

namespace XEM {

class OtherException : public Exception {

public:

	OtherException(std::string file, int line, OtherError error) throw ();
	OtherException(OtherError) throw ();
	Exception * clone() throw ();
	virtual const char* what() const throw ();
	virtual bool operator==(const Exception&) const throw ();
	virtual void run(std::ostream & flux = std::cout) const throw ();

	virtual ~OtherException() throw () {
	}


	static std::map<OtherError, const char*> create_map()
	{
		std::map<OtherError, const char*> m;

		m.insert(std::make_pair(badFormat, "Bad Format"));
		m.insert(std::make_pair(nullPointerError, "Internal error (Null pointer)"));
		m.insert(std::make_pair(wrongMatrixType, "Error : trying to apply a method on a wrong matrix type "));
		m.insert(std::make_pair(wrongConstructorType, "Error : when constructing an object by default"));
		m.insert(std::make_pair(nonImplementedMethod, "non implemented method"));
		m.insert(std::make_pair(badBinaryParameterClass, "Internal Mixmod Error: bad XEMBinrayParameter Class"));
		m.insert(std::make_pair(UnknownReason, "Error occurred due to unknown reason."));
		m.insert(std::make_pair(internalMixmodError, "Internal error in mixmod sofware"));
		m.insert(std::make_pair(FunctionNotYetImplemented, "Function that is called is not yet implemented"));
		m.insert(std::make_pair(AllModelsGotErros, "All models got errors"));
		m.insert(std::make_pair(AllTriesGotErros, "All tries got errors"));
		m.insert(std::make_pair(xmlFeaturesNotAvailable, "XML features are not available"));
		return m;
	}
	static std::map<OtherError, const char*> mapErrorMsg;
    OtherError getErrorType() const throw();
    
protected:

	OtherError _errorType;
};

 inline OtherError OtherException::getErrorType() const throw() {return _errorType;}
}

#endif /* XEMOTHEREXCEPTION_H_ */
