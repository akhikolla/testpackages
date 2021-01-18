/***************************************************************************
                             SRC/mixmod/Kernel/IO/ProbaDescription.h  description
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

#ifndef XEMPROBADESCRIPTION_H
#define XEMPROBADESCRIPTION_H

#include "mixmod/Kernel/IO/Description.h"

namespace XEM {

// pre-declaration
class Model;
class Proba;

/** 
 \class XEMProbaDescription
 @author F. Langrognet
		@date 2011
		@brief XEMProbaDescription class derived from XEMDescription
 */
class ProbaDescription : public Description {

public:
	
	/// Default constructor
	ProbaDescription();

	///constructor by initilization
	ProbaDescription(int64_t nbSample, int64_t nbCluster, FormatNumeric::FormatNumericFile format, std::string filename, std::string infoName = "");

	///constructor after an estimation->run
	ProbaDescription(Model * model);


	///constructor by copy
	ProbaDescription(ProbaDescription & probaDescription);

	/// Destructor
	virtual ~ProbaDescription();


	/// Comparison operator
	bool operator==(ProbaDescription & probaDescription) const;

	///operator=    
	ProbaDescription & operator=(ProbaDescription & probaDescription);

	Proba * getProba();

	void saveNumericValues(std::string fileName = "");

	//void editProba(ostream & f) const;

private:

	Proba * _proba;
};

inline Proba * ProbaDescription::getProba() {
	return _proba;
}

}

#endif // XEMPROBADESCRIPTION_H
