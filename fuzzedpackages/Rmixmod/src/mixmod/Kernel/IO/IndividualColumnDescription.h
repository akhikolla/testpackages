/***************************************************************************
                             SRC/mixmod/Kernel/IO/IndividualColumnDescription.h  description
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
#ifndef XEMINDIVIDUALCOLUMNDESCRIPTION_H
#define XEMINDIVIDUALCOLUMNDESCRIPTION_H

#include "mixmod/Kernel/IO/ColumnDescription.h"

namespace XEM {

//Individual Description

struct IndividualDescription {

	int64_t num;
	std::string name;
};

/** 
 \class XEMIndividualColumnDescription
 @author F. Langrognet
		@date 2011
		@brief XEMIndividualColumnDescription class derived from XEMColumnDescription
 */
class IndividualColumnDescription : public ColumnDescription {

public:
	
	/// Default constructor
	IndividualColumnDescription();

	/// initialization constructor
	IndividualColumnDescription(int64_t index);

	/// Destructor
	virtual ~IndividualColumnDescription();

	std::string editType();

	ColumnDescription * clone()const;

	const std::vector<IndividualDescription> & getIndividualDescription()const;

	void setIndividualDescription(IndividualDescription & individualDescription, unsigned int index);

	void insertIndividualDescription(IndividualDescription individualDescription, unsigned int index);

private:
	
	std::vector<IndividualDescription> _individualDescription;
};

inline const std::vector<IndividualDescription> & IndividualColumnDescription::getIndividualDescription()const {
	return _individualDescription;
}

}

#endif // XEMINDIVIDUALCOLUMNDESCRIPTION_H
