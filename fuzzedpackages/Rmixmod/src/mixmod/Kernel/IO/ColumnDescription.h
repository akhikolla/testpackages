/***************************************************************************
                             SRC/mixmod/Kernel/IO/ColumnDescription.h  description
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
#ifndef XEMCOLUMNDESCRIPTION_H
#define XEMCOLUMNDESCRIPTION_H

#include "mixmod/Utilities/Util.h"

namespace XEM {

/** 
 \class XEMColumnDescription
 @author F. Langrognet
		@date 2011
		@brief XEMColumnDescription
 */
class ColumnDescription {

public:
	
	/// Default constructor
	ColumnDescription();

	/// Initialization constructor
	ColumnDescription(int64_t index);

	/// Destructor
	virtual ~ColumnDescription();

	virtual std::string editType() = 0;

	virtual ColumnDescription * clone()const = 0;

	///selector

	///get index of column
	const int64_t & getIndex()const;

	///get name of column 
	const std::string & getName()const;

	///set name of column
	void setName(std::string & strName);

protected:

	///index of column (0 to XEMDataDescription::_nbColumn-1)
	int64_t _index;

	///name of column (optional)
	std::string _name;
};

inline const int64_t & ColumnDescription::getIndex()const {
	return _index;
}

inline const std::string & ColumnDescription::getName()const {
	return _name;
}

inline void ColumnDescription::setName(std::string & strName) {
	_name = strName;
}

}

#endif // XEMCOLUMNDESCRIPTION_H
