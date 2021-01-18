/***************************************************************************
                             SRC/mixmod/Kernel/IO/Description.h  description
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

#ifndef XEMDESCRIPTION_H
#define XEMDESCRIPTION_H

#include "mixmod/Utilities/Util.h"

namespace XEM {

// pre-declaration
class ColumnDescription;

/** 
 \class XEMDescription
 @author F. Langrognet
		@date 2011
		@brief XEDescription class
 */
class Description {

public:
	/// Default constructor
	Description();

	///constructor by initilization
	Description(int64_t nbSample, int64_t nbColumn, 
			std::vector<ColumnDescription *> columnDescription, 
			FormatNumeric::FormatNumericFile format, 
			std::string filename, std::string InfoName = "");

	///constructor by copy
	Description(Description & description);

	/// Destructor
	virtual ~Description();

	///initialization ColumnDescription by default
	void initializationColumnDescription();

	///operator=    
	Description & operator=(const Description & description);

	/// get
	//get Name
	std::string getInfoName()const;

	//get NbSample
	int64_t getNbSample() const;

	//get FileName
	std::string getFileName()const;

	//get NbColumn
	int64_t getNbColumn()const;

	//get Format
	FormatNumeric::FormatNumericFile getFormat() const;

	//get ColumnDescription
	const ColumnDescription * getColumnDescription(int64_t index)const;

	const std::vector<ColumnDescription*> getAllColumnDescription()const;

	int64_t getPbDimension() const;

	///set
	//set InfoName
	void setInfoName(std::string & iInfoName);

	virtual void saveNumericValues(std::string fileName = "") = 0;

protected:

	std::string _infoName;
	int64_t _nbSample;
	int64_t _nbColumn;
	std::string _fileName; //numeric file name
	FormatNumeric::FormatNumericFile _format; //format of  numeric file
	std::vector<ColumnDescription*> _columnDescription; //each variable has a description
};

inline std::string Description::getInfoName()const {
	return _infoName;
}

inline int64_t Description::getNbSample()const {
	return _nbSample;
}

inline int64_t Description::getNbColumn()const {
	return _nbColumn;
}

inline std::string Description::getFileName()const {
	return _fileName;
}

inline const ColumnDescription * Description::getColumnDescription(int64_t index)const {
	if (index >= 0 && index <= _nbColumn)
		return _columnDescription[index];
	else THROW(InputException, wrongIndexInGetMethod);
}

inline const std::vector<ColumnDescription*> Description::getAllColumnDescription()const {
	return _columnDescription;
}

inline FormatNumeric::FormatNumericFile Description::getFormat()const {
	return _format;
}

inline void Description::setInfoName(std::string & iInfoName) {
	_infoName = iInfoName;
}

}

#endif // XEMDESCRIPTION_H
