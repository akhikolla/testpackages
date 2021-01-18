/***************************************************************************
                             SRC/mixmod/Kernel/IO/DataDescription.h  description
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

#ifndef XEMDATADESCRIPTION_H
#define XEMDATADESCRIPTION_H

#include "mixmod/Kernel/IO/Description.h"

namespace XEM {

// pre-declaration
class Data;
class GaussianData;
class BinaryData;
class CompositeData;

/** 
 \class XEMDataDescription
 @author F. Langrognet
		@date 2011
		@brief XEMDataDescription class derived from XEMDescription
 */
class DataDescription : public Description {

public:
	
	/// Default constructor
	DataDescription();

	///constructor by initilization
	DataDescription(int64_t nbSample, int64_t nbColumn, 
			std::vector<ColumnDescription *> columnDescription, 
			FormatNumeric::FormatNumericFile format, 
			std::string filename, std::string infoName = "");

	/// constructor with a gaussianData
	DataDescription(GaussianData * gData);

	/// constructor with a binaryData
	DataDescription(BinaryData * bData);

	/// constructor with composite data
	DataDescription(CompositeData * cData);

	///constructor by copy
	DataDescription(DataDescription & dataDescription);

	/// Destructor
	virtual ~DataDescription();

	///operator=    
	DataDescription & operator=(const DataDescription & dataDescription);

	Data * getData() const;

	/// is binary Data now changed to GetXEMDataType
	DataType getDataType() const;

	/// verify data validity
	bool verifyData() const;

	void saveNumericValues(std::string fileName = "");

	void releaseData();
	
private:

	Data * _data;

	/// Create and return XEMData *
	Data * createData() const;
};

inline Data * DataDescription::getData() const {
	return _data;
}

inline void DataDescription::saveNumericValues(std::string fileName) {
	THROW(OtherException, internalMixmodError);
}

inline void DataDescription::releaseData(){
	_data = NULL;
}

}

#endif // XEMDATADESCRIPTION_H
