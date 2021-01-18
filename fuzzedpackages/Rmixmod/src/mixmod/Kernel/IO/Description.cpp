/***************************************************************************
                             SRC/mixmod/Kernel/IO/Description.cpp  description
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

#include "mixmod/Kernel/IO/Description.h"
#include "mixmod/Kernel/IO/ColumnDescription.h"
#include "mixmod/Kernel/IO/WeightColumnDescription.h"
#include "mixmod/Kernel/IO/IndividualColumnDescription.h"
#include "mixmod/Kernel/IO/QuantitativeColumnDescription.h"

namespace XEM {

//------------
// Constructor by default
//------------
Description::Description() {
	_fileName = "";
	_format = FormatNumeric::defaultFormatNumericFile;
	_infoName = "";
	_nbSample = 0;
	_nbColumn = 0;
	initializationColumnDescription();
}

//------------
// Constructor by initialization
//------------
Description::Description(int64_t nbSample, int64_t nbColumn, 
		std::vector<ColumnDescription *> columnDescription, 
		FormatNumeric::FormatNumericFile format, std::string filename, std::string infoName)
{
	_fileName = filename;
	_format = format;
	_infoName = infoName;
	_nbSample = nbSample;
	_nbColumn = nbColumn;
	const unsigned int columnSize = columnDescription.size();
	if (columnSize != _nbColumn) {
		THROW(InputException, errorInColumnDescription);
	}
	_columnDescription.resize(_nbColumn);
	for (int64_t i = 0; i < _nbColumn; ++i) {
		_columnDescription[i] = columnDescription[i]->clone();
	}
}

//------------
// Constructor by copy
//------------
Description::Description(Description & description) {
	*this = description;
}

//------------
// operator =
//------------
Description & Description::operator=(const Description & description) {
	_fileName = description._fileName;
	_format = description._format;
	_infoName = description._infoName;
	_nbSample = description._nbSample;
	_nbColumn = description._nbColumn;
	_columnDescription.resize(_nbColumn);
	for (int64_t i = 0; i < _nbColumn; ++i) {
		const ColumnDescription * cd = description.getColumnDescription(i);
		_columnDescription[i] = cd->clone();
	}
	return *this;
}

//------------
// Destructor
//------------
Description::~Description() {
	if (_columnDescription.size() != 0) {
		for (unsigned int i = 0; i < _columnDescription.size(); ++i) {
			delete _columnDescription[i];
		}
	}
}

void Description::initializationColumnDescription() {
	_columnDescription.resize(_nbColumn);
	for (int64_t i = 0; i < _nbColumn; ++i) {
		//auto_ptr<XEMQuantitativeColumnDescription> a (new XEMQuantitativeColumnDescription());
		//_columnDescription[i] = a;
		_columnDescription[i] = new QuantitativeColumnDescription(i);
	}
}

//------------------------------------------
// return pbDimension (number of variables)
//
int64_t Description::getPbDimension() const {
	int64_t nbVariable = _nbColumn;
	for (int64_t i = 0; i < _nbColumn; i++) {
		if (typeid (*(_columnDescription[i])) == typeid (IndividualColumnDescription)) {
			nbVariable--;
		}
		if (typeid (*(_columnDescription[i])) == typeid (WeightColumnDescription)) {
			nbVariable--;
		}
	}
	return nbVariable;
}

}
