/***************************************************************************
                             SRC/mixmod/Kernel/IO/DataDescription.cpp  description
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

#include "mixmod/Kernel/IO/DataDescription.h"
#include "mixmod/Kernel/IO/BinaryData.h"
#include "mixmod/Kernel/IO/GaussianData.h"
#include "mixmod/Kernel/IO/CompositeData.h"
#include "mixmod/Kernel/IO/WeightColumnDescription.h"
#include "mixmod/Kernel/IO/QuantitativeColumnDescription.h"
#include "mixmod/Kernel/IO/QualitativeColumnDescription.h"
#include "mixmod/Kernel/IO/UnusedColumnDescription.h"

namespace XEM {

//------------
// Constructor by default
//------------
DataDescription::DataDescription() : Description() {
	_data = NULL;
}

//------------
// Constructor by initialization
//------------
DataDescription::DataDescription(int64_t nbSample, int64_t nbColumn, 
		std::vector<ColumnDescription*> columnDescription,
		FormatNumeric::FormatNumericFile format, std::string filename, std::string infoName) 
: Description(nbSample, nbColumn, columnDescription, format, filename, infoName) 
{
	_data = createData();
}

//------------
// Constructor
//------------
DataDescription::DataDescription(GaussianData * gData) {
	_fileName = "";
	_format = FormatNumeric::defaultFormatNumericFile;
	_infoName = "";
	_nbSample = gData->getNbSample();
	_nbColumn = gData->getPbDimension();
	_columnDescription.resize(_nbColumn);
	for (int64_t i = 0; i < _nbColumn; ++i) {
		_columnDescription[i] = new QuantitativeColumnDescription(i);
	}
	_data = gData->clone();

	if (!_data->hasDefaultWeight()) {
		_columnDescription.push_back(new WeightColumnDescription(_nbColumn));
	}
}

//------------
// Constructor
//------------
DataDescription::DataDescription(BinaryData * bData) {
	_fileName = "";
	_format = FormatNumeric::defaultFormatNumericFile;
	_infoName = "";
	_nbSample = bData->getNbSample();
	_nbColumn = bData->getPbDimension();
	_columnDescription.resize(_nbColumn);
	int64_t * tabModality = bData->getTabNbModality();
	for (int64_t i = 0; i < _nbColumn; ++i) {
		_columnDescription[i] = new QualitativeColumnDescription(i, tabModality[i]);
	}
	_data = bData->clone();

	if (!_data->hasDefaultWeight()) {
		_columnDescription.push_back(new WeightColumnDescription(_nbColumn));
	}
}

// constructor for XEMCompositeData
DataDescription::DataDescription(CompositeData * cData) {
	BinaryData * bData = cData->getBinaryData();
	GaussianData * gData = cData->getGaussianData();
	assert(bData != NULL);
	assert(gData != NULL);
	_fileName = "";
	_format = FormatNumeric::defaultFormatNumericFile;
	_infoName = "";
	_nbSample = cData->getNbSample();
	_nbColumn = cData->getPbDimension();
	_columnDescription.resize(_nbColumn);
	int64_t * tabModality = bData->getTabNbModality();

	//Column 0:bData->getPbDimension() is Binary Data and from
	//Column bData->getPbDimension()-1:_nbColumn is Gaussian Data
	for (int64_t i = 0; i < bData->getPbDimension(); ++i) {
		_columnDescription[i] = new QualitativeColumnDescription(i, tabModality[i]);
	}

	for (int64_t i = bData->getPbDimension(); i < _nbColumn; ++i) {
		_columnDescription[i] = new QuantitativeColumnDescription(i);
	}
	_data = cData->clone();
	if (!_data->hasDefaultWeight()) {
		_columnDescription.push_back(new WeightColumnDescription(_nbColumn));
	}
}

//------------
// Constructor by copy
//------------
DataDescription::DataDescription(DataDescription & dataDescription) {
	*this = dataDescription;
}

//------------
// operator =
//------------
DataDescription & DataDescription::operator=(const DataDescription & dataDescription) {
	_fileName = dataDescription._fileName;
	_format = dataDescription._format;
	_infoName = dataDescription._infoName;
	_nbSample = dataDescription._nbSample;
	_nbColumn = dataDescription._nbColumn;
	const Data * data = dataDescription.getData();
	if (data) {
		_data = data->clone();
	}
	else {
		_data = NULL;
	}
	_columnDescription.resize(_nbColumn);
	for (int64_t i = 0; i < _nbColumn; ++i) {
		const ColumnDescription * cd = dataDescription.getColumnDescription(i);
		_columnDescription[i] = cd->clone();
	}
	return *this;
}

//------------
// Destructor
//------------
DataDescription::~DataDescription() {
	if (_data) {
		delete _data;
	}
}

//----------------------------
// createData
//----------------------------
// Create and return XEMData *
Data* DataDescription::createData() const {

	Data* data = NULL;
	std::vector<int64_t> nbModality;
	int64_t nbQualitativeVariable = 0;
	int64_t nbQuantitativeVariable = 0;
	bool weightColumn = false;

	for (std::vector<ColumnDescription*>::const_iterator it=_columnDescription.begin(); it != _columnDescription.end(); it++) {
		ColumnDescription* columnDescription = *it;
		if (typeid (*columnDescription) == typeid (QualitativeColumnDescription)) {
			nbQualitativeVariable++;
			QualitativeColumnDescription* qualitativeColumnDescription =
				dynamic_cast<QualitativeColumnDescription*> (columnDescription);
			nbModality.push_back(qualitativeColumnDescription->getNbFactor());		
		}
		else if (typeid (*columnDescription) == typeid (QuantitativeColumnDescription))
			nbQuantitativeVariable++;
		else if (typeid (*columnDescription) == typeid (WeightColumnDescription)) {
			if (weightColumn)
				THROW(InputException, tooManyWeightColumnDescription);
			weightColumn = true;
		}
		else if (typeid (*columnDescription) == typeid (UnusedColumnDescription)) {
			//nothing to do: just to not forget this type
		}
	}

	if (nbQualitativeVariable == 0 && nbQuantitativeVariable == 0)
		THROW(InputException, badDataDescription);

	else if (nbQuantitativeVariable != 0 && nbQualitativeVariable == 0)
		data = new GaussianData(_nbSample, nbQuantitativeVariable);

	else if (nbQuantitativeVariable == 0 && nbQualitativeVariable != 0)
		data = new BinaryData(_nbSample, nbQualitativeVariable, nbModality);

	else {
		// nbQuantitativeVariable!=0 and nbQualitativeVariable!=0
		GaussianData* gData = new GaussianData(_nbSample, nbQuantitativeVariable);
		BinaryData* bData = new BinaryData(_nbSample, nbQualitativeVariable, nbModality);
		data = new CompositeData(gData, bData);
	}
	data->input(*this);

	return data;
}

DataType DataDescription::getDataType() const {
	int64_t nbQualitativeVariable = 0;
	int64_t nbQuantitativeVariable = 0;
	for (int64_t i = 0; i < _nbColumn; i++) {
		if (typeid (*(_columnDescription[i])) == typeid (QualitativeColumnDescription))
			nbQualitativeVariable++;
		if (typeid (*(_columnDescription[i])) == typeid (QuantitativeColumnDescription))
			nbQuantitativeVariable++;
	}

	if (nbQualitativeVariable == 0)
		return QuantitativeData;

	if (nbQuantitativeVariable == 0)
		return QualitativeData;

	return HeterogeneousData;
}

//---------------------
// verify data validity
//---------------------
bool DataDescription::verifyData() const {
	return _data->verify();
}

}
