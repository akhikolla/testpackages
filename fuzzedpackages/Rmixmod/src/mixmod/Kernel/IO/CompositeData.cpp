/***************************************************************************
                             SRC/mixmod/Kernel/IO/CompositeData.cpp  description
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
#include "mixmod/Kernel/IO/CompositeData.h"
#include "mixmod/Kernel/IO/CompositeSample.h"
#include "mixmod/Kernel/IO/Data.h"
#include "mixmod/Kernel/IO/BinaryData.h"
#include "mixmod/Kernel/IO/GaussianData.h"
#include "mixmod/Kernel/IO/QualitativeColumnDescription.h"
#include "mixmod/Kernel/IO/QuantitativeColumnDescription.h"
#include "mixmod/Kernel/IO/UnusedColumnDescription.h"
#include "mixmod/Kernel/IO/DataDescription.h"

namespace XEM {

CompositeData::CompositeData() {
	// TODO Auto-generated constructor stub
}

CompositeData::CompositeData(const CompositeData* cData) : Data(*cData) {
	_dataComponent.resize(2);
	_dataComponent[0] = (const_cast<CompositeData*> (cData)->getBinaryData())->clone();
	_dataComponent[1] = (const_cast<CompositeData*> (cData)->getGaussianData())->clone();
	_matrix = new Sample*[_nbSample];
	Sample ** bmatrix = _dataComponent[0]->_matrix;
	Sample ** gmatrix = _dataComponent[1]->_matrix;

	for (int i = 0; i < _nbSample; ++i) {
		_matrix[i] = new CompositeSample(bmatrix[i], gmatrix[i]);
	}
}

CompositeData::CompositeData(Data* bdata, Data* gdata)
: Data(bdata->_nbSample, bdata->_pbDimension + gdata->_pbDimension)
{
	if (typeid (*bdata) == typeid (gdata)) {
		THROW (InputException, badInputType);
	}
	if (typeid (*bdata) != typeid (BinaryData)) {
		//swap bdata and gdata
		Data * tmp = bdata;
		bdata = gdata;
		gdata = tmp;
	}

	assert(bdata->_nbSample == gdata->_nbSample);
	_dataComponent.resize(2);
	_dataComponent[0] = bdata;
	_dataComponent[1] = gdata;

  _matrix = new Sample*[_nbSample];
	Sample ** bmatrix = _dataComponent[0]->_matrix;
	Sample ** gmatrix = _dataComponent[1]->_matrix;

	for (int i = 0; i < _nbSample; ++i) {
		_matrix[i] = new CompositeSample(bmatrix[i], gmatrix[i]);
	}

}

CompositeData * CompositeData::clone() const {
	CompositeData * ctemp = new CompositeData(this);
	return ctemp;
}

CompositeData::~CompositeData() {
	for (unsigned int i = 0; i < _dataComponent.size(); ++i) {
		if (_dataComponent[i]) {
			delete _dataComponent[i];
			_dataComponent[i] = NULL;
		}
	}
  for (int64_t i = 0; i < _nbSample; i++) {
    delete _matrix[i];
  }
  delete[] _matrix;
}

void CompositeData::input(const DataDescription& dataDescription)
{
	std::vector<ColumnDescription*> gColumnDescription;
	std::vector<ColumnDescription*> bColumnDescription;

//for (std::vector<ColumnDescription*>::const_iterator it=dataDescription.getAllColumnDescription().begin(); it != dataDescription.getAllColumnDescription().end(); it++) {
  for (int64_t i = 0; i <dataDescription.getAllColumnDescription().size(); i++) {
//  ColumnDescription* columnDescription = *it;
    ColumnDescription* columnDescription =(dataDescription.getAllColumnDescription())[i];
		if (typeid (*columnDescription) == typeid (QualitativeColumnDescription)) {
			bColumnDescription.push_back(columnDescription);
			gColumnDescription.push_back(new UnusedColumnDescription(columnDescription->getIndex()));
		}
		else if (typeid (*columnDescription) == typeid (QuantitativeColumnDescription)) {
			gColumnDescription.push_back(columnDescription);
			bColumnDescription.push_back(new UnusedColumnDescription(columnDescription->getIndex()));
		}
		else {
			bColumnDescription.push_back(new UnusedColumnDescription(columnDescription->getIndex()));
			gColumnDescription.push_back(new UnusedColumnDescription(columnDescription->getIndex()));
		}
	}

	_dataComponent.resize(2);
	// bData
	DataDescription bDataDescription(
		dataDescription.getNbSample(), dataDescription.getNbColumn(), bColumnDescription,
		dataDescription.getFormat(), dataDescription.getFileName(), dataDescription.getInfoName());
	_dataComponent[0] = bDataDescription.getData();
	_dataComponent[0]->input(bDataDescription);
	bDataDescription.releaseData();
	//gdata
	DataDescription gDataDescription(
		dataDescription.getNbSample(), dataDescription.getNbColumn(), gColumnDescription,
		dataDescription.getFormat(), dataDescription.getFileName(), dataDescription.getInfoName());
	_dataComponent[1] = gDataDescription.getData();
	_dataComponent[1]->input(gDataDescription);
	gDataDescription.releaseData();

	_matrix = new Sample*[_nbSample];
	Sample** bmatrix = _dataComponent[0]->_matrix;
	Sample** gmatrix = _dataComponent[1]->_matrix;

	for (int i = 0; i < _nbSample; ++i)
		_matrix[i] = new CompositeSample(bmatrix[i], gmatrix[i]);
}

}
