/***************************************************************************
                             SRC/mixmod/Kernel/IO/LabelDescription.cpp  description
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

#include "mixmod/Kernel/IO/LabelDescription.h"
#include "mixmod/Kernel/IO/Label.h"
#include "mixmod/Kernel/Model/Model.h"
#include "mixmod/Kernel/IO/WeightColumnDescription.h"
#include "mixmod/Kernel/IO/QualitativeColumnDescription.h"
#include "mixmod/Kernel/IO/QuantitativeColumnDescription.h"
#include "mixmod/Kernel/IO/IndividualColumnDescription.h"
#include <algorithm>

namespace XEM {

//------------
// Constructor by default
//------------
LabelDescription::LabelDescription() : Description() {
	_label = NULL;
	_nbCluster = 0;
}

// ---------------------------
//constructor by initilization
// ---------------------------
LabelDescription::LabelDescription(
		int64_t nbSample,
		int64_t nbColumn,
		std::vector< ColumnDescription* > columnDescription,
		FormatNumeric::FormatNumericFile format,
		std::string filename,
		std::string infoName)
: Description(nbSample, nbColumn, columnDescription, format, filename, infoName)
{
	_label = createLabel();
	// get the number of cluster
	_nbCluster = *max_element(_label->getLabel().begin(), _label->getLabel().end());
}

//---------------------------------
// constructor from a vector of int
//----------------------------------
LabelDescription::LabelDescription(int64_t nbSample, std::vector<int64_t> vLabel)
: Description()
{
	// get the number of cluster
	_nbCluster = *max_element(vLabel.begin(), vLabel.end());
	_infoName = "Label";
	_nbSample = nbSample;
	_nbColumn = 1;
	_fileName = "";
	_format = FormatNumeric::txt;
	_columnDescription.resize(1);
	_columnDescription[0] = new QualitativeColumnDescription(0, _nbCluster);
	std::string name("Label");
	_columnDescription[0]->setName(name);

  if (_nbSample != vLabel.size())
    THROW (InputException, badNumberOfValuesInLabelInput);

	_label = new Label(_nbSample);
	_label->setLabel(vLabel, _nbSample);
}

//-------------------------------------
// Constructor after an estimation->run
//--------------------------------------
LabelDescription::LabelDescription(Model * estimation) : Description() {
	if (estimation) {
		_infoName = "Label";
		_nbSample = estimation->getNbSample();
		_nbColumn = 1;
		_fileName = "";
		_format = FormatNumeric::txt;
		_columnDescription.resize(1);
		_columnDescription[0] = new QualitativeColumnDescription(0, estimation->getNbCluster());
		std::string name("Label");
		_columnDescription[0]->setName(name);
		_label = new Label(estimation);
		_nbCluster = estimation->getNbCluster();
	}
	else {
		THROW(OtherException, nullPointerError);
	}
}

//------------
// Constructor by copy
//------------
LabelDescription::LabelDescription(LabelDescription & labelDescription) {
	(*this) = labelDescription;
}

//------------
// operator =
//------------
LabelDescription & LabelDescription::operator=(LabelDescription & labelDescription) {
	_fileName = labelDescription._fileName;
	_format = labelDescription._format;
	_infoName = labelDescription._infoName;
	_nbSample = labelDescription._nbSample;
	_nbColumn = labelDescription._nbColumn;
	_columnDescription.resize(_nbColumn);
	_nbCluster = labelDescription._nbCluster;
	_label = new Label(*(labelDescription.getLabel()));
	return *this;
}

//---------------------
/// Comparison operator
//---------------------
bool LabelDescription::operator ==(const LabelDescription & labelDescription) const {
	if (_fileName != labelDescription._fileName) return false;
	if (_format != labelDescription._format) return false;
	if (_infoName != labelDescription._infoName) return false;
	if (_nbSample != labelDescription._nbSample) return false;
	if (_nbColumn != labelDescription._nbColumn) return false;
	for (int64_t i = 0; i < _nbColumn; ++i) {
		if (_columnDescription[i]->getName()
				!= labelDescription.getColumnDescription(i)->getName())
		{
			return false;
		}
	}
	if (_nbCluster != labelDescription._nbCluster) return false;
	if (!(_label == labelDescription.getLabel())) return false;
	return true;
}

//------------
// Destructor
//------------
LabelDescription::~LabelDescription() {
	if (_label) delete _label;
}

//--------
// ostream
//--------
void LabelDescription::saveNumericValues(std::string fileName) {
	//if (_fileName==""){
	std::ofstream fo(fileName.c_str(), ios::out);
	_label->edit(fo);
	_fileName = fileName;
	//}
	/* else : if _fileName!="", labelDescription has been created by a XML file.
	In this case, the numeric file already exists.
	 */
}

Label* LabelDescription::createLabel() {
	Label * label = new Label();

	int64_t nbQualitativeVariable = 0;
	int64_t nbIndividualVariable = 0;

	for (int64_t i = 0; i < _nbColumn; i++) {
		if (typeid (*(_columnDescription[i])) == typeid (QualitativeColumnDescription)) {
			nbQualitativeVariable++;
		}
		if (typeid (*(_columnDescription[i])) == typeid (QuantitativeColumnDescription)) {
			THROW(InputException, badLabelDescription);
		}
		if (typeid (*(_columnDescription[i])) == typeid (WeightColumnDescription)) {
			THROW(InputException, tooManyWeightColumnDescription);
		}
		if (typeid (*(_columnDescription[i])) == typeid (IndividualColumnDescription)) {
			nbIndividualVariable++;
		}
	}

	if (nbQualitativeVariable != 1 || nbIndividualVariable > 1) {
		THROW(InputException, badLabelDescription);
	}
	label-> input(*this);
	return label;
}

}
