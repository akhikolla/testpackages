/***************************************************************************
                             SRC/mixmod/Kernel/IO/ProbaDescription.cpp  description
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

#include "mixmod/Kernel/IO/ProbaDescription.h"
#include "mixmod/Kernel/IO/Proba.h"
#include "mixmod/Kernel/Model/Model.h"
#include "mixmod/Kernel/IO/QuantitativeColumnDescription.h"
#include <sstream>

namespace XEM {

//------------
// Constructor by default
//------------
ProbaDescription::ProbaDescription() : Description() {
	_proba = NULL;
}

// ---------------------------
//constructor by initilization
// ---------------------------
ProbaDescription::ProbaDescription(int64_t nbSample, int64_t nbCluster, 
		FormatNumeric::FormatNumericFile format, std::string filename, std::string infoName) 
{
	_infoName = "infoName";
	_nbSample = nbSample;
	_nbColumn = nbCluster;
	_fileName = filename;
	_format = format;
	_columnDescription.resize(nbCluster);
	for (int64_t i = 0; i < nbCluster; i++) {
		_columnDescription[i] = new QuantitativeColumnDescription(i);
		std::string name("Proba cluster=");
		std::ostringstream sNum;
		sNum << (i + 1);
		name.append(sNum.str());
		_columnDescription[i]->setName(name);
	}
	_proba = new Proba(_nbSample, nbCluster);
	std::ifstream fi(filename.c_str(), ios::in);
	if (!fi.is_open()) {
		THROW(InputException, wrongLabelFileName);
	}
	_proba->input(fi);
}

//------------
// Constructor after an estimation->run
//------------
ProbaDescription::ProbaDescription(Model * model) : Description() {
	if (model) {
		_infoName = "Probability";
		_nbSample = model->getNbSample();
		_nbColumn = model->getNbCluster();
		_fileName = "";
		_format = FormatNumeric::txt;
		_columnDescription.resize(_nbColumn);
		for (int64_t iCol = 0; iCol < _nbColumn; iCol++) {
			_columnDescription[iCol] = new QuantitativeColumnDescription(iCol);
			std::string name("Probability for cluster ");
			std::ostringstream sNum;
			sNum << (iCol + 1);
			name.append(sNum.str());
			_columnDescription[iCol]->setName(name);
		}
		_proba = new Proba(model);
	}
	else {
		THROW(OtherException, nullPointerError);
	}
}

//------------
// Constructor by copy
//------------
ProbaDescription::ProbaDescription(ProbaDescription & probaDescription) {
	(*this) = probaDescription;
}

//------------
// operator ==
//------------
bool ProbaDescription::operator==(ProbaDescription & probaDescription) const {
	if (_fileName != probaDescription._fileName) return false;
	if (_format != probaDescription._format) return false;
	if (_infoName != probaDescription._infoName) return false;
	if (_nbSample != probaDescription._nbSample) return false;
	if (_nbColumn != probaDescription._nbColumn) return false;
	for (int64_t i = 0; i < _nbColumn; ++i) {
		if (_columnDescription[i]->getName() 
				!= probaDescription.getColumnDescription(i)->getName())
		{
			return false;
		}
	}
	if (!(_proba == probaDescription.getProba())) return false;
	return true;
}

//------------
// operator =
//------------
ProbaDescription & ProbaDescription::operator=(ProbaDescription & probaDescription) {
	_fileName = probaDescription._fileName;
	_format = probaDescription._format;
	_infoName = probaDescription._infoName;
	_nbSample = probaDescription._nbSample;
	_nbColumn = probaDescription._nbColumn;
	_columnDescription.resize(_nbColumn);
	for (int64_t i = 0; i < _nbColumn; ++i) {
		const ColumnDescription * cd = probaDescription.getColumnDescription(i);
		_columnDescription[i] = cd->clone();
	}
	_proba = new Proba(*(probaDescription.getProba()));
	return *this;
}

//------------
// Destructor
//------------
ProbaDescription::~ProbaDescription() {
	if (_proba) {
		delete _proba;
	}
}

//--------
// ostream
//--------
void ProbaDescription::saveNumericValues(std::string fileName) {
	//if (_fileName==""){
	std::ofstream fo(fileName.c_str(), ios::out);
	_proba->edit(fo);
	_fileName = fileName;
	//}
	/* else : if _fileName!="", probaDescription has been created by a XML file.
	In this case, the numeric file already exists. 
	 */
}

}
