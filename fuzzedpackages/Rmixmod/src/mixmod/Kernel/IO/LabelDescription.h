/***************************************************************************
                             SRC/mixmod/Kernel/IO/LabelDescription.h  description
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

#ifndef XEMLABELDESCRIPTION_H
#define XEMLABELDESCRIPTION_H

#include "mixmod/Kernel/IO/Description.h"

namespace XEM {

// pre-declaration
class Model;
class Label;

/** 
 \class XEMLabelDescription
 @author F. Langrognet
		@date 2011
		@brief XEMLabelDescription class derived from XEMDescription
 */
class LabelDescription : public Description {

public:
	/// Default constructor
	LabelDescription();

	///constructor by initilization
	LabelDescription(int64_t nbSample, int64_t nbColumn, 
			std::vector< ColumnDescription* > columnDescription, 
			FormatNumeric::FormatNumericFile format, 
			std::string filename, std::string infoName = "");

	/// constructor from a vector of int
	LabelDescription(int64_t nbSample, std::vector<int64_t> vLabel);

	///constructor after an estimation->run
	LabelDescription(Model * estimation);

	///constructor by copy
	LabelDescription(LabelDescription & labelDescription);

	/// Destructor
	virtual ~LabelDescription();

	/// Comparison operator
	bool operator ==(const LabelDescription & labelDescription) const;

	///operator=    
	LabelDescription & operator=(LabelDescription & labelDescription);

	const Label * getLabel() const;
	Label * getLabel();

	const int64_t getNbCluster() const;

	void saveNumericValues(std::string fileName = "");

private:

	Label * _label;
	Label * createLabel();
	int64_t _nbCluster;
};

inline const Label * LabelDescription::getLabel() const {
	return _label;
}

inline Label * LabelDescription::getLabel() {
	return _label;
}

inline const int64_t LabelDescription::getNbCluster() const {
	return _nbCluster;
}

}

#endif // XEMLABELDESCRIPTION_H
