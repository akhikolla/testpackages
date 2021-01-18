/***************************************************************************
                             SRC/mixmod/Kernel/IO/QualitativeColumnDescription.h  description
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

#ifndef XEMQUALITATIVECOLUMNDESCRIPTION_H
#define XEMQUALITATIVECOLUMNDESCRIPTION_H

#include "mixmod/Kernel/IO/ColumnDescription.h"

namespace XEM {

//Variable Description
struct VariableDescription {

	int64_t num;
	std::string name;
};

/** 
 \class XEMQualitativeColumnDescription
 @author F. Langrognet
		@date 2011
		@brief XEMQualitativeColumnDescription class derived from XEMColumnDescription
 */
class QualitativeColumnDescription : public ColumnDescription {

public:
	
	/// Default constructor
	QualitativeColumnDescription();

	/// initialization constructor
	QualitativeColumnDescription(int64_t index, int64_t nbFactor);

	/// Destructor
	virtual ~QualitativeColumnDescription();

	std::string editType();

	ColumnDescription * clone()const;

	const int64_t & getNbFactor()const;

	const std::vector<VariableDescription> & getVariableDescription()const;

	void setVariableDescription(VariableDescription & variableDescription, unsigned int index);

private:
	
	int64_t _nbFactor;
	std::vector<VariableDescription> _variableDescription;
};

inline const int64_t & QualitativeColumnDescription::getNbFactor()const {
	return _nbFactor;
}

inline const std::vector<VariableDescription> & QualitativeColumnDescription::getVariableDescription()const {
	return _variableDescription;
}

}

#endif // XEMQUALITATIVECOLUMNDESCRIPTION_H
