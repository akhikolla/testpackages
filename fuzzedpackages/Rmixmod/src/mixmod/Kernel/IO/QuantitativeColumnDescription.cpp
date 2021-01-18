/***************************************************************************
                             SRC/mixmod/Kernel/IO/QuantitativeColumnDescription.cpp  description
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

#include "mixmod/Kernel/IO/QuantitativeColumnDescription.h"

namespace XEM {

QuantitativeColumnDescription::QuantitativeColumnDescription() : ColumnDescription() {
}

QuantitativeColumnDescription::QuantitativeColumnDescription(int64_t index) 
: ColumnDescription(index) {
}

QuantitativeColumnDescription::~QuantitativeColumnDescription() {
}

std::string QuantitativeColumnDescription::editType() {
	return "Quantitative";
}

ColumnDescription * QuantitativeColumnDescription::clone()const {
	QuantitativeColumnDescription * QCD = new QuantitativeColumnDescription();
	QCD->_index = _index;
	QCD->_name = _name;
	return QCD;
}

}
