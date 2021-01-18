/***************************************************************************
                             SRC/mixmod/Kernel/IO/QuantitativeColumnDescription.h  description
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

#ifndef XEMQUANTITATIVECOLUMNDESCRIPTION_H
#define XEMQUANTITATIVECOLUMNDESCRIPTION_H

#include "mixmod/Kernel/IO/ColumnDescription.h"

namespace XEM {

/** 
 \class XEMQuantitativeColumnDescription
 @author F. Langrognet
		@date 2011
		@brief XEMQuantitativeColumnDescription class derived from XEMColumnDescription
 */
class QuantitativeColumnDescription : public ColumnDescription {

public:
	
	/// Default constructor
	QuantitativeColumnDescription();

	/// initialization constructor
	QuantitativeColumnDescription(int64_t index);

	/// Destructor
	virtual ~QuantitativeColumnDescription();

	std::string editType();

	ColumnDescription * clone()const;
};

}

#endif // XEMQUANTITATIVECOLUMNDESCRIPTION_H
