/***************************************************************************
                             SRC/mixmod/Kernel/IO/BinarySample.h  description
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
#ifndef XEMBINARYSample_H
#define XEMBINARYSample_H

#include "mixmod/Kernel/IO/Sample.h"

namespace XEM {

/**
  @brief Base class for Sample
  @author F Langrognet 
 */

class BinarySample : public Sample {

public:

	/// Constructor
	BinarySample();

	/// Constructor
	BinarySample(int64_t pbDimension);

	/// Constructor
	BinarySample(BinarySample * iSample);

	/// Constructor
	BinarySample(int64_t pbDimension, int64_t * tabValue);

	/// Destructor
	virtual ~BinarySample();

	/// Set value vector of sample
	void setDataTabValue(int64_t * tabValue);

	/// Set one value of sample
	void setDataValue(int64_t idxDim, int64_t iValue);

	/// get value vector of sample
	int64_t * getTabValue() const;

	/// get one value of sample
	int64_t getDataValue(int64_t idxDim) const;


protected:

	/// Vector of sample value
	int64_t * _value;
};

//---------------
// inline methods
//---------------

inline int64_t * BinarySample::getTabValue() const {
	return _value;
}

inline int64_t BinarySample::getDataValue(int64_t idxDim) const {
	return _value[idxDim];
}

}

#endif
