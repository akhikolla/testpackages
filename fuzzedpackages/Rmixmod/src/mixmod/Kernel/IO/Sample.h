/***************************************************************************
                             SRC/mixmod/Kernel/IO/Sample.h  description
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
#ifndef XEMSample_H
#define XEMSample_H

#include <stdint.h>
#include "mixmod/Utilities/Util.h"

namespace XEM {

class GaussianSample;
class BinarySample;

/**
  @brief Base class for Sample
  @author F Langrognet 
 */

class Sample {

public:

	/// Constructor
	Sample();

	/// Constructor
	Sample(Sample * iSample);

	/// Constructor
	Sample(int64_t pbDimension);

	/// Destructor
	virtual ~Sample();

	virtual GaussianSample* getGaussianSample() const {
		return (GaussianSample*)this;
	}

	virtual BinarySample* getBinarySample() const {
		return (BinarySample*)this;
	}
	
	/// Selector
	int64_t getPbDimension();

protected:

	/// Problem dimension
	int64_t _pbDimension;
};

inline int64_t Sample::getPbDimension() {
	return _pbDimension;
}

}

#endif
