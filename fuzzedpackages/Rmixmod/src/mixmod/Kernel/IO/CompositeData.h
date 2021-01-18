/***************************************************************************
                             SRC/mixmod/Kernel/IO/CompositeData.h  description
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
#ifndef XEMCOMPOSITEDATA_H_
#define XEMCOMPOSITEDATA_H_
/**@file XEMCompositeData.h
 * @brief Composite data class for heterogeneous clustering.
 * @author Parmeet Bhatia
 */
#include "mixmod/Kernel/IO/Data.h"
#include "mixmod/Kernel/IO/GaussianData.h"
#include "mixmod/Kernel/IO/BinaryData.h"

namespace XEM {

class CompositeData : public Data {

public:
	
	//default constructor
	CompositeData();
	//copy constructor
	CompositeData(const CompositeData*);
	//Initialization constructor
	CompositeData(Data*, Data*);
	virtual CompositeData * clone() const;

	virtual void input(std::ifstream&) {
		THROW(OtherException, FunctionNotYetImplemented);
	}

	virtual void input(const DataDescription&);

	virtual void output(std::ostream&) {
		THROW(OtherException, FunctionNotYetImplemented);
	}

	virtual Sample** cloneMatrix() {
		THROW(OtherException, FunctionNotYetImplemented);
	}

	/**type-cast overloading to return CompositeData::_dataComponent[1]*/
	inline operator GaussianData*() {
		return (GaussianData*) _dataComponent[1];
	}

	/**type-cast overloading to return CompositeData::_dataComponent[0]*/
	inline operator BinaryData*() {
		return (BinaryData*) _dataComponent[0];
	}

	/** get Gaussian data */
	virtual GaussianData* getGaussianData() {
		return (GaussianData*) _dataComponent[1];
	}

	/** get Gaussian data */
	virtual BinaryData* getBinaryData() {
		return (BinaryData*) _dataComponent[0];
	}
	/**Virtual Destructor*/
	virtual ~CompositeData();
	
protected:
	
	vector<Data*> _dataComponent;
};

}

#endif /* XEMCOMPOSITEDATA_H_ */
