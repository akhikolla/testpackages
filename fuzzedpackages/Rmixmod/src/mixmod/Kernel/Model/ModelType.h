/***************************************************************************
                             SRC/mixmod/Kernel/Model/ModelType.h  description
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
#ifndef XEMMODELTYPE_H
#define XEMMODELTYPE_H

#include "mixmod/Utilities/Util.h"

namespace XEM {

/**
  @brief Base class for ModelType(s)
  @author F Langrognet
 */

class ModelType {

public:

	/// Default constructor
	ModelType();

	// constructor
	ModelType(ModelName name, int64_t nbSubDimensionFree = 0);

	// copy constructor
	ModelType(const ModelType & iModelType);

	/// Destructor
	~ModelType();

	/// Comparison operator
	bool operator ==(const ModelType & modelType) const;

	/// Input model type
	void input(std::ifstream & fi, int64_t nbCluster);

	/// name of the model
	ModelName _nameModel;

	ModelType* clone();

	//// list of number of subDimensionEqual
	//int64_t _nbSubDimensionEqual;
	//// list of number of subDimensionFree
	//int64_t _nbSubDimensionFree;

	/// list of subDimensionEqual
	int64_t _subDimensionEqual;

	/// _nbSubDimensionFree : size of array _tabSubDimensionFree
	int64_t _nbSubDimensionFree;

	/// array of subDimensionFree
	int64_t * _tabSubDimensionFree;

	/// getModelName
	const ModelName & getModelName() const;

	/// getSubDimensionEqual
	const int64_t & getSubDimensionEqual() const;

	/// getTabSubDimensionFree
	const int64_t * getTabSubDimensionFree() const;

	///getTabSubDimensionFreeI
	const int64_t & getTabSubDimensionFreeI(int64_t index) const;

	/// setSubDimensionFree
	void setTabSubDimensionFree(int64_t iTabSubDimensionFree, int64_t position);

	/// setSubDimensionEqual
	void setSubDimensionEqual(int64_t iSubDimensionEqual);

	/// <<
	friend std::ostream & operator<<(std::ostream & fo, ModelType & modelType);

	// print out model type
	void print(std::ostream & flux) const;
	// print out model type short cut
	void printShortcut(std::ostream & flux) const;
	/// editModelType
	void edit(std::ostream & oFile);
};

inline const ModelName & ModelType::getModelName() const {
	return _nameModel;
}

inline const int64_t & ModelType::getSubDimensionEqual() const {
	return _subDimensionEqual;
}

inline const int64_t * ModelType::getTabSubDimensionFree() const {
	return _tabSubDimensionFree;
}

inline const int64_t & ModelType::getTabSubDimensionFreeI(int64_t index) const {
	return _tabSubDimensionFree[index];
}


}

#endif
