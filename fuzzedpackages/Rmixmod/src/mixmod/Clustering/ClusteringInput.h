/***************************************************************************
                             SRC/mixmod/Clustering/ClusteringInput.h  description
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
#ifndef XEMCLUSTERINGINPUT_H
#define XEMCLUSTERINGINPUT_H

#include "mixmod/Kernel/IO/Input.h"

namespace XEM {

// pre-declaration
class ClusteringStrategy;

/**
 \class XEMClusteringInput
 @author F. Langrognet
		@date 2012
		@brief XEMClusteringInput derived from XEMInput
 */

class ClusteringInput : public Input {

public:

	/// Default Constructor
	ClusteringInput();

	/// Copy Constructor
	ClusteringInput(const ClusteringInput & CInput);

	/// Initialisation constructor
	ClusteringInput(const std::vector<int64_t> & iNbCluster,
			const DataDescription & iDataDescription);

	/// Destructor
	virtual ~ClusteringInput();

	// getStrategy
	ClusteringStrategy * getStrategy() const;

	// setStrategy
	void setStrategy(ClusteringStrategy * strat);

	/// setCriterionName
	virtual void setCriterion(std::vector<CriterionName> const & criterionName);

	/// setCriterionName
	virtual void setCriterion(const CriterionName criterionName, unsigned int index);

	///insertCriterionName[i]
	virtual void insertCriterion(const CriterionName criterionName, unsigned int index);

	// add a new criterion
	void addCriterion(const CriterionName criterionName);

	// print input
	virtual void edit(std::ostream & out ) const;

	/// model add, insert, set
	//By default, HD models are not allowed


	/// setModelType
	virtual void setModelType(const ModelType * modelType, unsigned int index);

	/// insertModelType
	virtual void insertModelType(const ModelType * modelType, unsigned int index);

	/// add new model type (at the end)
	virtual void addModelType(const ModelType * modelType);

	/// add new model (modelName -> modelType)
	virtual void addModel(ModelName const modelName);

	/// setModel (modelName -> modelType)
	virtual void setModel(std::vector<ModelName> const & modelName);

protected:

	/// verif
	virtual bool verif();
	// Clustering strategy
	ClusteringStrategy * _strategy;
};

// getStrategy
inline ClusteringStrategy * ClusteringInput::getStrategy() const {
	return _strategy;
}

}

#endif
