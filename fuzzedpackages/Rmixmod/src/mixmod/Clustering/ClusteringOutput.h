/***************************************************************************
                             SRC/mixmod/Clustering/ClusteringOutput.h  description
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
#ifndef XEMCLUSTERINGOUTPUT_H
#define XEMCLUSTERINGOUTPUT_H

#include "mixmod/Utilities/Util.h"

namespace XEM {

// pre-declaration
class Model;
class Criterion;
class ClusteringModelOutput;

/** 
 \class XEMClusteringOutput
 Main class for Clustering Output
 @author F. Langrognet
		@date 2012
		@brief XEMClusteringOutput class
 */
class ClusteringOutput {

public:

	/// Default Constructor
	ClusteringOutput(std::vector<CriterionName> const & criterionName);

	/// Copy Constructor
	ClusteringOutput(const ClusteringOutput & cOutput);

	ClusteringOutput* clone();

	// Add one ClusteringModelOutput (incremental construction). Count is needed because vector.puch_back won't work with openMP
	void addEstimation(ClusteringModelOutput* cmoutput, int64_t count);

	// Resize ClusteringModelOutput vector
	void clusteringModelOutputResize(int64_t size);

	/// Initialization constructor
	ClusteringOutput(std::vector<Model*> const & estimations, std::vector<CriterionName> const & criterionName);

	/// Destructor
	virtual ~ClusteringOutput();

	/// Comparison operator
	bool operator ==(const ClusteringOutput & output) const;

	bool operator !=(const ClusteringOutput & output) const;

	bool atLeastOneEstimationNoError() const;

	const int getNbEstimationWithNoError() const;

	/// sort vector of XEMClusteringModelOutput (with the ith criterion value)
	void sort(CriterionName criterionName);

	void editFile() const;

	/// return the index'th' ClusteringModelOutput
	/// Note : index is between 0 and size(ClusteringModelOutput)-1
	ClusteringModelOutput * getClusteringModelOutput(const int64_t index) const;

	int64_t getNbClusteringModelOutput() const;

	std::vector<ClusteringModelOutput*> const &  getClusteringModelOutput() const;

	void setClusteringModelOutput(std::vector<ClusteringModelOutput *> & clusteringModelOutput);

	const int getCriterionSize() const;
	const CriterionName & getCriterionName(const int index) const;
	const std::vector<CriterionName> & getCriterionName() const;

private:
	
	// Vector containing output for each model
	std::vector<ClusteringModelOutput*> _clusteringModelOutput;
	// vector containing criterion name
	// that will be useful to deal with output in mixmodGUI
	std::vector<CriterionName> const & _criterionName;
};

inline std::vector<ClusteringModelOutput*> const & ClusteringOutput::getClusteringModelOutput() const {
	return _clusteringModelOutput;
}

inline ClusteringModelOutput * ClusteringOutput::getClusteringModelOutput(const int64_t index) const {
	return _clusteringModelOutput[index];
}

inline int64_t ClusteringOutput::getNbClusteringModelOutput() const {
	return _clusteringModelOutput.size();
}

inline const int ClusteringOutput::getCriterionSize() const {
	return _criterionName.size();
}

inline const CriterionName & ClusteringOutput::getCriterionName(const int index) const {
	return _criterionName[index];
}

inline const std::vector<CriterionName> & ClusteringOutput::getCriterionName() const {
	return _criterionName;
}

}

#endif
