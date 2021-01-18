/***************************************************************************
                             SRC/mixmod/Clustering/ClusteringModelOutput.h  description
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
#ifndef XEMCLUSTERINGMODELOUTPUT_H
#define XEMCLUSTERINGMODELOUTPUT_H

#include "mixmod/Kernel/IO/ModelOutput.h"

namespace XEM {

/** 
 \class XEMClusteringModelOutput
 @author F. Langrognet
		@date 2012
		@brief XEMClusteringModelOutput derived from XEMModelOutput
 */
class ClusteringModelOutput : public ModelOutput {

public:

	/// Default Constructor
	ClusteringModelOutput();

	/// Initialization Constructor 1
	ClusteringModelOutput(Model * estimation);

	/// Initialization Constructor 2
	ClusteringModelOutput(ModelType & modelType, int64_t nbCluster, 
			std::vector<CriterionOutput*> & criterionOutput, double likelihood, 
			ParameterDescription & parameterDescription, LabelDescription & labelDescription,  
			ProbaDescription & probaDescription);

	/// Initialization Constructor 3
	ClusteringModelOutput(ModelType & modelType, int64_t nbCluster, Exception& error);

	/// Copy Constructor
	ClusteringModelOutput(const ClusteringModelOutput & cModelOutput);

	/// Destructor
	virtual ~ClusteringModelOutput();
};

}

#endif
