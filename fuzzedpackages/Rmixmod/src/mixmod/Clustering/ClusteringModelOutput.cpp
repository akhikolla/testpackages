/***************************************************************************
                             SRC/mixmod/Clustering/ClusteringModelOutput.cpp  description
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
#include "mixmod/Clustering/ClusteringModelOutput.h"

namespace XEM {

//--------------------
// Default Constructor
//--------------------
ClusteringModelOutput::ClusteringModelOutput() {
}

//-----------------
//  Copy constructor
//-----------------
ClusteringModelOutput::ClusteringModelOutput(const ClusteringModelOutput & cModelOutput) {
	THROW(OtherException, internalMixmodError);
}

//-----------------
//  Initialization Constructor
//-----------------
ClusteringModelOutput::ClusteringModelOutput(Model * estimation) : ModelOutput(estimation) {
}

//-----------------
//  Initialization Constructor
//-----------------
ClusteringModelOutput::ClusteringModelOutput( ModelType & modelType,
		int64_t nbCluster,
		std::vector< CriterionOutput* >& criterionOutput,
		double likelihood,
		ParameterDescription& parameterDescription,
		LabelDescription& labelDescription,
		ProbaDescription& probaDescription)
: ModelOutput(modelType, nbCluster, criterionOutput, likelihood, parameterDescription, labelDescription, probaDescription) {
}

//-----------------
//  Initialization Constructor
//-----------------
ClusteringModelOutput::ClusteringModelOutput(ModelType& modelType, int64_t nbCluster, Exception& error)
: ModelOutput(modelType, nbCluster, error) {
}

//-----------
// Destructor
//-----------
ClusteringModelOutput::~ClusteringModelOutput() {
}

}
