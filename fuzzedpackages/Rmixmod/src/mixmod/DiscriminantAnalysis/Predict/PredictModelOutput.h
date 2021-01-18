/***************************************************************************
                             SRC/mixmod/DiscriminantAnalysis/Predict/PredictModelOutput.h  description
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
#ifndef XEMPREDICTMODELOUTPUT_H
#define XEMPREDICTMODELOUTPUT_H

#include "mixmod/Kernel/IO/ModelOutput.h"

namespace XEM {

// pre-declaration
class Model;

/** 
 \class XEMPredictModelOutput
 @author F. Langrognet - R Lebret
		@date 2012
		@brief XEMPredictModelOutput class derived from XEMModelOutput
 */
class PredictModelOutput : public ModelOutput {

public:

	/// Default Constructor
	PredictModelOutput();

	/// Initialization Constructor 1
	PredictModelOutput(Model * estimation);

	/// Initialization Constructor 2
	PredictModelOutput(ModelType & modelType, int64_t nbCluster, 
			std::vector<CriterionOutput*> & criterionOutput, double likelihood, 
			ParameterDescription & parameterDescription, LabelDescription & labelDescription,
			ProbaDescription & probaDescription);

	/// Initialization Constructor 3
	PredictModelOutput(ModelType & modelType, int64_t nbCluster, Exception& error);

	/// Copy Constructor
	PredictModelOutput(const PredictModelOutput & cModelOutput);

	/// Destructor
	virtual ~PredictModelOutput();
};

}

#endif
