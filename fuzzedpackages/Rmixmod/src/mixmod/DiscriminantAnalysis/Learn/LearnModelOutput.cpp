/***************************************************************************
                             SRC/mixmod/DiscriminantAnalysis/Learn/LearnModelOutput.cpp  description
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

#include "mixmod/DiscriminantAnalysis/Learn/LearnModelOutput.h"
#include "mixmod/Kernel/IO/LabelDescription.h"
#include "mixmod/Kernel/Model/BinaryModel.h"

namespace XEM {

//--------------------
// Default Constructor
//--------------------
LearnModelOutput::LearnModelOutput() {
    _CVLabel = nullptr;
}

//-----------------
//  Copy constructor
//-----------------
LearnModelOutput::LearnModelOutput(const LearnModelOutput & cModelOutput) {
	THROW(OtherException, internalMixmodError);
}

//-----------------
//  Initialization Constructor
//-----------------
LearnModelOutput::LearnModelOutput(Model * estimation) : ModelOutput(estimation) {
    _CVLabel = nullptr;
}

//-----------------
//  Initialization Constructor
//-----------------
LearnModelOutput::LearnModelOutput(ModelType & modelType, 
		int64_t nbCluster,
		std::vector< CriterionOutput* >& criterionOutput, 
		double likelihood, 
		ParameterDescription& parameterDescription,
		LabelDescription& labelDescription, 
		ProbaDescription& probaDescription)
: ModelOutput(modelType, nbCluster, criterionOutput, likelihood, 
		parameterDescription, labelDescription, probaDescription) {
  _CVLabel = nullptr;
}

//-----------------
//  Initialization Constructor
//-----------------
LearnModelOutput::LearnModelOutput(ModelType& modelType, int64_t nbCluster, Exception& error) 
: ModelOutput(modelType, nbCluster, error) {
    _CVLabel = nullptr;
}

//-----------
// Destructor
//-----------
LearnModelOutput::~LearnModelOutput() {
  if (_CVLabel) delete _CVLabel;
}

/// set CV Labels
void LearnModelOutput::setCVLabel(Model * estimation, std::vector<int64_t> & cvLabel) {
	if (isBinary(estimation->getModelType()->_nameModel) && DATA_REDUCE) {
		////
		//binary case

		// cmake a copy of cvLabel
		std::vector<int64_t> cvLabelCopy(cvLabel);

		const std::vector<int64_t> & correspondenceOriginDataToReduceData = 
				dynamic_cast<BinaryModel*> (estimation)->getCorrespondenceOriginDataToReduceData();
		// get the true number of sample
		const int64_t nbSample = correspondenceOriginDataToReduceData.size();

		// resize cvLabel
		cvLabel.resize(nbSample);

		// convert labelReduce, to label
		for (int64_t i = 0; i < nbSample; i++) {
			cvLabel[i] = cvLabelCopy[correspondenceOriginDataToReduceData[i]];
		}
	}
	// create a new instance of XEMLabelDescription
	_CVLabel = new LabelDescription(cvLabel.size(), cvLabel);
}

}
