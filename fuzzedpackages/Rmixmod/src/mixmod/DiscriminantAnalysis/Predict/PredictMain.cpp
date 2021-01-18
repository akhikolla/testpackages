/***************************************************************************
                             SRC/mixmod/DiscriminantAnalysis/Predict/PredictMain.cpp  description
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

#include "mixmod/DiscriminantAnalysis/Predict/PredictMain.h"
#include "mixmod/DiscriminantAnalysis/Predict/PredictInput.h"
#include "mixmod/DiscriminantAnalysis/Predict/PredictOutput.h"
#include "mixmod/DiscriminantAnalysis/Predict/PredictStrategy.h"
#include "mixmod/Kernel/IO/BinaryData.h"
#include "mixmod/Kernel/IO/LabelDescription.h"
#include "mixmod/Kernel/Model/Model.h"
#include "mixmod/Kernel/Model/BinaryModel.h"
#include "mixmod/Utilities/Error.h"
#include "mixmod/Kernel/IO/Partition.h"
#include <ctime>

namespace XEM {

//------------
// Constructor
//------------
PredictMain::PredictMain() {
	THROW(OtherException, internalMixmodError);
}

//------------
// Constructor
//------------
PredictMain::PredictMain(PredictInput * input, PredictOutput * output)
: _input(input), _output(output) {
}

//-----------
// Destructor
//-----------
PredictMain::~PredictMain() {
	if (_output) {
		delete _output;
	}
}

//---
//run
//---
void PredictMain::run(IoMode iomode, int verbose, int massiccc) {

//TODO: naming conventions ? (this doesn't look good)
	IOMODE = iomode;
	VERBOSE = verbose;
  MASSICCC = massiccc;

	if (!_input) {
		THROW(OtherException, nullPointerError);
	}
	if (!_input->isFinalized()) {
		THROW(InputException, inputNotFinalized);
	}

	ModelType * modelType = _input->getModelType()[0];
	int64_t nbCluster = _input->getNbCluster(0);
	Data * data = (_input->getDataDescription()).getData();

	// define a new estimation
	Model * estimation;

	// create model for binary data
	if (_input->getDataType() == QualitativeData && DATA_REDUCE) {

		Data * workingData = data;
		Partition * inputKnownPartition = NULL;
		Partition * workingKnownPartition = NULL;

		std::vector<int64_t> correspondenceOriginDataToReduceData;

		//--------------------------------
		//Reduce Data
		//------------
		BinaryData * bData = dynamic_cast<BinaryData*> (data);

		// initPartition
		Partition * inputInitPartition = NULL;
		Partition * workingInitPartition = NULL;

		try {
			//TODO RD : data ne doit pas forcément etre recréé
			workingData = bData->reduceData(
					correspondenceOriginDataToReduceData, inputKnownPartition, 
					inputInitPartition, workingKnownPartition, workingInitPartition);
			/* TODO ?:
			 if inputKnownPartition : delete workingKnownPartition 
			 if inputInitPartition : delete workingStrategy, workingInitPartition
       
			 */
		}
		catch (Exception& errorType) {
			workingData = NULL;
			throw;
		}
		// fin de ReduceData

		// create new estimation
		estimation = new BinaryModel(modelType, nbCluster, workingData, 
				workingKnownPartition, correspondenceOriginDataToReduceData);
	}
		// create model for quantitative data
	else if (_input->getDataType() == QuantitativeData) {
		// create new estimation
		estimation = new Model(modelType, nbCluster, data, NULL);
	}
		// create model for heterogeneous data
	else {
		// create new estimation
		estimation = new Model(modelType, nbCluster, data, NULL);
	}

	// create new strategy
	PredictStrategy strategy(_input->getClassificationRule());

	try {
		strategy.run(estimation);
	}
	catch (Exception&errorType) {
		
		if (VERBOSE == 1) {
			Error error(errorType);
			error.run();
		}
		
		// set error for that model
		estimation->setError(errorType);
	}
	catch (...) {
		Exception * unknown = new OtherException(UnknownReason);
		estimation->setError(*unknown);
		delete unknown;
	}
	// create output
	_output = new PredictOutput(estimation);

// release memory
	delete estimation;
}

//------------------------
// return pointer to Input
//------------------------
Input * PredictMain::getInput() {
	if (_input) {
		return _input;
	}
	else {
		THROW(OtherException, nullPointerError);
	}
}

}
