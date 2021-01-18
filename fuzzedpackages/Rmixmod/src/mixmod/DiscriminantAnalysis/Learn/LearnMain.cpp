/***************************************************************************
                             SRC/mixmod/DiscriminantAnalysis/Learn/LearnMain.cpp  description
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

#include "mixmod/DiscriminantAnalysis/Learn/LearnMain.h"
#include "mixmod/DiscriminantAnalysis/Learn/LearnInput.h"
#include "mixmod/DiscriminantAnalysis/Learn/LearnOutput.h"
#include "mixmod/DiscriminantAnalysis/Learn/LearnModelOutput.h"
#include "mixmod/DiscriminantAnalysis/Learn/LearnStrategy.h"
#include "mixmod/Utilities/Random.h"
#include "mixmod/Kernel/IO/BinaryData.h"
#include "mixmod/Kernel/IO/LabelDescription.h"
#include "mixmod/Kernel/Model/Model.h"
#include "mixmod/Kernel/Model/BinaryModel.h"
#include "mixmod/Utilities/Error.h"
#include "mixmod/Kernel/IO/Partition.h"
#include "mixmod/Kernel/Criterion/BICCriterion.h"
#include "mixmod/Kernel/Criterion/CVCriterion.h"
#include "mixmod/Kernel/IO/ModelOutput.h"
#include <ctime>

namespace XEM {

//------------
// Constructor
//------------
LearnMain::LearnMain() {
	THROW(OtherException, internalMixmodError);
}

//------------
// Constructor
//------------
LearnMain::LearnMain(LearnInput * input, LearnOutput * output)
: _input(input), _output(output) {
}

//-----------
// Destructor
//-----------
LearnMain::~LearnMain() {
	if (_output) {
		delete _output;
	}
}

//---
//run
//---
void LearnMain::run(int seed, IoMode iomode, int verbose, int massiccc) {
	//--------------------------------------------------------------------
	// Important : Call randomize function to change random seed or not ...
	//---------------------------------------------------------------------
	
//TODO: naming conventions ? (this doesn't look good)
	IOMODE = iomode;
	VERBOSE = verbose;
  MASSICCC = massiccc;

	initRandomize(seed);

	if (!_input) {
		THROW(OtherException, nullPointerError);
	}
	if (!_input->isFinalized()) {
		THROW(InputException, inputNotFinalized);
	}

	// get the number of models
	int64_t nbModelType = _input->getModelType().size();
	// get the list of number of clusters
	std::vector<int64_t> vNbCluster = _input->getNbCluster();
	// get the size of the cluster's list
	//  int64_t nbNbCluster = vNbCluster.size();
	// compute the total number of estimation to run
	int64_t nbEstimation = nbModelType;

	// get number of cluster
	int64_t nbCluster = vNbCluster[0];

	// get data
	Data * inputData = (_input->getDataDescription()).getData();

	// get labels
	Partition * inputKnownPartition = new Partition(
			_input->getKnownLabelDescription()->getLabel(), vNbCluster[0]);

	// define a new Learn Strategy
	LearnStrategy learnStrategy;

	// Define a vector of models
	std::vector<Model*> estimations(nbEstimation);

	// create model for binary data
	if (_input->getDataType() == QualitativeData && DATA_REDUCE) {

		Data * workingData = inputData;
		Partition * workingKnownPartition = inputKnownPartition;
		std::vector<int64_t> correspondenceOriginDataToReduceData;
		//------------
		//Reduce Data
		//------------
		BinaryData * bData = dynamic_cast<BinaryData*> (inputData);

		// initPartition
		Partition * inputInitPartition = NULL;
		Partition * workingInitPartition = NULL;

		try {
			//TODO RD : data ne doit pas forcément etre recréé
			workingData = bData->reduceData(correspondenceOriginDataToReduceData, 
					inputKnownPartition, inputInitPartition, 
					workingKnownPartition, workingInitPartition);
			/* TODO ?:
			 if inputKnownPartition : delete workingKnownPartition 
			 if inputInitPartition : delete workingStrategy, workingInitPartition
       
			 */
		}
		catch (Exception&errorType) {
			workingData = NULL;
			throw;
		}
		// fin de ReduceData
		
		// create tabEstmation[iEstimation]
		//--------------------------------
		for (int iModelType = 0; iModelType < nbModelType; iModelType++) {
			ModelType * modelType = _input->getModelType()[iModelType];
			estimations[iModelType] = new BinaryModel(modelType, nbCluster, workingData, 
					workingKnownPartition, correspondenceOriginDataToReduceData);
		}
	}
		// create model for quantitative data
	else if (_input->getDataType() == QuantitativeData) {
		// loop over all models
		for (int iModelType = 0; iModelType < nbModelType; iModelType++) {
			ModelType * modelType = _input->getModelType()[iModelType];
			// create new estimation
			estimations[iModelType] = 
					new Model(modelType, nbCluster, inputData, inputKnownPartition);
		}
	}
		//create model for heterogeneous data
	else {
		// loop over all models
		for (int iModelType = 0; iModelType < nbModelType; iModelType++) {
			ModelType * modelType = _input->getModelType()[iModelType];
			// create new estimation
			estimations[iModelType] = 
					new Model(modelType, nbCluster, inputData, inputKnownPartition);
		}
	}
	// release memory
	delete inputKnownPartition;

	//-------------------
	// 2. run Estimations
	//-------------------
	// loop over the number of estimation to do
	int64_t iEstimation = 0;
	while (iEstimation < nbEstimation) {
		try {
			learnStrategy.run(estimations[iEstimation]);
		}
		catch (Exception&errorType) {
			
			if (VERBOSE == 1) {
				Error error(errorType);
				error.run();
			}
			
			// set error for that model
			estimations[iEstimation]->setError(errorType);
		}
		catch (...) {
			Exception * unknown = new OtherException(UnknownReason);
			estimations[iEstimation]->setError(*unknown);
			delete unknown;
		}
		iEstimation++;
	}
	//----------------
	// create Output
	//-----------------
	_output = new LearnOutput(estimations);
	//-------------------
	// 3. run Selections
	//-------------------
	// get criterion names
	std::vector<CriterionName> const & criterion = _input->getCriterionName();

	// Initialize timer info	
	time_t startTime;
	ofstream progressFile;
	if (MASSICCC == 1) {
		time(&startTime);
	}

	// loop over all models
	for (unsigned int iModel = 0; iModel < nbEstimation; iModel++) {
		// check whether error occurred for that model
		if ((estimations[iModel]->getErrorType()) == NOERROR) {
			// loop over criterion name
			for (unsigned int iCriterion = 0; iCriterion < criterion.size(); iCriterion++) {
				switch (criterion[iCriterion]) {
				case BIC:
				{
					// create BIC criterion
					BICCriterion bic(estimations[iModel]);
					// compute criterion outputs
					bic.run(_output->getLearnModelOutput(iModel)->getCriterionOutput(BIC));
					break;
				}
				case CV:
				{
					// create CV criterion

					CVCriterion cv(estimations[iModel], _input->getNbCVBlock());
					// compute criterion outputs
					cv.run(_output->getLearnModelOutput(iModel)->getCriterionOutput(CV));
					// add label from CV
					_output->getLearnModelOutput(iModel)
							->setCVLabel(estimations[iModel], cv.getCVLabel());
					break;
				}
				case ICL: THROW(InputException, badCriterion);
				case NEC: THROW(InputException, badCriterion);
				case UNKNOWN_CRITERION_NAME: THROW(OtherException, internalMixmodError);
				default: THROW(OtherException, internalMixmodError);
				}

        //Write progress in file
        if (MASSICCC == 1) {
          progressFile.open ("progress.json");
          progressFile << "{ \"Progress\" :  " << ((double)iCriterion + 1 + (double)iModel * (double)criterion.size())/((double)criterion.size() * (double)nbEstimation) * 100.0;
          time_t currTime;
          time(&currTime);
          double timePerModel = difftime(currTime, startTime) / ((double)iCriterion + 1 + (double)iModel * (double)criterion.size());
          progressFile << ", \"Estimated remaining time\" : " << timePerModel * (criterion.size() - iCriterion - 1) +  (nbEstimation - iModel - 1) * criterion.size() * timePerModel << " } ";
          progressFile.close();
        }

      } // end iCriterion
    }
		else { // set criterion error
			// loop over criterion name
			for (unsigned int iCriterion = 0; iCriterion < criterion.size(); iCriterion++) {
				_output->getLearnModelOutput(iModel)->setCriterionOutput(CriterionOutput(
						criterion[iCriterion], 0.0, estimations[iModel]->getErrorType()));
			} // end iCriterion
		}
    //Compute Entropy if needed (for massiccc's visualization)
    if (MASSICCC == 1 || MASSICCC == 10 || MASSICCC == 11) {
      vector< vector<double> > entropyM = estimations[iModel]->getEntropyMatrix();
      int nbVariables = estimations[iModel]->getData()->getPbDimension();
      ofstream entropyFile;
      entropyFile.open ("entropy" + std::to_string(iModel + 1) + ".txt");
      for (int j = 0; j < nbVariables; j++) {
        for (int l = 0; l < nbCluster; l++) {
          entropyFile << entropyM[j][l] << "  ";
        }
        entropyFile << endl;
      } 
      entropyFile.close();
    }

  }//end iModel

	// release memory
	// NOTE [bauder, 2013-02-24]: not everything should be released here, 
	// since line 185 in examples/.../discriminant_analysis_3.cpp
	// ( "XEMParameterDescription * paramPredict = new XEMParameterDescription ( lMOutput );" )
	// attempts to use some feature of these estimations.
	// However, fixing memory management here seems rather subtle.
	// To be fixed soon [see also mixmodLib#5279]
	//   for (unsigned int iModel=0; iModel<nbEstimation; iModel++){
	//     delete estimations[iModel];
	//     estimations[iModel] = 0;
	//   }
}

//------------------------
// return pointer to Input
//------------------------
Input * LearnMain::getInput() {
	if (_input) {
		return _input;
	}
	else {
		THROW(OtherException, nullPointerError);
	}
}

}
