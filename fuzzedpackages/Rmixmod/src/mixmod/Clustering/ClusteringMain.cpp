/***************************************************************************
                             SRC/mixmod/Clustering/ClusteringMain.cpp  description
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

#include "mixmod/Clustering/ClusteringMain.h"
#include "mixmod/Clustering/ClusteringInput.h"
#include "mixmod/Clustering/ClusteringOutput.h"
#include "mixmod/Clustering/ClusteringModelOutput.h"
#include "mixmod/Clustering/ClusteringStrategy.h"
#include "mixmod/Clustering/ClusteringStrategyInit.h"
#include "mixmod/Utilities/Random.h"
#include "mixmod/Kernel/IO/BinaryData.h"
#include "mixmod/Kernel/IO/Label.h"
#include "mixmod/Kernel/Model/BinaryModel.h"
#include "mixmod/Kernel/IO/LabelDescription.h"
#include "mixmod/Kernel/IO/Partition.h"
#include "mixmod/Utilities/Error.h"
#include "mixmod/Kernel/Criterion/BICCriterion.h"
#include "mixmod/Kernel/Criterion/ICLCriterion.h"
#include "mixmod/Kernel/Criterion/NECCriterion.h"
#include "mixmod/Kernel/IO/ModelOutput.h"
#include <ctime>

namespace XEM {

ClusteringMain::ClusteringMain() {
	THROW(OtherException, internalMixmodError);
}

ClusteringMain::ClusteringMain(ClusteringInput * input, ClusteringOutput * output)
: _input(input), _output(output) {
}

ClusteringMain::~ClusteringMain() {
	// Deleting NULL or nullptr is fine.
	// No check, because two calls of the destructor on the same object is a bug.
  delete _output;
}

void ClusteringMain::run(int seed, IoMode iomode, int verbose, int massiccc) {
	//TODO: naming conventions ? (this doesn't look good)
	IOMODE = iomode;
	VERBOSE = verbose;
  MASSICCC = massiccc;

  bool atLeastOneModelOk = false;

	// Call randomize() function to use a random seed, or antiRandomize() for deterministic seed.
	initRandomize(seed);

	if (!_input)
		THROW(OtherException, nullPointerError);
	// TODO: finalization should occur (only) here.
	if (! _input->isFinalized())
		THROW(InputException, inputNotFinalized);

	//----------------------------------------------
	// 0. Data preparation, variables initialization
	//----------------------------------------------

	// Check if the user provided a partition (hard or fuzzy)
	// TODO [bauder]: fuzzy partitions should be allowed, and used as t_ik initial values.
	//                In this case, we would start with an M-step and then iterate as usual.
	Partition* workingKnownPartition = _input->getKnownPartition();
	if (_input->getKnownLabelDescription()) {
		const Label * lab = _input->getKnownLabelDescription()->getLabel();
		workingKnownPartition = new Partition(lab, _input->getNbCluster()[0]);
		if (_input->getKnownPartition()) {
			// knowPartition and knownLabelDescription can't be both not NULL !
			THROW(OtherException, internalMixmodError);
		}
	}

  // Initialize timer info	
  int64_t count = 0; // = nbModel_i + 1 + nbCluster_i * nbModel. Used for timer info with openMP.
  time_t startTime;
  ofstream progressFile;
  ofstream progressFile2;
  if (MASSICCC == 1) {
    time(&startTime);
  }

	// Initialize output [HACK: using pointer because of an error when retrieving CriterionOutput later]
	//std::vector<CriterionName>* vCriterion = new std::vector<CriterionName>();
	std::vector<CriterionName>* vCriterion = new std::vector<CriterionName>();
	const int nbCriterion = _input->getCriterionName().size();
	for (int iCriterion = 0; iCriterion<nbCriterion; iCriterion++)
		vCriterion->push_back(_input->getCriterionName()[iCriterion]);
	_output = new ClusteringOutput(*vCriterion);

	// Main loop : build and run every model, and incrementally fill output
	// NOTE: potential parallelization here (OpenMP at least)


  const int nbModel = _input->getModelType().size();
  const int nbnbCluster = _input->getNbCluster().size();
  _output->clusteringModelOutputResize(nbModel * nbnbCluster);

#ifdef _OPENMP
#pragma omp parallel
#endif
{
  Data* workingData = _input->getDataDescription().getData()->clone();
  ClusteringStrategy* workingStrategy = _input->getStrategy()->clone();

	// Prepare reduced data if qualitative data and DATA_REDUCE==true
	std::vector<int64_t> correspondenceOriginDataToReduceData;
	if (_input->getDataType() == QualitativeData && DATA_REDUCE) {
		BinaryData * bData = dynamic_cast<BinaryData*> (_input->getDataDescription().getData());

		// initPartition
		Partition* inputInitPartition = NULL;
		Partition* workingInitPartition = NULL;
		if (_input->getStrategy()->getStrategyInit()->getStrategyInitName() == USER_PARTITION) {
			inputInitPartition = _input->getStrategy()->getStrategyInit()->getPartition(0);
		}

		try {
			// NOTE: workingInitPartition is always NULL
			workingData = bData->reduceData(
					correspondenceOriginDataToReduceData, workingKnownPartition,
					inputInitPartition, workingKnownPartition, workingInitPartition);
			workingStrategy = new ClusteringStrategy(*(_input->getStrategy()));
		}
		catch (Exception& errorType) {
			throw;
		}
	}


#ifdef _OPENMP
#pragma omp for collapse(2)
#endif
  for (int nbCluster_i = 0; nbCluster_i < nbnbCluster; nbCluster_i++) {
    for (int nbModel_i = 0; nbModel_i < nbModel; nbModel_i++) {
      int64_t nbCluster = _input->getNbCluster()[nbCluster_i];
      //Write progress in algo, not in this loop if both equals 1
      if (MASSICCC == 1 && nbnbCluster == 1 && nbModel == 1) {
        MASSICCC = 10;
      }
      ModelType* modelType = _input->getModelType()[nbModel_i];
      

			//-----------------
			// 1. Create model
			//-----------------

			if (VERBOSE == 1)
				std::cout << "Model name : "
					<< ModelNameToString(modelType->getModelName()) << std::endl;

			Model* model = NULL;
			switch (_input->getDataType()) {
			case QualitativeData:
				model = new BinaryModel(modelType, nbCluster, workingData,
						workingKnownPartition, correspondenceOriginDataToReduceData);
				break;
			case QuantitativeData:
				model = new Model(
						modelType, nbCluster, workingData, workingKnownPartition);
				break;
			case HeterogeneousData:
				model = new Model(
						modelType, nbCluster, workingData, workingKnownPartition);
				break;
			}

			//-----------------------------
			// 2. run parameters estimation
			//-----------------------------

			try {
				workingStrategy->run(model);
			}
			catch (Exception& e) {
				if (VERBOSE == 1) {
					Error error(e);
					error.run();
				}
				model->setError(e);
			}
			catch (...) {
				Exception * unknown = new OtherException(UnknownReason);
				model->setError(*unknown);
        delete unknown;
			}

			//--------------------
			// 3. compute criteria
			//--------------------

			ClusteringModelOutput* cmoutput;
			cmoutput = new ClusteringModelOutput(model);
			_output->addEstimation(cmoutput, nbModel_i + nbCluster_i * nbModel);
			if (model->getErrorType() == NOERROR) {
			  atLeastOneModelOk = true;
				// Loop over criterion name
				for (std::vector<CriterionName>::const_iterator it= vCriterion->begin(); it != vCriterion->end(); it++) {
					CriterionName criterion = *it;
					switch (criterion) {
					case BIC:
					{
						BICCriterion bic(model);
						bic.run(cmoutput->getCriterionOutput(BIC));
						break;
					}
					case CV:
						THROW(InputException, DAInput);
						break;
					case ICL:
					{
						ICLCriterion icl(model);
						icl.run(cmoutput->getCriterionOutput(ICL));
						break;
					}
					case NEC:
					{
						NECCriterion nec(model);
						nec.run(cmoutput->getCriterionOutput(NEC));
						break;
					}
					case UNKNOWN_CRITERION_NAME:
						THROW(OtherException, internalMixmodError);
						break;
					default:
						THROW(OtherException, internalMixmodError);
					}
				}
			}
			else {
				// Set criterion error (over all criterion names)
				for (std::vector<CriterionName>::const_iterator it= vCriterion->begin(); it != vCriterion->end(); it++) {
					CriterionName criterion = *it;
					cmoutput->setCriterionOutput(
						CriterionOutput(criterion, 0.0, model->getErrorType()));
				}
			}
      
      //Compute Entropy if needed (for massiccc's visualization)
      if (MASSICCC == 1 || MASSICCC == 10 || MASSICCC == 11) {
        vector< vector<double> > entropyM = model->getEntropyMatrix();
        int nbVariables = model->getData()->getPbDimension();
        ofstream entropyFile;
        entropyFile.open ("entropy" + std::to_string(nbModel_i + 1 + nbCluster_i * nbModel) + ".txt");
        for (int j = 0; j < nbVariables; j++) {
          for (int l = 0; l < nbCluster; l++) {
            entropyFile << entropyM[j][l] << "  ";
          }
          entropyFile << endl;
        } 
        entropyFile.close();
      }

			// model cannot be null (especially when we'll use enum class)
      delete model;
			//Write progress in file
      if (MASSICCC == 1) {
        time_t currTime;
        time(&currTime);
        double timePerModel = difftime(currTime, startTime) / ((double)count + 1);
        progressFile.open ("progress.json");
        progressFile << "{ \"Progress\" :  " << ((double)count + 1.0)/((double)nbModel * (double)nbnbCluster) * 100.0;
        progressFile << ", \"Estimated remaining time\" : " << timePerModel * (nbModel * nbnbCluster - count - 1) << " } ";
        progressFile.close();

      }
      count++;
    }
  }


	// Conditionally release memory
  if (_input->getKnownLabelDescription()) delete workingKnownPartition;
  if (_input->getDataType() == QualitativeData && DATA_REDUCE) {
    delete workingStrategy;
  }
  delete workingData;
//delete workingStrategy; //Invalid read dans le cas USER_PARTITION && example... 

//delete vCriterion;

  if (!atLeastOneModelOk){
    THROW(OtherException, AllModelsGotErros);
  }

} // End of openMP parallel section
}

}
