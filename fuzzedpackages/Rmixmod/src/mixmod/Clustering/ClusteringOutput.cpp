/***************************************************************************
                             SRC/mixmod/Clustering/ClusteringOutput.cpp  description
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
#include "mixmod/Clustering/ClusteringOutput.h"
#include "mixmod/Clustering/ClusteringModelOutput.h"
#include "mixmod/Kernel/Model/Model.h"
#include "mixmod/Kernel/Criterion/BICCriterion.h"
#include "mixmod/Kernel/Criterion/ICLCriterion.h"
#include "mixmod/Kernel/Criterion/NECCriterion.h"
#include "mixmod/Kernel/IO/ProbaDescription.h"
#include "mixmod/Kernel/IO/Proba.h"
#include "mixmod/Kernel/IO/LabelDescription.h"
#include "mixmod/Kernel/IO/Label.h"
#include "mixmod/Kernel/Parameter/GaussianGeneralParameter.h"
#include "mixmod/Kernel/Parameter/BinaryParameter.h"
#include "mixmod/Kernel/Parameter/CompositeParameter.h"
#include "mixmod/Matrix/Matrix.h"

#include <algorithm>
#include <cmath>

namespace XEM {

//--------------------------------------
// Constructor from a vector of criteria
//--------------------------------------
ClusteringOutput::ClusteringOutput(std::vector<CriterionName> const & criterionName)
: _criterionName(criterionName) {
}

//-----------------
//  Copy constructor
//-----------------
ClusteringOutput::ClusteringOutput(const ClusteringOutput & cOutput)
: _clusteringModelOutput(cOutput.getClusteringModelOutput())
, _criterionName(cOutput.getCriterionName()) {
}

ClusteringOutput* ClusteringOutput::clone() {
	return new ClusteringOutput(*this);
}

// Add one ClusteringModelOutput (incremental construction)
void ClusteringOutput::addEstimation(ClusteringModelOutput* cmoutput, int64_t count) {
//_clusteringModelOutput.push_back(cmoutput); //Can't use push_back with openMP. Erase this line if there is no problems
  _clusteringModelOutput[count] = cmoutput; 
}

// Resize ClusteringModelOutput
void ClusteringOutput::clusteringModelOutputResize(int64_t size) {
  _clusteringModelOutput.resize(size);
}


//---------------------------
// Initialization Constructor (unused. TODO: remove)
//---------------------------
ClusteringOutput::ClusteringOutput(std::vector<Model*> const & estimations,
		std::vector<CriterionName> const & criterionName)
: _clusteringModelOutput(estimations.size())
, _criterionName(criterionName)
{
	// get a constant of the number of estimations
	const int64_t sizeEstimation = estimations.size();
	// loop over the estimations
	for (unsigned int i = 0; i < sizeEstimation; i++) {
		_clusteringModelOutput[i] = new ClusteringModelOutput(estimations[i]);
	}
}

//-----------
// Destructor
//-----------
ClusteringOutput::~ClusteringOutput() {
	for (unsigned int i = 0; i < _clusteringModelOutput.size(); i++) {
		if (_clusteringModelOutput[i] != NULL) {
			delete _clusteringModelOutput[i];
			_clusteringModelOutput[i] = NULL;
		}
	}
}

//---------------------
/// Comparison operator
//---------------------
bool ClusteringOutput::operator ==(const ClusteringOutput & cOutput) const {

  // NRtests must be made only with Eigen for full precision
  double precision = 0.0;
  if (XEMmathLib != 1) {
    cout << "****************************************************************" << endl;
    cout << "WARNING : NRtests must be launched with Eigen for full precision" << endl;
    cout << "Build will fail on Continuous Integration Server because of Eigen Unit Tests" << endl;
    cout << "Let XEMmathLib value to 1 in SelectLibrary.h before commiting anything" << endl;
    cout << "You can still check other libraries results/performances with the following precision" << endl;
    precision = 1.e-6; // Seems to be a good precision to check other libs
    cout << "Precision = " << precision << endl;
    cout << "****************************************************************" << endl;
  }
  cout.precision (std::numeric_limits<double>::digits10 + 1);
	int64_t nbSample = _clusteringModelOutput[0]->getProbaDescription()->getProba()->getNbSample();

	for (uint64_t k = 0; k < _clusteringModelOutput.size(); k++) {
		ClusteringModelOutput* cOutputThis = _clusteringModelOutput[k];
		ClusteringModelOutput* cOutputOther = cOutput._clusteringModelOutput[k];

		if ((cOutputThis->getStrategyRunError() == NOERROR &&
				cOutputOther->getStrategyRunError() != NOERROR) ||
			(cOutputThis->getStrategyRunError() != NOERROR &&
				cOutputOther->getStrategyRunError() == NOERROR))
		{
			cout << "UNEQUAL: one model failed, the other one did not" << endl;
			return false;
		}

		if (cOutputThis->getStrategyRunError() != NOERROR &&
				cOutputOther->getStrategyRunError() != NOERROR)
		{
			// Skip if both models have error
			continue;
		}

		// compare likelihood, probabilities, parameters, labels
		// NOTE: criterion values are never set
		//       (see XEMModelOutput::XEMModelOutput(XEMModel * estimation))
		//       ==> we don't check this - TODO later.

		//Checking likelihood
  	if (fabs(cOutputThis->getLikelihood() - cOutputOther->getLikelihood()) > precision) {
  		cout << "UNEQUAL: likelihood differ res" << k+1 << "  : computed " << cOutputThis->getLikelihood()
  				<< ", expected " << cOutputOther->getLikelihood() << endl;
  		return false;
  	}

		//Checking criterionValues
    int64_t j = 0;
    for (int64_t i = 0; i < maxNbCriterion; i++) {
      for (int64_t j = 0; j < maxNbCriterion; j++) {
        if (cOutputThis->getCriterionOutput(i).getCriterionName() == cOutputOther->getCriterionOutput(j).getCriterionName()) {
          if (fabs(cOutputThis->getCriterionOutput(i).getValue() - cOutputOther->getCriterionOutput(j).getValue()) > precision) {
            cout << "UNEQUAL: criterion differ res" << k+1 << "  : computed " << cOutputThis->getCriterionOutput(i).getValue()
              << ", expected " << cOutputOther->getCriterionOutput(j).getValue() << endl;
            return false;
          }
        }
      }
    }

		//precomputation: remap clusters
		//(may be needed if e.g. (111222333) and (222333111) are obtained)
		int64_t nbClusters = cOutputThis->getNbCluster();
		vector<int64_t> clustersCorrespondence(nbClusters);
		for (int64_t i = 0; i < nbClusters; i++)
			clustersCorrespondence[i] = -1;
		vector<int64_t> labelsThis = cOutputThis->getLabelDescription()->getLabel()->getLabel();
		vector<int64_t> labelsOther = cOutputOther->getLabelDescription()->getLabel()->getLabel();
		for (uint64_t i = 0; i < labelsThis.size(); i++) {
			if (clustersCorrespondence[labelsThis[i] - 1] < 0) {
				//WARNING: labels start at 1, not 0.
				clustersCorrespondence[labelsThis[i] - 1] = labelsOther[i] - 1;

				if (labelsThis[i] != labelsOther[i])
					cout << "WARNING: classes misalignment res" << k+1 << endl;
			}
		}


		// [TEMPORARY: THIS SHOULD NEVER HAPPEN]
		// sanity check: if some cluster correspondence is unassigned, we have an issue...
  	for (int64_t i = 0; i < nbClusters; i++) {
  		if (clustersCorrespondence[i] < 0) {
  			cout << "ERROR: component " << (i+1) << " is empty in res" << k+1 << endl;
  			return false;
  		}
  	}


		// NOTE [bauder]: the clean way from here would be to add a 'permutation' parameter
		// to every sub-class comparison operator. However, it would be quite intrusive since
		// there are many of them. So, for the moment, only essential checks are done 'by hand'.

		//Checking labels
  	for (int64_t i = 0; i < nbSample; i++) {
  		if (clustersCorrespondence[labelsThis[i] - 1] != labelsOther[i] - 1) {
  			cout << "UNEQUAL: labels differ resLabel" << k+1 << "   (at least) at row " << (i+1) << endl;
  			return false;
  		}
  	}

		//Checking probabilities
  	int64_t nbCluster = cOutputThis->getProbaDescription()->getProba()->getNbCluster();
  	vector<vector<double> > probaThis =
  			cOutputThis->getProbaDescription()->getProba()->getProba();
  	vector<vector<double> > probaOther =
  			cOutputOther->getProbaDescription()->getProba()->getProba();
  	for (int64_t i = 0; i < nbSample; i++) {
  		for (int64_t j = 0; j < nbCluster; j++) {
  			if (fabs(probaThis[i][j] - probaOther[i][clustersCorrespondence[j]]) > precision) {
  				cout << "UNEQUAL: probabilities differ resProba" << k+1<< "   (at least) at row " << (i+1)
  						<< ", cluster " << (j+1) << ": computed " << probaThis[i][j]
  						<< ", expected " << probaOther[i][clustersCorrespondence[j]] << endl;
  				return false;
  			}
  		}
  	}

    //Binary models
    if (isBinary(cOutputThis->getParameterDescription()->getModelType()->getModelName())) {
      BinaryParameter * bParameterThis = 
        dynamic_cast<BinaryParameter*> (cOutputThis->getParameterDescription()->getParameter());
      BinaryParameter * bParameterOther = 
        dynamic_cast<BinaryParameter*> (cOutputOther->getParameterDescription()->getParameter());
      int64_t pbDimension = bParameterThis->getPbDimension();
      for (int64_t j = 0; j < nbCluster; j++) {
        //Checking tabProportion
        if (fabs((bParameterThis->getTabProportion())[j] - (bParameterOther->getTabProportion())[clustersCorrespondence[j]]) > precision) {
          cout << "UNEQUAL: proportions differ resParam" << k+1 
            << ", cluster " << (j+1) << ": computed " << (bParameterThis->getTabProportion())[j]
            << ", expected " << (bParameterOther->getTabProportion())[clustersCorrespondence[j]] << endl;
            return false;
        }
        //Checking tabCenter
        for (int64_t h = 0; h < pbDimension; h++) {
          if (fabs((bParameterThis->getTabCenter())[j][h] - (bParameterOther->getTabCenter())[clustersCorrespondence[j]][h]) > precision) {
            cout << "UNEQUAL: means differ resParam" << k+1
              << ", cluster " << (j+1) << " Variable " << (h+1) << ": computed " << (bParameterThis->getTabCenter())[j][h]
              << ", expected " << (bParameterOther->getTabCenter())[clustersCorrespondence[j]][h] << endl;
              return false;
          }
        }
      }
      //Checking tabScatter
      double *** scatterThis = bParameterThis->scatterToArray();
      double *** scatterOther = bParameterOther->scatterToArray();
      for (int64_t j = 0; j < nbCluster; j++) {
        for (int64_t h = 0; h < pbDimension; h++) {
          int64_t tabNbModality = bParameterThis->getTabNbModality()[h];
          for (int64_t m = 0; m < tabNbModality; m++) {
            if (fabs(scatterThis[j][h][m] - scatterOther[clustersCorrespondence[j]][h][m]) > precision) {
              cout << "UNEQUAL: Scatter differ resParam" << k+1 
                << ", cluster " << (j+1) << " Scatter[" << h+1 << "][" << m+1 << "] " << ": computed " << scatterThis[j][h][m]
                << ", expected " << scatterOther[clustersCorrespondence[j]][h][m] << endl;
                return false;
            }
          }
        }
      }
      for (int64_t j = 0; j < nbCluster; j++) {
        for (int64_t h = 0; h < pbDimension; h++) {
          delete[] scatterThis[j][h];
          delete[] scatterOther[j][h];
        }
        delete[] scatterThis[j];
        delete[] scatterOther[j];
      }
      delete[] scatterThis;
      delete[] scatterOther;
      scatterThis = NULL;
      scatterOther = NULL;
    }

    //Gaussian models
    else if (isEDDA(cOutputThis->getParameterDescription()->getModelType()->getModelName())) {
      GaussianEDDAParameter * gParameterThis = 
        dynamic_cast<GaussianEDDAParameter*> (cOutputThis->getParameterDescription()->getParameter());
      GaussianEDDAParameter * gParameterOther = 
        dynamic_cast<GaussianEDDAParameter*> (cOutputOther->getParameterDescription()->getParameter());
      int64_t pbDimension = gParameterThis->getPbDimension();
      for (int64_t j = 0; j < nbCluster; j++) {
        //Checking tabProportion
        if (fabs((gParameterThis->getTabProportion())[j] - (gParameterOther->getTabProportion())[clustersCorrespondence[j]]) > precision) {
          cout << "UNEQUAL: proportions differ resParam" << k+1
            << ", cluster " << (j+1) << ": computed " << (gParameterThis->getTabProportion())[j]
            << ", expected " << (gParameterOther->getTabProportion())[clustersCorrespondence[j]] << endl;
            return false;
        }
        //Checking tabMean
        for (int64_t h = 0; h < pbDimension; h++) {
          if (fabs((gParameterThis->getTabMean())[j][h] - (gParameterOther->getTabMean())[clustersCorrespondence[j]][h]) > precision) {
            cout << "UNEQUAL: means differ resParam" << k+1
              << ", cluster " << (j+1) << " Variable " << (h+1) << ": computed " << (gParameterThis->getTabMean())[j][h]
              << ", expected " << (gParameterOther->getTabMean())[clustersCorrespondence[j]][h] << endl;
              return false;
          }
        }
        //Checking tabSigma
        double ** storeThis = gParameterThis->getTabSigma()[j]->storeToArray();
        double ** storeOther = gParameterOther->getTabSigma()[clustersCorrespondence[j]]->storeToArray();
        for (int64_t l = 0; l < pbDimension; l++) {
          for (int64_t m = 0; m < pbDimension; m++) {
            if (fabs(storeThis[l][m] - storeOther[l][m]) > precision) {
              cout << "UNEQUAL: Sigma differ resParam" << k+1
                << ", cluster " << (j+1) << " Sigma[" << l+1 << "][" << m+1 << "] " << ": computed " << storeThis[l][m]
                << ", expected " << storeOther[l][m] << endl;
                return false;
            }
          }
        }
        for (int64_t l = 0; l < pbDimension; l++) {
          delete[] storeThis[l];
          delete[] storeOther[l];
        }
        delete[] storeThis;
        delete[] storeOther;
        storeThis = NULL;
        storeOther = NULL;
      }
    }
    //Heterogeneous models
    else if (isHeterogeneous(cOutputThis->getParameterDescription()->getModelType()->getModelName())) { 
      CompositeParameter * cParameterThis = 
        dynamic_cast<CompositeParameter*> (cOutputThis->getParameterDescription()->getParameter());
      CompositeParameter * cParameterOther = 
        dynamic_cast<CompositeParameter*> (cOutputOther->getParameterDescription()->getParameter());
      //Binary Parameter
      BinaryParameter * bParameterThis = cParameterThis->getBinaryParameter(); 
      BinaryParameter * bParameterOther = cParameterOther->getBinaryParameter();
      int64_t pbDimension = bParameterThis->getPbDimension();
      for (int64_t j = 0; j < nbCluster; j++) {
        //Checking tabProportion
        if (fabs((bParameterThis->getTabProportion())[j] - (bParameterOther->getTabProportion())[clustersCorrespondence[j]]) > precision) {
          cout << "UNEQUAL: proportions differ resParam" << k+1 
            << ", cluster " << (j+1) << ": computed " << (bParameterThis->getTabProportion())[j]
            << ", expected " << (bParameterOther->getTabProportion())[clustersCorrespondence[j]] << endl;
            return false;
        }
        //Checking tabCenter
        for (int64_t h = 0; h < pbDimension; h++) {
          if (fabs((bParameterThis->getTabCenter())[j][h] - (bParameterOther->getTabCenter())[clustersCorrespondence[j]][h]) > precision) {
            cout << "UNEQUAL: means differ resParam" << k+1
              << ", cluster " << (j+1) << " Variable " << (h+1) << ": computed " << (bParameterThis->getTabCenter())[j][h]
              << ", expected " << (bParameterOther->getTabCenter())[clustersCorrespondence[j]][h] << endl;
              return false;
          }
        }
      }
      //Checking tabScatter
      double *** scatterThis = bParameterThis->scatterToArray();
      double *** scatterOther = bParameterOther->scatterToArray();
      for (int64_t j = 0; j < nbCluster; j++) {
        for (int64_t h = 0; h < pbDimension; h++) {
          int64_t tabNbModality = bParameterThis->getTabNbModality()[h];
          for (int64_t m = 0; m < tabNbModality; m++) {
            if (fabs(scatterThis[j][h][m] - scatterOther[clustersCorrespondence[j]][h][m]) > precision) {
              cout << "UNEQUAL: Scatter differ resParam" << k+1 
                << ", cluster " << (j+1) << " Scatter[" << h+1 << "][" << m+1 << "] " << ": computed " << scatterThis[j][h][m]
                << ", expected " << scatterOther[clustersCorrespondence[j]][h][m] << endl;
                return false;
            }
          }
        }
      }
      for (int64_t j = 0; j < nbCluster; j++) {
        for (int64_t h = 0; h < pbDimension; h++) {
          delete[] scatterThis[j][h];
          delete[] scatterOther[j][h];
        }
        delete[] scatterThis[j];
        delete[] scatterOther[j];
      }
      delete[] scatterThis;
      delete[] scatterOther;
      scatterThis = NULL;
      scatterOther = NULL;

      //Gaussian Parameter
      GaussianEDDAParameter * gParameterThis = 
        dynamic_cast<GaussianEDDAParameter*> (cParameterThis->getGaussianParameter());
      GaussianEDDAParameter * gParameterOther = 
        dynamic_cast<GaussianEDDAParameter*> (cParameterOther->getGaussianParameter());
      pbDimension = gParameterThis->getPbDimension();
      for (int64_t j = 0; j < nbCluster; j++) {
        //Checking tabProportion
        if (fabs((gParameterThis->getTabProportion())[j] - (gParameterOther->getTabProportion())[clustersCorrespondence[j]]) > precision) {
          cout << "UNEQUAL: proportions differ resParam" << k+1
            << ", cluster " << (j+1) << ": computed " << (gParameterThis->getTabProportion())[j]
            << ", expected " << (gParameterOther->getTabProportion())[clustersCorrespondence[j]] << endl;
            return false;
        }
        //Checking tabMean
        for (int64_t h = 0; h < pbDimension; h++) {
          if (fabs((gParameterThis->getTabMean())[j][h] - (gParameterOther->getTabMean())[clustersCorrespondence[j]][h]) > precision) {
            cout << "UNEQUAL: means differ resParam" << k+1
              << ", cluster " << (j+1) << " Variable " << (h+1) << ": computed " << (gParameterThis->getTabMean())[j][h]
              << ", expected " << (gParameterOther->getTabMean())[clustersCorrespondence[j]][h] << endl;
              return false;
          }
        }
        //Checking tabSigma
        double ** storeThis = gParameterThis->getTabSigma()[j]->storeToArray();
        double ** storeOther = gParameterOther->getTabSigma()[clustersCorrespondence[j]]->storeToArray();
        for (int64_t l = 0; l < pbDimension; l++) {
          for (int64_t m = 0; m < pbDimension; m++) {
            if (fabs(storeThis[l][m] - storeOther[l][m]) > precision) {
              cout << "UNEQUAL: Sigma differ resParam" << k+1
                << ", cluster " << (j+1) << " Sigma[" << l+1 << "][" << m+1 << "] " << ": computed " << storeThis[l][m]
                << ", expected " << storeOther[l][m] << endl;
                return false;
            }
          }
        }
        for (int64_t l = 0; l < pbDimension; l++) {
          delete[] storeThis[l];
          delete[] storeOther[l];
        }
        delete[] storeThis;
        delete[] storeOther;
        storeThis = NULL;
        storeOther = NULL;
      }
    } 
  }

  return true;
}

bool ClusteringOutput::atLeastOneEstimationNoError() const {
	bool res = false;
	unsigned int i = 0;
	while (res == false && i < _clusteringModelOutput.size()) {
		if ( (_clusteringModelOutput[i]->getStrategyRunError()) == NOERROR ) {
			res = true;
		}
		i++;
	}
	return res;
}

const int ClusteringOutput::getNbEstimationWithNoError() const {
	int n = getNbClusteringModelOutput();
	for ( unsigned int i = 0; i < _clusteringModelOutput.size(); i++ ) {
		if (!((_clusteringModelOutput[i]->getStrategyRunError()) == NOERROR )) {
			--n;
		}
	}
	return n;
}

//-----
// sort
//-----
void ClusteringOutput::sort(CriterionName criterionName) {
	std::sort(_clusteringModelOutput.begin(), _clusteringModelOutput.end(), SortByCriterion(criterionName));
}

void ClusteringOutput::editFile() const {
	//TODO
}

void ClusteringOutput::setClusteringModelOutput(std::vector<ClusteringModelOutput *> & clusteringModelOutput) {
	for (unsigned int i = 0; i < _clusteringModelOutput.size(); i++) {
		delete _clusteringModelOutput[i];
	}
	_clusteringModelOutput = clusteringModelOutput;
}

}
