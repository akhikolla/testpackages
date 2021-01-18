/***************************************************************************
                             SRC/mixmod/Kernel/Parameter/GaussianParameter.h  description
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
#ifndef XEMGaussianParameter_H
#define XEMGaussianParameter_H

#include "mixmod/Utilities/Util.h"
#include "mixmod/Kernel/Parameter/Parameter.h"

namespace XEM {

// pre-declaration
class Model;
class Matrix;
class DiagMatrix;

/**
  @brief Base class for XEMGaussianParameter(s)
  @author F. Langrognet
 */
class GaussianParameter : public Parameter {

public:

	//----------------------------
	// constructors / desctructors
	// ---------------------------

	/// Default constructor
	GaussianParameter();

	/// Constructor
	// called by XEMModel
	GaussianParameter(Model * iModel, ModelType * iModelType);

	/// Constructor
	// called by XEMGaussianEDDAParameter if initialization is USER
	GaussianParameter(int64_t iNbCluster, int64_t iPbDimension, ModelType * iModelType);

	/// Constructor
	GaussianParameter(const GaussianParameter * iParameter);

	/// Destructor
	virtual ~GaussianParameter();

	/// Comparison operator
	virtual bool operator ==(const GaussianParameter & param) const;

	virtual Parameter* clone()const = 0;

	/// reset to default values
	virtual void reset();
	
	
	//----------
	// selectors
	//----------

	/// get TabMean
	double ** getTabMean() const;

	
	//----------------
	// compute methods
	//----------------

	///computeDiagGlobalDataVariance
	void computeGlobalDiagDataVariance(DiagMatrix * matrixDiagDataVar);

	/// Compute table of cluster scattering matrices Wk and W
	virtual void computeTabWkW();

	/// compute label of idxSample
	int64_t computeClassAssigment(int64_t idxSample);

	/// Compute table of means of the samples for each cluster
	void computeTabMean();

	/// Compute table of means of the samples for each cluster in initUSER_PARTITION case
	/// outputs :
	/// -  nbInitializedCluster
	/// - tabNotInitializedCluster (array of size _nbCluster)
	void computeTabMeanInitUSER_PARTITION(int64_t & nbInitializedCluster, bool * tabNotInitializedCluster, Partition * initPartition);

	/// Compute table of sigmas of the samples for each cluster
	// NB : compute also lambda, shape, orientation, wk, w
	//virtual void computeTabSigma();

	//Compute normal probability density function
	//       for iSample the sample and kCluster th cluster
	virtual double getPdf(int64_t iSample, int64_t kCluster) const = 0;

	// compute normal probability density function
	// for all i=1,..,n and k=1,..,K
	virtual void getAllPdf(double ** tabFik, double * tabProportion) const = 0;

	// compute normal probability density function
	// for the line x within the kCluster cluster
	virtual double getPdf(Sample * x, int64_t kCluster) const = 0;


	//-----------
	// Algorithms
	//-----------

	/// Maximum a posteriori step method
	void MAPStep();

	/// Maximization step method
	virtual void MStep();


	// SELECTORS for square matrices

	/// get TabSigma
	//XEMMatrix ** getTabSigma();

	/** @brief Selector
		@return Table of cluster scattering matrices Wk
	 */
	Matrix ** getTabWk() const;

	/** @brief Selector
		@return Scattering matrix W
	 */
	Matrix * getW() const;

	// edit (for debug)
	virtual void edit() = 0;

	/// Edit
	virtual void edit(std::ofstream & oFile, bool text = false) = 0;

	virtual void input(std::ifstream & fi) = 0;

	virtual double getLogLikelihoodOne() const = 0;

	/// recopie sans faire construction / destruction
	// utilisé par SMALL_EM, CEM_INIT
	virtual void recopy(Parameter * otherParameter) = 0;

	virtual void updateForCV(Model * originalModel, CVBlock & CVBlock);

	//init
	//----
	void updateForInitRANDOMorUSER_PARTITION(Sample ** tabSampleForInit, bool * tabClusterToInitialize);

	/// init user
	virtual void initUSER(Parameter * iParam) = 0;

protected:

	// Square matrices
	// Table of covariance Matrix of each cluster
	//XEMMatrix ** _tabSigma;

	// Table of inverse of covariance matrix of each cluster
	//XEMMatrix ** _tabInvSigma; 

	/// Table of cluster scattering matrix Wk
	Matrix ** _tabWk;

	/// Scattering matrix W
	Matrix * _W;

	// 1/det(Sigma)
	//double * _tabInvSqrtDetSigma;

	/// Table of means vector of each cluster
	double ** _tabMean;

	// called by constructor
	// update _freeProportion
	void initFreeProportion(ModelType * iModelType);

	/// compute Mean when only one cluster
	/// called by initRANDOM, getLogLikelihoodOne
	void computeMeanOne(double * Mean, double * weight, double** y_Store, 
			int64_t nbSample, double totalWeight) const;

	void putIdentityInDiagonalMatrix(double * mat_store);

	void putIdentityInMatrix(double * mat_store);

	void initDiagonalMatrixToZero(double * A_store);

	double determinantDiag(double * mat_store, Exception& errorType);
};

//---------------
// inline methods
//---------------

inline double ** GaussianParameter::getTabMean() const {
	return this->_tabMean;
}

/*inline XEMMatrix ** XEMGaussianParameter::getTabSigma(){
  return _tabSigma;
}*/

inline Matrix ** GaussianParameter::getTabWk() const {
	return _tabWk;
}

inline Matrix * GaussianParameter::getW() const {
	return _W;
}

}

#endif
