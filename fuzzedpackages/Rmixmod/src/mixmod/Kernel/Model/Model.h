/***************************************************************************
                             SRC/mixmod/Kernel/Model/Model.h  description
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
#ifndef XEMMODEL_H
#define XEMMODEL_H

#include "mixmod/Utilities/Util.h"
#include "mixmod/Utilities/Error.h"
#include "mixmod/Kernel/IO/Data.h"
#include "mixmod/Kernel/Parameter/Parameter.h"

namespace XEM {

/**
	  @brief Base class for Model(s)
	  @author F Langrognet
 */

// pre-declaration
class Parameter;
class GaussianData;
class BinaryData;
class ClusteringStrategyInit;
class Partition;
class Sample;
class ModelType;
class LabelDescription;

class Model {

public:

	/// Default constructor
	Model();

	//clone the model
	virtual Model * clone();
	/// Constructor
	Model(Model * iModel);

	/// Constructor
	Model(ModelType * modelType, int64_t nbCluster, Data *& data, Partition * knownPartition);

	/// Destructor
	virtual ~Model();

	void updateForCV(Model * originalModel, CVBlock & CVBlock);
	
	
	//-------
	// select
	//-------

	/** @brief Selector
		@return The parameters
	 */
	Parameter * getParameter();
	GaussianParameter * getGaussianParameter();
	BinaryParameter * getBinaryParameter();

	/** @brief Selector
		@return The current number for cluster
	 */
	int64_t getNbCluster();

	/** @brief Selector
		@return The current data
	 */
	Data * getData();

	/**
	 * @brief Return Gaussian data
	 * @return
	 */
	GaussianData* getGaussianData();
	/**
	 * @brief Return Binary data
	 * @return
	 */
	BinaryData* getBinaryData();

	/** @brief Selector
		@return The type of Error
	 */
	Exception& getErrorType() const;

	/** @brief Selector
	 @return The type of the model
	 */
	ModelType * const & getModelType() const;

	/** @brief Selector
		@return The number of samples
	 */
	int64_t getNbSample();

	/** @brief Selector
		@return Table of Fik of each cluster : probabilitites: _fik = pk * f(xi,muk,Sk)
	 */
	double ** getTabFik();

	/// return _tabSumF
	double * getTabSumF();

	/** @brief Selector
	    @return Table of Tik of each cluster :
	            conditional probabilities that xi arises from the k-th 
	            mixture component, 0 <= tik[i]k0] <= 1
	 */
	double ** getTabTik();

	/** @brief Selector
		@return Table of Zik zik[i][k0] = 1 if xi arises from the k0-th mixture component, 0 else
	 */
	int64_t ** getTabZikKnown();

	double ** getTabCik();

	/// getTabZikKnown
	bool * getTabZiKnown();

	/** @brief Selector
		@return Table of number of elements in each cluster
	 */
	double * getTabNk();

	bool getDeleteData();

	
	//---------
	// compute
	//--------

	/// compute _fik
	void computeFik();

	/// Compute the number of points in each class
	void computeNk();

	/** @brief Compute the log-likelihood
		@return The log-likelihood
	 */
	double getLogLikelihood(bool fikMustBeComputed);


	/** @brief Compute the log-likelihood with one cluster
		@return The log-likelihood
	 */
	double getLogLikelihoodOne();

	/** @brief Compute the entropy
		@return The entropy
	 */
	double getEntropy();

	/** @brief Compute the entropy matrix (for massiccc's visualization)
		@return The entropy matrix
	 */
	vector< vector<double> > getEntropyMatrix();

	/** @brief Compute the completed log-likelihood
		@return The completed log-likelihood
	 */
	double getCompletedLogLikelihood();

	/** get completed LL (if CEM) or LL (elseif)*/
	double getCompletedLogLikelihoodOrLogLikelihood();

	/// return the number of free parameters
	int64_t getFreeParameter();

	/** @brief Selector
		@return Log of the weight total
	 */
	double getLogN();

	/// getLabel and partition
	/// label=1...nbSample
	void getLabelAndPartitionByMAPOrKnownPartition(int64_t * label, int64_t ** partition);

	/// get label of the ith individual (i=0 .... nbSample-1) by MAP (or known label)
	/// return value in [0 nbCluster-1]
	int64_t getLabelByMAPOrKnownPartition(int64_t i);

	/// get knownLabel of the ith individual (i=0 .... nbSample-1)
	/// return value in [0 nbCluster-1]
	/// throw an error if the label is unknown
	int64_t getKnownLabel(int64_t i);

	/// getPostProba
	double ** getPostProba();


	//--------
	// compute
	//--------

	/** @brief Compute the label of the i0-th point of the sample
		@return The label of i0 (i0=0 -> _nBSample -1)
	 */
	int64_t computeLabel(int64_t i0);

	/** @brief Compute the label of new point x
		@return The label of x
	 */
	int64_t computeLabel(Sample * x);


	//------
	// algo
	//------

	/// Maximum a posteriori step method
	void MAPstep();

	/// Expectation step method
	void Estep();

	/// Maximization step method
	void Mstep();

	/// Stochastic classification step method
	void Sstep();

	/// Classification step method
	void Cstep();


	//-----
	// init
	//-----

	/// Random center initialization of the parameters of the model
	void initRANDOM(int64_t nbTry);

	/// random step for init RANDOM or USER_PARTITION
	void randomForInitRANDOMorUSER_PARTITION(
			bool * tabIndividualCanBeUsedForInitRandom, bool * tabClusterToInitialize);

	/// User initialization of the parameters of the model
	void initUSER(Parameter * initParameter);

	/// User partition initialization of the parameters of the model
	void initUSER_PARTITION(Partition * initPartition, int64_t nbTryInInit = defaultNbTryInInit);

	// set name of the algorithm
	void setParameter(Parameter * parameter);

	// set name of the algorithm
	void setAlgoName(AlgoName algoName);

	AlgoName getAlgoName();

	// set an error for the model
	void setError(Exception& errorType);

	/// Fix label Known
	void FixKnownPartition(Partition *& y);

	// edit debug information
	void editDebugInformation();
	void editFik();
	void editCik();
	void editTik();
	void editNk();

protected:

	/// type of the model
	ModelType * _modelType;

	/// Number of clusters
	int64_t _nbCluster;

	/// Number of samples
	int64_t _nbSample;

	/// Current data
	Data * _data;
	bool _deleteData;

	/// parameter of model
	Parameter * _parameter;

	/// Probabilities: _fik = pk * f(xi,muk,Sk)
	/// dim : _nbSample * _nbCluster
	double ** _tabFik;

	/// table of sum of _tabFik for all k (dim : _nbSample)
	double * _tabSumF;

	/// Conditional probabilities that x(i) arises from the k-th mixture component, 
	/// 0 <= tik[i][k0] <= 1. dim : _nbSample * _nbCluster
	double ** _tabTik;

	/// zikKnown : _tabZikKonwn[i][k0] = 1 if xi arises from the k0-th mixture component, 0 else.
	/// dim : _nbSample * _nbCluster
	int64_t ** _tabZikKnown;

	/** classification array for individual i and class k
	  // if zikKnown 
	  //		cik = zikKnown
	  //	else :
	  //		cik = tik if EM
	  //		cik = zik by MAP rule if CEM or MAP
	  //		cik = 'random' if SEM
	 */
	double ** _tabCik;


	/// is the label zik known (fixed)
	bool * _tabZiKnown;


	/// Number of points in each class
	double * _tabNk;

	// name of the algorithm
	AlgoName _algoName;

	// Error handler
	Error _error;
};

//--------------
//inline methods
//--------------

inline bool * Model::getTabZiKnown() {
	return _tabZiKnown;
}

inline int64_t ** Model::getTabZikKnown() {
	return _tabZikKnown;
}

inline double ** Model::getTabCik() {
	return _tabCik;
}

inline double ** Model::getTabTik() {
	return _tabTik;
}

inline double ** Model::getTabFik() {
	return _tabFik;
}

inline double * Model::getTabSumF() {
	return _tabSumF;
}

inline double * Model::getTabNk() {
	return _tabNk;
}

inline int64_t Model::getNbCluster() {
	return _nbCluster;
}

inline Data * Model::getData() {
	return _data;
}

inline GaussianData* Model::getGaussianData() {
	return _data->getGaussianData();
}

inline BinaryData* Model::getBinaryData() {
	return _data->getBinaryData();
}

inline Parameter * Model::getParameter() {
	return _parameter;
}

inline GaussianParameter * Model::getGaussianParameter() {
	return _parameter->getGaussianParameter();
}

inline BinaryParameter * Model::getBinaryParameter() {
	return _parameter->getBinaryParameter();
}

inline int64_t Model::getNbSample() {
	return _nbSample;
}

inline double ** Model::getPostProba() {
	return _tabTik;
}

inline ModelType * const & Model::getModelType() const {
	return _modelType;
}

inline Exception& Model::getErrorType() const {
	return _error.getError();
}

}

#endif
