/***************************************************************************
                             SRC/mixmod/Kernel/Parameter/BinaryParameter.h  description
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
#ifndef XEMBinaryParameter_H
#define XEMBinaryParameter_H

#include "mixmod/Utilities/Util.h"
#include "mixmod/Kernel/Parameter/Parameter.h"

namespace XEM {

class GaussianParameter;
class BinaryParameter;

/**
  @brief Base class for XEMBinaryParameter(s)
  @author F. Langrognet
 */

class BinaryParameter : public Parameter {

public:

	//----------------------------
	// constructors / desctructors
	// ---------------------------

	/// Default constructor
	BinaryParameter();

	/// Constructor
	// called by XEMModel (via XEMBinary...Parameter)
	BinaryParameter(Model * iModel, ModelType * iModelType, int64_t * tabNbModality);

	// constructor
	// called if USER initialisation
	BinaryParameter(int64_t iNbCluster, int64_t iPbDimension, 
			ModelType * iModelType, int64_t * tabNbModality);

	/// Constructor
	BinaryParameter(const BinaryParameter * iParameter);

	/// Destructor
	virtual ~BinaryParameter();

	/// Comparison operator
	virtual bool operator ==(const BinaryParameter & param) const;


	/// reset to default values
	virtual void reset();

	/// create the same parameter than this but after updating because without xi0
	Parameter * createParameter(Model * iModel, int64_t i0, int64_t ki0);


	//----------
	// selectors
	//----------

	/// get TabCenter
	int64_t ** getTabCenter() const;

	/// get _tabNbModality
	int64_t * getTabNbModality() const;

	/// get total number of modality
	int64_t getTotalNbModality() const;

	//----------------
	// compute methods
	//----------------

	void getAllPdf(double ** tabFik, double * tabProportion) const;

	/** @brief Compute probability density
	@param iSample  Probability for sample iSample
	@param kCluster Probability in class kCluster
	 */

	virtual double getPdf(int64_t iSample, int64_t kCluster) const = 0;

	/** @brief Compute log probability density
	@param iSample  Probability for sample iSample
	@param kCluster Probability in class kCluster
	 */
	virtual long double getLogPdf(int64_t iSample, int64_t kCluster) const = 0;

	/** Compute normal probability density function
		 for x vector and kCluster th cluster
	 */
	//double getPdf(RowVector x, int64_t kCluster);
	virtual double getPdf(Sample * x, int64_t kCluster) const = 0;

	/// getlogLikelihoodOne (one cluster)
	virtual double getLogLikelihoodOne() const = 0;

	/// compute Tik for xi (i=0 -> _nbSample-1) when underflow
	virtual void computeTikUnderflow(int64_t i, double ** tabTik);

	/// Compute table of centers of the samples for each cluster
	void computeTabCenter();

	/// Compute scatter(s) 
	virtual void computeScatter() = 0;

	/// Compute random scatter(s)
	virtual void computeRandomScatter() = 0;

	///recopy scatter from param (used for init  : USER)
	virtual void recopyScatter(Parameter * iParam) = 0;


	//---------------
	// initialization
	//---------------

	/// init user
	void initUSER(Parameter * iParam);

	/// updateForInitRANDOMorUSER_PARTITION
	void updateForInitRANDOMorUSER_PARTITION(
			Sample ** tabSampleForInit, bool * tabClusterToInitialize);

	/// initialize attributes before an InitRandom 
	void initForInitRANDOM();

	/// initialize attributes for init USER_PARTITION
	/// outputs :
	/// -  nbInitializedCluster
	/// - tabNotInitializedCluster (array of size _nbCluster)
	void initForInitUSER_PARTITION(int64_t & nbInitializedCluster, 
			bool * tabNotInitializedCluster, Partition * initPartition);

	/// computeTabCenterInitUSER_PARTITIONoutputs :
	/// -  nbInitializedCluster
	/// - tabNotInitializedCluster (array of size _nbCluster)
	void computeTabCenterInitUSER_PARTITION(int64_t & nbInitializedCluster, 
			bool * tabNotInitializedCluster, Partition * initPartition);


	//-----------
	// Algorithms
	//-----------

	/// Maximum a posteriori step method
	void MAPStep();

	/// Expectation step method
	//void EStep();

	/// Maximization step method
	void MStep();


	//---------------
	// input / output
	//---------------

	// edit (for debug)
	void edit();

	/// editScatter (for debug)
	virtual void editScatter(int64_t k) = 0;

	/// Edit
	void edit(std::ofstream & oFile, bool text = false);

	/// editScatter 
	virtual void editScatter(std::ofstream & oFile, int64_t k, bool text = false) = 0;

	// Read Parameters in input file
	void input(std::ifstream & fi);

	// Read Parameters in input containers
	void input(
			double * proportions, 
			double ** centers, 
			double *** scatters);

	// Read Scatter in input file
	virtual void inputScatter(std::ifstream & fi, int64_t k) = 0;
	virtual void inputScatter(double *** scatters) = 0;

	/// recopie sans faire construction / destruction
	// utilise par SMALL_EM, CEM_INIT, SEM ...
	void recopy(Parameter * otherParameter);

	///create Scatter from "Binary Parameter Ekjh"
	virtual void createScatter(double *** scatter) = 0;
	virtual double *** scatterToArray() const = 0;

	void updateForCV(Model * originalModel, CVBlock & CVBlock);

protected:

	///compute TabCenter if there is only One Cluster
	// _tabCenter will not be changed
	void getTabCenterIfOneCluster(int64_t * tabCenter, double * tabNbSampleInMajorModality, 
			double ** tabNbSamplePerModality = NULL) const;

	/// Table of centers vector of each cluster
	int64_t ** _tabCenter;

	/// Table of modality
	int64_t * _tabNbModality;

	/// Total number of modality
	int64_t _totalNbModality;
};

/// compute Pdf in case nbCluster=1 (Scatter is a scalar)
double computePdfOneCluster(Sample * x, int64_t * Center, double Scatter, int64_t * tabNbModality);

/// compute Pdf in case nbCluster=1 (Scatter is a array of double, depends on variables)
double computePdfOneCluster(Sample * x, int64_t * Center, 
		double * Scatter, int64_t * tabNbModality);

/// compute Pdf in case nbCluster=1 
// (Scatter is a array of double*double, depends on variables and modalities)
double computePdfOneCluster(Sample * x, int64_t * Center, 
		double ** Scatter, int64_t * tabNbModality);

//---------------
// inline methods
//---------------

inline int64_t ** BinaryParameter::getTabCenter() const {
	return _tabCenter;
}

inline int64_t * BinaryParameter::getTabNbModality() const {
	return _tabNbModality;
}

inline int64_t BinaryParameter::getTotalNbModality() const {
	return _totalNbModality;
}

}

#endif
