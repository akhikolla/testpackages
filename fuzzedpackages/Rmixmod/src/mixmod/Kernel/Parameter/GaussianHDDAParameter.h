/***************************************************************************
                             SRC/mixmod/Kernel/Parameter/GaussianHDDAParameter.h  description
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
#ifndef XEMGAUSSIANHDDAPARAMETER_H
#define XEMGAUSSIANHDDAPARAMETER_H

#include "mixmod/Kernel/Parameter/GaussianParameter.h"

namespace XEM {

// pre-declaration
class DiagMatrix;
class GeneralMatrix;
class SymmetricMatrix;

/**
  @brief Derived class of XEMGaussianParameter for HDDA Gaussian Model(s)
  @author F Langrognet
 */
class GaussianHDDAParameter : public GaussianParameter {

public:

	/// Default constructor
	GaussianHDDAParameter();

	/// Constructor
	// called by XEMModel
	GaussianHDDAParameter(Model * iModel, ModelType * iModelType);

	/// Constructor
	// called by XEMStrategyType
	GaussianHDDAParameter(int64_t iNbCluster, int64_t iPbDimension, 
			ModelType * iModelType, std::string & iFileName);

	/// Constructor
	GaussianHDDAParameter(const GaussianHDDAParameter * iParameter);

	/// Destructor
	virtual ~GaussianHDDAParameter();

	/// reset to default values
	virtual void reset();

	/** @brief Selector
		@return A copy of the model
	 */
	Parameter * clone() const;

	/** @brief Selector
		 @return Table of shape matrix for each cluster
	 */
	DiagMatrix ** getTabShape() const;

	/** @brief Selector
		@return Table of orientation matrix for each cluster
	 */
	GeneralMatrix ** getTabQ() const;

	/** @brief Selector
		@return Control the shape of the density in the subspace Ei
	 */
	double** getTabA() const;

	/** @brief Selector
		@return Control the shape of the density in the subspace orthogonal to Ei
	 */
	double* getTabB() const;

	/** @brief Selector
		@return Dimension of each subspace
	 */
	int64_t * getTabD() const;

	SymmetricMatrix ** getTabGammak() const;

	double ** getGamma() const;

	void MStep();

	
	// init
	//-----

	/// initialize attributes before an InitRandom  
	virtual void initForInitRANDOM();

	/// initialize attributes for init USER_PARTITION
	/// outputs :
	/// -  nbInitializedCluster
	/// - tabNotInitializedCluster (array of size _nbCluster)
	void initForInitUSER_PARTITION(int64_t & nbInitializedCluster, 
			bool * tabNotInitializedCluster, Partition * initPartition);

	/// User initialisation of the parameters of the model
	virtual void initUSER(Parameter* iParam);

	double getLogLikelihoodOne() const;

	void edit();

	void edit(std::ofstream & oFile, bool text = false);

	void recopy(Parameter * otherParameter);

	double getPdf(int64_t iSample, int64_t kCluster) const;

	void getAllPdf(double ** tabFik, double * tabProportion) const;

	double getPdf(Sample * x, int64_t kCluster)const;

	double* computeLoglikelihoodK(double** K);

	void input(std::ifstream & fi);

protected:
	
	/// Table of shape matrix of each cluster
	DiagMatrix** _tabShape;

	/// Table of orientation matrix of each cluster
	GeneralMatrix ** _tabQk;

	int64_t __storeDim;
	///Table of parameters Akj of each cluster
	double ** _tabAkj;
	///Table of parameters Bk of each cluster
	double * _tabBk;
	///Table of sub dimension of each cluster
	int64_t * _tabDk;

	/// _tabGammak = _Gammak * _Gammak' replaces matrix Wk when tabNk smaller than _pbDimension
	SymmetricMatrix ** _tabGammak; // matrice nk * p
	/// Array of individuals * weight - mean[k] in class k
	double ** _Gammak;

	int64_t getFreeParameter()const;

	void computeTabWkW();

	///compute function of cost for each tabQk_k
	double ** computeCost(GeneralMatrix ** tabQ)const;

	///compute parameters for the model AkjBkQk
	void computeAkjBkQk();
	///compute parameters for the model AkjBQk
	void computeAkjBQk();
	///compute parameters for the model AjBkQk
	void computeAjBkQk();
	///compute parameters for the model AjBQk
	void computeAjBQk();

	///compute parameters for the model AkBkQk
	void computeAkBkQk();
	///compute parameters for the model AkBQk
	void computeAkBQk();

	/// compute the intrinsic dimension when non given
	void computeTabDk();
};

inline DiagMatrix** GaussianHDDAParameter::getTabShape()const {
	return _tabShape;
}

inline GeneralMatrix ** GaussianHDDAParameter::getTabQ() const {
	return _tabQk;
}

inline double** GaussianHDDAParameter::getTabA() const {
	return _tabAkj;
}

inline double* GaussianHDDAParameter::getTabB() const {
	return _tabBk;
}

inline int64_t * GaussianHDDAParameter::getTabD() const {
	return _tabDk;
}

inline SymmetricMatrix** GaussianHDDAParameter::getTabGammak() const {
	return _tabGammak;
}

inline double ** GaussianHDDAParameter::getGamma() const {
	return _Gammak;
}

}

#endif
