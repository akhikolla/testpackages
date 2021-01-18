/***************************************************************************
                             SRC/mixmod/Kernel/Parameter/GaussianGeneralParameter.h  description
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
#ifndef XEMGaussianGeneralParameter_H
#define XEMGaussianGeneralParameter_H

#include "mixmod/Kernel/Parameter/GaussianEDDAParameter.h"

namespace XEM {

// pre-declaration
class DiagMatrix;
class GeneralMatrix;

/**
  @brief Derived class of XEMGaussianParameter for Spherical Gaussian Model(s)
  @author F Langrognet
 */

class GaussianGeneralParameter : public GaussianEDDAParameter {

public:

	/// Default constructor
	GaussianGeneralParameter();

	/// Constructor
	// called by XEMModel
	GaussianGeneralParameter(Model * iModel, ModelType * iModelType);

	// Constructor
	// called by XEMStrategyType if initialization is USER
	GaussianGeneralParameter(
			int64_t iNbCluster, 
			int64_t iPbDimension, 
			ModelType * iModelType, 
			std::string & iFileName);
  // For heterogeneous models
	GaussianGeneralParameter(
			int64_t iNbCluster, 
			int64_t iPbDimension, 
			ModelType * iModelType, 
			std::string & iFileName,
      int64_t iNbVariable_binary,
      std::vector< int64_t > inbFactor);
	GaussianGeneralParameter(
			int64_t iNbCluster, 
			int64_t iPbDimension, 
			ModelType * iModelType, 
			double * proportions, 
			double ** means, 
			double *** variances);
	
	/// Constructor (copy)
	GaussianGeneralParameter(const GaussianGeneralParameter * iParameter);

	/// Destructor
	virtual ~GaussianGeneralParameter();

	/// reset to default values
	virtual void reset();

	/** @brief Selector
		@return A copy of the model
	 */
	Parameter * clone() const;

	void initUSER(Parameter * iParam);

	/// Compute table of sigmas of the samples of each cluster
	// NB : compute also lambda, shape, orientation, wk, w
	void computeTabSigma();

	/// Flury Algorithm
	/// return the value of Flury function
	double flury(double F);


	//     SELECTORS
	// ------ / -------- //
	double * getTabLambda() const;

	/** @brief Selector
		@return Table of shape matrix for each cluster
	 */
	DiagMatrix ** getTabShape() const;

	/** @brief Selector
		@return Table of orientation matrix for each cluster
	 */
	GeneralMatrix ** getTabOrientation() const;

	double getLogLikelihoodOne() const;

protected:
	
	/// Table of volume of each cluster
	double * _tabLambda; /* Volume      */

	/// Table of shape matrix of each cluster
	DiagMatrix ** _tabShape; /* Shape       */

	// Table of orientation matrix of each cluster
	GeneralMatrix ** _tabOrientation; /* Orientation */

	int64_t __storeDim;

	// model dependant methods for computing _tabSigma
	void computeTabSigma_L_C();
	void computeTabSigma_Lk_Ck();
	void computeTabSigma_L_Ck();
	void computeTabSigma_L_Dk_A_Dk();
	void computeTabSigma_Lk_Dk_A_Dk();
	void computeTabSigma_Lk_C();
	void computeTabSigma_L_D_Ak_D();
	void computeTabSigma_Lk_D_Ak_D();

	//void recopySymmetricMatrixInMatrix(SymmetricMatrix & sym, Matrix& mat, double facteur);
	int64_t getFreeParameter() const;
};

inline DiagMatrix ** GaussianGeneralParameter::getTabShape() const {
	return _tabShape;
}

inline GeneralMatrix ** GaussianGeneralParameter::getTabOrientation() const {
	return _tabOrientation;
}

inline double * GaussianGeneralParameter::getTabLambda() const {
	return _tabLambda;
}

}

#endif
