/***************************************************************************
                             SRC/mixmod/Kernel/Parameter/GaussianDiagParameter.h  description
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
#ifndef XEMGaussianDiagParameter_H
#define XEMGaussianDiagParameter_H

#include "mixmod/Kernel/Parameter/GaussianEDDAParameter.h"

namespace XEM {

// pre-declaration
class DiagMatrix;

/**
  @brief Derived class of XEMGaussianParameter for Spherical Gaussian Model(s)
  @author F Langrognet
 */

class GaussianDiagParameter : public GaussianEDDAParameter {

public:

	/// Default constructor
	GaussianDiagParameter();

	/// Constructor
	// called by XEMModel
	GaussianDiagParameter(Model * iModel, ModelType * iModelType);
    GaussianDiagParameter(int64_t iNbCluster, int64_t iPbDimension, ModelType * iModelType);
	/// Constructor (copy)
	GaussianDiagParameter(const GaussianDiagParameter * iParameter);

	/// Destructor
	virtual ~GaussianDiagParameter();

	/// reset to default values
	virtual void reset();

	/** @brief Selector
		@return A copy of the model
	 */
	Parameter * clone() const;

	/// initialisation USER
	void initUSER(Parameter * iParam);

	/// Compute table of sigmas of the samples of each cluster
	// NB : compute also lambda, shpae, orientation, wk, w
	void computeTabSigma();


	//     SELECTORS
	// ------ / -------- //

	double * getTabLambda() const;

	//double ** getTabShape();
	DiagMatrix** getTabShape() const;

	double getLogLikelihoodOne() const;

	int64_t getFreeParameter() const;

protected:
	
	/// Table of volume of each cluster
	double * _tabLambda; // Volume

	//double ** _tabShape;  
	DiagMatrix** _tabShape;
};

inline double * GaussianDiagParameter::getTabLambda() const {
	return _tabLambda;
}

inline DiagMatrix** GaussianDiagParameter::getTabShape() const {
	return _tabShape;
}

}

#endif
