/***************************************************************************
                             SRC/mixmod/Kernel/Parameter/GaussianSphericalParameter.h  description
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
#ifndef XEMGaussianSphericalParameter_H
#define XEMGaussianSphericalParameter_H

#include "mixmod/Kernel/Parameter/GaussianEDDAParameter.h"

namespace XEM {

/**
  @brief Derived class of XEMGaussianParameter for Spherical Gaussian Model(s)
  @author F Langrognet
 */

class GaussianSphericalParameter : public GaussianEDDAParameter {

public:

	/// Default constructor
	GaussianSphericalParameter();

	/// Constructor
	// called by XEMModel
	GaussianSphericalParameter(Model * iModel, ModelType * iModelType);
    GaussianSphericalParameter(int64_t iNbCluster, int64_t iPbDimension, ModelType * iModelType);
	/// Constructor (copy)
	GaussianSphericalParameter(const GaussianSphericalParameter * iParameter);

	/// Destructor
	virtual ~GaussianSphericalParameter();

	/** @brief Selector
		@return A copy of the model
	 */
	Parameter * clone() const;

	/// initialisation USER
	void initUSER(Parameter * iParam);

	/// Compute table of sigmas of the samples of each cluster
	// NB : compute also lambda, shpae, orientation, wk, w
	void computeTabSigma();

	double getLogLikelihoodOne() const;

	int64_t getFreeParameter() const;
};

}

#endif
