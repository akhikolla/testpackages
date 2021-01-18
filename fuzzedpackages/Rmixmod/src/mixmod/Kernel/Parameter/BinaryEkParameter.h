/***************************************************************************
                             SRC/mixmod/Kernel/Parameter/BinaryEkParameter.h  description
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
#ifndef XEMBinaryEkParameter_H
#define XEMBinaryEkParameter_H

#include "mixmod/Kernel/Parameter/BinaryParameter.h"

namespace XEM {

/**
@author F LANGROGNET
 */
class BinaryEkParameter : public BinaryParameter {

public:
	
	/// Default constructor
	BinaryEkParameter();

	/// Constructor
	// called by XEMModel
	BinaryEkParameter(Model * iModel, ModelType * iModelType, int64_t * tabNbModality);

	/// Constructor
	BinaryEkParameter(const BinaryEkParameter * iParameter);

	/// Destructor
	~BinaryEkParameter();

	/// Comparison operator
	virtual bool operator ==(const BinaryEkParameter & param) const;

	/// reset to default values
	virtual void reset();

	/// clone
	Parameter * clone() const;

	/// selector :  return scatter value
	double * getScatter() const;

	/// getFreeParameter
	int64_t getFreeParameter() const;

	double getPdf(int64_t iSample, int64_t kCluster) const;

	long double getLogPdf(int64_t iSample, int64_t kCluster) const;

	/** Compute normal probability density function
		 for x vector and kCluster th cluster
	 */
	double getPdf(Sample * x, int64_t kCluster) const;

	/// getlogLikelihoodOne (one cluster)
	double getLogLikelihoodOne() const;

	/// Compute scatter(s) 
	void computeScatter();

	/// Compute random scatter(s)
	void computeRandomScatter();

	///recopy scatter from param (used for init  : USER)
	void recopyScatter(Parameter * iParam);

	///create Scatter from "Binary Parameter Ekjh"
	void createScatter(double *** scatter);

	/// editScatter (for debug)
	void editScatter(int64_t k);

	/// editScatter 
	void editScatter(std::ofstream & oFile, int64_t k, bool text = false);

	// Read Scatter in input file
	void inputScatter(std::ifstream & fi, int64_t k);
	void inputScatter(double *** scatters);

	double *** scatterToArray() const;

private:
	
	/// scatter
	double * _scatter;
};

inline double * BinaryEkParameter::getScatter() const {
	return _scatter;
}

}

#endif
