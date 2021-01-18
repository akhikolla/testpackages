/***************************************************************************
                             SRC/mixmod/DiscriminantAnalysis/Predict/PredictOutput.h  description
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
#ifndef XEMPREDICTOUTPUT_H
#define XEMPREDICTOUTPUT_H

#include "mixmod/Utilities/Util.h"

namespace XEM {

// pre-declaration
class Model;
class PredictModelOutput;

/** 
 \class XEMPredictOutput
 Main class for Predict Output (2nd step of discriminant analysis)
 @author F. Langrognet - R Lebret
		@date 2012
		@brief XEMPredictOutput class
 */
class PredictOutput {

public:

	/// Default Constructor
	PredictOutput();

	/// Copy Constructor
	PredictOutput(const PredictOutput & lOutput);

	/// Initialisation constructor
	PredictOutput(Model * estimation);

	/// Destructor
	virtual ~PredictOutput();

	/// Comparison operator
	bool operator ==(const PredictOutput & output) const;

	/// return the index'th' ClusteringModelOutput
	/// Note : index is between 0 and size(ClusteringModelOutput)-1
	PredictModelOutput * getPredictModelOutput(unsigned int index) const;

	std::vector<PredictModelOutput *> const & getPredictModelOutput() const;
    void setPredictModelOutput(std::vector<PredictModelOutput *> & predictModelOutput);
private:
	
	// Vector containing output for each model
	std::vector<PredictModelOutput *> _predictModelOutput;
};

inline std::vector<PredictModelOutput *> const & PredictOutput::getPredictModelOutput() const {
	return _predictModelOutput;
}

inline PredictModelOutput * PredictOutput::getPredictModelOutput(unsigned int index) const {
	if ( index <= _predictModelOutput.size() ) {
		return _predictModelOutput[index];
	}
	else {
		THROW(InputException, wrongCriterionPositionInGet);
	}
}

}

#endif
