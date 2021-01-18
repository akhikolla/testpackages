/***************************************************************************
                             SRC/mixmod/DiscriminantAnalysis/Learn/LearnOutput.h  description
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
#ifndef XEMLEARNOUTPUT_H
#define XEMLEARNOUTPUT_H

#include "mixmod/Utilities/Util.h"

namespace XEM {

// pre-declaration
class Criterion;
class Model;
class LearnMain;
class LearnModelOutput;

/** 
 \class XEMLearnOutput
 Main class for Learn Output (1rst step of discriminant analysis)
 @author F. Langrognet - R Lebret
		@date 2012
		@brief XEMLearnOutput class
 */
class LearnOutput {

public:

	/// Default Constructor
	LearnOutput();

	/// Copy Constructor
	LearnOutput(const LearnOutput & lOutput);

	/// Initialisation constructor
	LearnOutput( std::vector<Model*> const & estimations );

	/// Destructor
	virtual ~LearnOutput();

	/// Comparison operator
	bool operator ==(const LearnOutput & output) const;

	bool atLeastOneEstimationNoError() const;
	const int getNbEstimationWithNoError() const;

	/// sort the model output
	void sort(CriterionName criterionName);

	void editFile() const;

	/// return the index'th' ClusteringModelOutput
	/// Note : index is between 0 and size(ClusteringModelOutput)-1
	LearnModelOutput * getLearnModelOutput(unsigned int index) const;

	std::vector<LearnModelOutput*> const & getLearnModelOutput() const;

	int64_t getNbLearnModelOutput() const;

	void setLearnModelOutput(std::vector<LearnModelOutput *> & learnModelOutput);

private:
	
	// Vector containing output for each model
	std::vector<LearnModelOutput *> _learnModelOutput;
};

inline  std::vector<LearnModelOutput*> const & LearnOutput::getLearnModelOutput() const {
	return _learnModelOutput;
}

inline  LearnModelOutput *  LearnOutput::getLearnModelOutput(unsigned int index) const {
	if (index < _learnModelOutput.size()) {
		return _learnModelOutput[index];
	}
	else {
		THROW(InputException, wrongIndexInGetMethod);
	}
}

inline int64_t LearnOutput::getNbLearnModelOutput() const {
	return _learnModelOutput.size();
}

}

#endif
