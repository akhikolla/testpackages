/***************************************************************************
                             SRC/mixmod/DiscriminantAnalysis/Learn/LearnModelOutput.h  description
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
#ifndef XEMLEARNMODELOUTPUT_H
#define XEMLEARNMODELOUTPUT_H

#include "mixmod/Kernel/IO/ModelOutput.h"

namespace XEM {

// pre-declaration
class Model;

/** 
 \class XEMLearnModelOutput
 @author F. Langrognet - R Lebret
		@date 2012
		@brief XEMLearnModelOutput class
 */
class LearnModelOutput : public ModelOutput {

public:

	/// Default Constructor
	LearnModelOutput();

	/// Initialization Constructor 1
	LearnModelOutput(Model * estimation);

	/// Initialization Constructor 2
	LearnModelOutput(ModelType & modelType, int64_t nbCluster, 
			std::vector<CriterionOutput*> & criterionOutput, double likelihood, 
			ParameterDescription & parameterDescription, LabelDescription & labelDescription,
			ProbaDescription & probaDescription);

	/// Initialization Constructor 3
	LearnModelOutput(ModelType & modelType, int64_t nbCluster, Exception& error);

	/// Copy Constructor
	LearnModelOutput(const LearnModelOutput & cModelOutput);


	/// Destructor
	virtual ~LearnModelOutput();


	/// set CV Labels
	void setCVLabel(Model * estimation, std::vector<int64_t> & cvLabel);

	/// get CV Label
	const LabelDescription * getCVLabel() const;

private:
	
	// label from CV criterion
	LabelDescription * _CVLabel;
};

inline const LabelDescription * LearnModelOutput::getCVLabel() const {
	return _CVLabel;
}

}

#endif
