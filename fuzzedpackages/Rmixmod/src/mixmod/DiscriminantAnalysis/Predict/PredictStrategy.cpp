/***************************************************************************
                             SRC/mixmod/DiscriminantAnalysis/Predict/PredictStrategy.cpp  description
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

#include "mixmod/DiscriminantAnalysis/Predict/PredictStrategy.h"
#include "mixmod/Kernel/Algo/MAPAlgo.h"
#include "mixmod/Kernel/Parameter/Parameter.h"
#include "mixmod/Kernel/Model/Model.h"

namespace XEM {

//-----------
//Constructor
//-----------
PredictStrategy::PredictStrategy(Parameter * classificationRule)
: _classificationRule(classificationRule) {
}

//-----------
//Copy constructor
//-----------
PredictStrategy::PredictStrategy(const PredictStrategy & strategy)
: _classificationRule(strategy.getClassificationRule()) {
}

//----------
//Destructor
//----------
PredictStrategy::~PredictStrategy() {
}

//---
//run
//---
void PredictStrategy::run(Model * model) {
	// 1rst step of Discriminant analysis
	// use the initUSER() to add parameter form the Learn Step
	model->initUSER(_classificationRule);
	// only MAP step is done
	MAPAlgo algo;
	algo.run(model);
}

// verify method
bool PredictStrategy::verify() {
	return true;
}

}
