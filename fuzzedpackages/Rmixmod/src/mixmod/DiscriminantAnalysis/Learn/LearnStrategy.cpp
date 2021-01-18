/***************************************************************************
                             SRC/mixmod/DiscriminantAnalysis/Learn/LearnStrategy.cpp  description
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

#include "mixmod/DiscriminantAnalysis/Learn/LearnStrategy.h"
#include "mixmod/Kernel/Algo/MAlgo.h"
#include "mixmod/Kernel/Algo/MAPAlgo.h"
#include "mixmod/Kernel/Model/Model.h"

namespace XEM {

//-----------
//Constructor
//-----------
LearnStrategy::LearnStrategy() {
}

//-----------
//Copy constructor
//-----------
LearnStrategy::LearnStrategy(const LearnStrategy & strategy) {
}

//----------
//Destructor
//----------
LearnStrategy::~LearnStrategy() {
}

//---
//run
//---
void LearnStrategy::run(Model * model) {
	// 1rst step of Discriminant analysis : USER_PARTITION with complete partition and M algorithm : only M step is done
	MAlgo MAlgo;
	MAlgo.run(model);
	// 2nd step : run the MAP Algo
	MAPAlgo MAPAlgo;
	MAPAlgo.run(model);
}

// verify method
bool LearnStrategy::verify() {
	return true;
}

}
