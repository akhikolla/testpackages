/***************************************************************************
                             SRC/mixmod/Kernel/Algo/MAPAlgo.cpp  description
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
#include "mixmod/Kernel/Algo/MAPAlgo.h"
#include "mixmod/Kernel/Model/Model.h"

namespace XEM {

//------------
// Constructor
//------------
MAPAlgo::MAPAlgo() {
	_algoStopName = NBITERATION;
	_nbIteration = 1;
}

/// Copy constructor
MAPAlgo::MAPAlgo(const MAPAlgo & mapAlgo) : Algo(mapAlgo) {
}

//------------
// Constructor
//------------
MAPAlgo::MAPAlgo(AlgoStopName algoStopName, double epsilon, int64_t nbIteration)
: Algo(algoStopName, epsilon, nbIteration) {
}

//-----------
// Destructor
//-----------
MAPAlgo::~MAPAlgo() {
}

// clone
//------
Algo * MAPAlgo::clone() {
	return (new MAPAlgo(*this));
}

//---
//run
//---
void MAPAlgo::run(Model *& model) {
	_indexIteration = 0;
	//  model = model;
	model->setAlgoName(MAP);
	model->MAPstep(); // MAP Step
	//cout << "\nMAP algorithm \n";
}

}
