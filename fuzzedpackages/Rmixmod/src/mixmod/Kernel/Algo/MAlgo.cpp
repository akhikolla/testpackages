/***************************************************************************
                             SRC/mixmod/Kernel/Algo/MAlgo.cpp  description
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
#include "mixmod/Kernel/Algo/MAlgo.h"
#include "mixmod/Kernel/Model/Model.h"

namespace XEM {

//------------
// Constructor
//------------
MAlgo::MAlgo() {
	_algoStopName = NBITERATION;
	_nbIteration = 1;
}

/// Copy constructor
MAlgo::MAlgo(const MAlgo & mAlgo) : Algo(mAlgo) {
}

//------------
// Constructor
//------------
MAlgo::MAlgo(AlgoStopName algoStopName, double epsilon, int64_t nbIteration)
: Algo(algoStopName, epsilon, nbIteration) {
}

//-----------
// Destructor
//-----------
MAlgo::~MAlgo() {
}

// clone
//------
Algo * MAlgo::clone() {
	return (new MAlgo(*this));
}

//---
//run
//---
void MAlgo::run(Model *& model) {
	_indexIteration = 0;
	//  model = model;
	model->setAlgoName(M);
	model->Mstep(); // M Step
	model->Estep(); // E step to update Tik
	//cout << "\nMaximization algorithm \n";
}

}
